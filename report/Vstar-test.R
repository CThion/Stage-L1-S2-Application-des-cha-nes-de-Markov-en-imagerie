#========================================================================================================
imageGenerator=function(forme, N){
  #fonction permettant de g�n�rer une matrice-image de forme choisie, et de taille N*N
  
  img <- matrix(-1, nrow= N, ncol= N)
  
  #----rectangle----
  if(forme== "rectangle"){
    img[(N/4):(N*3/4),(N/4):(N*3/4)] <- 1}
  
  #----losange----
  if(forme== "losange"){
    for (i in 1:N){
      for (j in 1:N){
        img[i,j]=-1+2*(abs(i-j)<=N/4)*(abs(i+j-N)<=N/4)  
      }
    }#for
  }#if
  
  #----damier----
  #----triangle----
  #----cercle----
  #----bande----
  if(forme=="bande"){
    for (j in 1:(N/2)){ img[,2*j]=rep(1,N) }
  }#if
  
  
  return(img)
}#func



#========================================================================================================
bruiteur=function(image, br, p){
  #fonction qui revoit l'image donn�e en argument par la m�me image mais bruit� par la m�thode "br"
  imgBr <- image
  N <- dim(image)[1]
  bruit <- 0
  p <- 0.15
  #----binomiale----
  if(br== "rbinom"){bruit <- matrix(2*rbinom(N^2,1,1-p)-1,nrow=N,ncol=N)}
  
  #----poisson----
  
  #----normale----
  
  imgBr <- image*bruit
  return(imgBr)
}#func



#========================================================================================================
correspondance=function(img0, img1, N){
	#----pourcentage de correspondance----
	nok <- 0 #nombre de pixels divergents
	correspondance <- 0 #pourcentage sur nb pixels total
	for (i in (1:N)){
		for (j in (1:N)){
			if(img0[i,j] != img1[i,j]){nok <- nok+1}
		}#for
	}#for
	correspondance <- ((N^2)-nok)/(N^2)*100
	return(correspondance)
}#func



#========================================================================================================
vbranch=function(img, N, is, js, iv, jv, L){
	#la taille de l'image N doit �tre pr�cis�e (� remplacer par m�thode info image)
	#L d�termine la longueur maximale d'une branche
	#(iv, jv) appartient � [(0,1),(0,-1),(1,0),(-1,0)] ce qui donne le sens de la branche � partir de s

	v = c()
	vl = c()
	if(is<1 | is>N | js<1 | js>N){return("les coordonn�es donn�es sont en dehors de l'image")}

	#----r�cup�ration des L sites de la branche potentielles sans condition----
	for(k in(1:L)){
		i <- is+k*iv
		j <- js+k*jv
		if(i>=1 && i<=N){ if(j>=1 && j<=N){ #teste si les coordonn�es du candidat sont dans l'image
			vl <- c(vl,img[i,j])
		}}#if&if
	}#for
	if(is.null(vl[1])){return(0)} #si on est sur le bord
	
	#----�puration selon la condition de valeur de site----
	lead <- vl[1]
	for(k in 1:length(vl)){ #parcours les L sites
		if(vl[k]==lead){ v <- c(v,vl[k]) }#condition pour faire partie de la branche.
	}#for
	
	#print(v)
	#print(length(v))
	return(length(v)*lead) #"poid" de la branche. length(v) fait office de facteur, comme beta classiquement
}#func
#vbranch(img, 50, 1,1, 0, -1, 50 )



#========================================================================================================
Vstar =function(img, N, i, j, L){
	s <- 0
	s <- s + vbranch(img, N, i, j, -1, 0, L)
	s <- s + vbranch(img, N, i, j, 1, 0, L)
	s <- s + vbranch(img, N, i, j, 0, -1, L)
	s <- s + vbranch(img, N, i, j, 0, 1, L)
	return(s)
}#func



#========================================================================================================
metropolisIsing=function(forme, N, beta, B, alpha, L, n, p){
	#fonction qui r�alise l'algorithme de Metropolis sur une image img0 qui est bruit� en img1, et doit � la fin �tre restaur�e en img2
	
	#----cr�ation de l'image parfaite img0----
	img0 = imageGenerator(forme, N)
	
	#----cr�ation d'une image auxiliaire img1, avec bord plus long et bruitage----
	img1 = matrix(0, nrow=N+2, ncol=N+2)
	img1[2:(N+1),2:(N+1)] = bruiteur(img0, "rbinom", p)
	img2 = img1 #image tampon
	img3 = img1 #m�moire du bruit initiale

	for(k in 1:n){
	
		#----choix d'un site---
		i <- 1+ceiling(N*runif(1)) # indice entre 2 et N+1
		j <- 1+ceiling(N*runif(1))
		#i <- 2 + k%%N #de 2 � N+1
		#j <- 2 + k%/%(N-1)%%N #de 2 � N+1
	
		#----tirage nouvelle valeure candidate----
		img2[i,j] = -img1[i,j]
	
		#----voisinage----
		#v = img1[i+1,j] + img1[i-1,j] + img1[i,j+1] + img1[i,j-1]
		v = Vstar(img1, N+2, i, j, L)
	
		#----potentiels----
		#U1 = -B*img1[i,j]-beta*img1[i,j]*v-alpha*img3[i,j]
		#U2 = -B*img2[i,j]-beta*img2[i,j]*v-alpha*img3[i,j]
		U1 = -B*img1[i,j]-img1[i,j]*(beta*v + alpha*img3[i,j])#potentiel actuel
		U2 = -B*img2[i,j]-img2[i,j]*(beta*v + alpha*img3[i,j])#potentiel candidat
	
		#----test----
		dU = U2-U1
		p <- exp(-U2)
		
		if(runif(1)<=p){img1[i,j]= -img1[i,j]}
		else{img2[i,j] <- img1[i,j]} #sinon on ne change pas
	}#for
	
	#----affichage----
	img3 <- img3[2:(N+1),2:(N+1)] #on retire les contours
	img1 <- img1[2:(N+1),2:(N+1)]
	
	par(mfrow=c(1,3))
	image(1:N,1:N,img0)#image initiale
	image(1:N+2,1:N+2,img3)#image bruit�e
	image(1:N+2,1:N+2,img1)#image restaur�e
  
	#----observation des performance de restauration----
	corresp <- correspondance(img0, img1, N)
	return(corresp) 
}#func


#metropolisIsing(forme, N, beta, B, alpha, L, n, p)
#metropolisIsing("bande", 64, 1, 0, 0.3, 5, 2*10^4, 0.15)



#========================================================================================================
recuitMetro=function(forme, N, beta, B, alpha, L, n, p){
  #fonction qui r�alise l'algorithme de Metropolis sur une image img0 qui est bruit� en img1, et doit � la fin �tre restaur�e en img2
  
  #----cr�ation de l'image parfaite img0----
  img0 = imageGenerator(forme, N)
  
  #----cr�ation d'une image auxiliaire img1, avec bord plus long et bruitage----
  img1 = matrix(0, nrow=N+2, ncol=N+2)
  img1[2:(N+1),2:(N+1)] = bruiteur(img0, "rbinom", p)
  img2 = img1 #image tampon
  img3 = img1 #m�moire du bruit initiale
  
  for(k in 1:n){
    
    #---recuit des parametres----
    t <- 1/log(k)
    #Br <- B/t
    #betar <- beta/t
    alphar <- alpha/t    
    
    #----choix d'un site---
    i <- 1+ceiling(N*runif(1)) # indice entre 2 et N+1
    j <- 1+ceiling(N*runif(1))
    #i <- 2 + k%%N #de 2 � N+1
    #j <- 2 + k%/%(N-1)%%N #de 2 � N+1
    
    #----tirage nouvelle valeure candidate----
    img2[i,j] = -img1[i,j]
    
    #----voisinage----
    #v = img1[i+1,j] + img1[i-1,j] + img1[i,j+1] + img1[i,j-1]
    v = Vstar(img1, N+2, i, j, L)
    
    #----potentiels----
    #U1 = -B*img1[i,j]-beta*img1[i,j]*v-alpha*img3[i,j]
    #U2 = -B*img2[i,j]-beta*img2[i,j]*v-alpha*img3[i,j]
    #U1 = -Br*img1[i,j]-img1[i,j]*(betar*v + alphar*img3[i,j])#potentiel actuel
    #U2 = -Br*img2[i,j]-img2[i,j]*(betar*v + alphar*img3[i,j])#potentiel candidat
    U1 = -B/t*img1[i,j]-img1[i,j]*(beta/t*v + alpha/t*img3[i,j])#potentiel actuel
    U2 = -B/t*img2[i,j]-img2[i,j]*(beta/t*v + alpha/t*img3[i,j])#potentiel candidat
    
    #----test----
    dU = U2-U1
    p <- exp(-U2)
    
    if(runif(1)<=p){img1[i,j]= -img1[i,j]}
    else{img2[i,j] <- img1[i,j]} #sinon on ne change pas
  }#for
  
  #----affichage----
  img3 <- img3[2:(N+1),2:(N+1)] #on retire les contours
  img1 <- img1[2:(N+1),2:(N+1)]
  
  par(mfrow=c(1,3))
  image(1:N,1:N,img0)#image initiale
  image(1:N+2,1:N+2,img3)#image bruit�e
  image(1:N+2,1:N+2,img1)#image restaur�e
  
  #----observation des performance de restauration----
  corresp <- correspondance(img0, img1, N)
  return(corresp) 
}#func

recuitMetro("rectangle", 64, 1, 0, 0.3, 10^5, 0.15)



#========================================================================================================
recuitMetroHV=function(forme, N, p, alpha, n){
  return(1)
}#func




#========================================================================================================
recuitSimuleGuyader = function(n, beta, N, p, alpha, stop){
  
  #----cr�ation de sigma0 le rectangle initial----
  sigma0 <- matrix(-1,nrow=N,ncol=N)  
  sigma0[(N/4):(N*3/4),(N/4):(N*3/4)] <- 1
  
  #----cr�ation de sigma le rectangle bruit�----
  bruit <- matrix(2*rbinom(N^2,1,1-p)-1,nrow=N,ncol=N) 
  sigma <- sigma0*bruit  
  
  #----cr�ation image auxiliaire, avec contour supl�mentaire pour calculer U plus facilement----
  Saux <- matrix(0, nrow=N+2, ncol=N+2) 
  Saux[2:(N+1),2:(N+1)] <- sigma
  sigmaaux <- Saux #sugmaaux = matrice auxiliaire pour le calcule de U, garde la valeur initiale de Saux (qui lui change)
  
  
  #--------le recuit simul�--------
  if(stop==0){ #cas o� on s'arr�te simplement apr�s n it�rations
    for (k in (1:n)) {
      #T <- 1/sqrt(sqrt((1+log(k))))
      T <- 1/k #T tend lentement vers 0
      iaux <- 1+ceiling(N*runif(1)) # indice entre 2 et N+1
      jaux <- 1+ceiling(N*runif(1))
      #iaux <- 2 + k%%N #de 2 � N+1
      #jaux <- 2 + k%/%(N-1)%%N #de 2 � N+1
      
      
      #----voisinage----
      #s <- Vstar(Saux, N, iaux, jaux)
      s <- Saux[iaux-1,jaux]+Saux[iaux+1,jaux]+Saux[iaux,jaux-1]+Saux[iaux,jaux+1]
      #print("s")
      #print(s)
      #----proba pour un recuit simul�, avec rapport de m�tropolis----
      r <- exp(-2*Saux[iaux,jaux]*(2*alpha*T*sigmaaux[iaux,jaux]+beta*s)/T)
      
      #----changement ou non de spin selon la proba r---
      if(runif(1)<r){Saux[iaux,jaux] <- -Saux[iaux,jaux]}
    }#for
  }#if
  
  if(stop==1){#cas ou on s'arr�te apr�s n it�rations sans changement de spin
    X<-0#nombre d'it�rations cons�cutives sans changement de spin
    it<-0
    while(X<=n) {
      it<-it+1
      #T <- 1/sqrt(sqrt((1+log(k))))
      T <- 1/k #T tend lentement vers 0
      iaux <- 1+ceiling(N*runif(1)) # indice entre 2 et N+1
      jaux <- 1+ceiling(N*runif(1))
      #iaux <- 2 + k%%N #de 2 � N+1
      #jaux <- 2 + k%/%(N-1)%%N #de 2 � N+1
      
      
      #----voisinage----
      s <- Saux[iaux-1,jaux]+Saux[iaux+1,jaux]+Saux[iaux,jaux-1]+Saux[iaux,jaux+1]
      
      #----proba pour un recuit simul�----
      r <- exp(-2*Saux[iaux,jaux]*(2*alpha*T*sigmaaux[iaux,jaux]+beta*s)/T)
      
      #----changement ou non de spin selon la proba r---
      if(runif(1)<r){
        Saux[iaux,jaux] <- -Saux[iaux,jaux]
        X <- 0#on vient de changer un spin, donc la condition d'arr�t retombe � 0
      }#if
      else{X <-X+1 } 
    }#while
    print(it)
  }#if
  
  S <- Saux[2:(N+1),2:(N+1)] #image finale (on retire le contour auxiliaire)
  
  
  #----pourcentage de correspondance----
  nok <- 0 #nombre de pixels divergents
  correspondance <- 0 #pourcentage sur nb pixels total
  for (i in (1:N)){
    for (j in (1:N)){
      if(S[i,j] != sigma0[i,j]){nok <- nok+1}
    }#for
  }#for
  correspondance <- ((N^2)-nok)/(N^2)*100
  #print(nok)
  #print(correspondance)
  
  
  #----affichage----
  par(mfrow=c(1,3))
  image(1:N,1:N,sigma0)
  image(1:N,1:N,sigma)
  image(1:N,1:N,S)
  return(correspondance)
}#function











