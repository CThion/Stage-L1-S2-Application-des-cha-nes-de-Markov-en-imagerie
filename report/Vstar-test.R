#========================================================================================================




imageGenerator=function(forme, N){
	#fonction permettant de générer une matrice-image de forme choisie, et de taille N*N
	img <- matrix(-1, nrow= N, ncol= N)	
	#----rectangle----
	if(forme== "rectangle"){
		img[(N/4):(N*3/4),(N/4):(N*3/4)] <- 1}
	#----losange----
	if(forme== "losange"){
		for (i in 1:N){for (j in 1:N){
				sigma0[i,j]=-1+2*(abs(i-j)<=N/4)*(abs(i+j-N)<=N/4)
	}}}#if&for&for
	#----damier----
	#----triangle----
	#----cercle----
	#----bande----
	if(forme== "bande"){
		for (j in 1:(N/2)){ img[,2*j]=rep(1,N)}
	}

	return(img)
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
	#la taille de l'image N doit être précisée (à remplacer par méthode info image)
	#L détermine la longueur maximale d'une branche
	#(iv, jv) appartient à [(0,1),(0,-1),(1,0),(-1,0)] ce qui donne le sens de la branche à partir de s

	v = c()
	vl = c()
	if(is<1 | is>N | js<1 | js>N){return("les coordonnées données sont en dehors de l'image")}

	#----récupération des L sites de la branche potentielles sans condition----
	for(k in(1:L)){
		i <- is+k*iv
		j <- js+k*jv
		if(i>=1 && i<=N){ if(j>=1 && j<=N){ #teste si les coordonnées du candidat sont dans l'image
			vl <- c(vl,img[i,j])
		}}#if&if
	}#for
	if(is.null(vl[1])){return(0)} #si on est sur le bord
	
	#----épuration selon la condition de valeur de site----
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


Vstar =function(img, N, i, j, L)a{
	s <- 0
	s <- s + vbranch(img, N, i, j, -1, 0, L)
	s <- s + vbranch(img, N, i, j, 1, 0, L)
	s <- s + vbranch(img, N, i, j, 0, -1, L)
	s <- s + vbranch(img, N, i, j, 0, 1, L)
	return(s)
}#func


#========================================================================================================

metropolisIsing=function(forme, N, beta, B, alpha, L, n, p){
	#fonction qui réalise l'algorithme de Metropolis sur une image img0 qui est bruité en img1, et doit à la fin être restaurée en img2
	
	#----création de l'image parfaite img0----
	img0 = imageGenerator(forme, N)
	image(img0)
	
	#----création d'une image auxiliaire img1, avec bord plus long et bruitage----
	img1 = matrix(0, nrow=N+2, ncol=N+2)
	bruit = matrix(2*rbinom(N^2,1,1-p)-1,nrow=N,ncol=N)
	img1[2:(N+1),2:(N+1)] = img0*bruit
	image(img1)
	img2 = img1 #image tamponz
	img3 = img1 #mémoire du bruit initiale

	for(k in 1:n){
	
		#----choix d'un site---
		i <- 1+ceiling(N*runif(1)) # indice entre 2 et N+1
		j <- 1+ceiling(N*runif(1))
		#i <- 2 + k%%N #de 2 à N+1
		#j <- 2 + k%/%(N-1)%%N #de 2 à N+1
	
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
		#if(dU<0){img1[i,j]= -img1[i,j]}
		#else{
		p <- exp(-U2)
		if(runif(1)<=p){img1[i,j]= -img1[i,j]}
		else{img2[i,j] <- img1[i,j]} #sinon on ne change pas
		#}#else
	}#for
	#----affichage----
	par(mfrow=c(1,3))
	image(1:N,1:N,img0)
	img3 <- img3[2:(N+1),2:(N+1)]
	image(1:N+2,1:N+2,img3)
	img1 <- img1[2:(N+1),2:(N+1)]
	image(1:N+2,1:N+2,img1)


	#image(img1)
	corresp <- correspondance(img0, img1, N)
	return(corresp) #utile pour MCdeMC_Vstar
	#calcul de la correspondance
}#func

#image(imageGenerator("bande",50))
#metropolisIsing(forme, N, beta, B, alpha, L, n, p)
metropolisIsing("bande", 64, 1, 0, 0.3, 5, 2*10^4, 0.15)

