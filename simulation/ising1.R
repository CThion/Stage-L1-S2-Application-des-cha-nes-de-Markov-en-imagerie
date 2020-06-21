library(stats)
library(ggplot2)
library(plot3D)
library('plot.matrix')

#====================VARIABLES====================

N = 100 #nombre ligne et colonnes

B = 2
beta = 15



globalU=function(B, beta, N){
  
  grid = matrix(2*rbinom(N*N, 1, 0.5)-1, nrow=N, ncol=N)  #génération d'une configuration aléatoire
  
  #====================Calcul de l'énergie====================
  
  energie_ij = matrix(0,nrow=N,ncol=N)
  
  #énergie pour les quatres coins
  
  energie_ij[1,1] = grid[1,1]*(grid[1,2]+grid[2,1])
  energie_ij[N,1] = grid[N,1]*(grid[N-1,1]+grid[N,2])
  energie_ij[1,N] = grid[1,N]*(grid[1,N-1]+grid[2,N])
  energie_ij[N,N] = grid[N,N]*(grid[N,N-1]+grid[N-1,N])
  
  #---- énergie pour ligne du haut et du bas (bordure) (3 voisins)
  for (i in 2:(N-1)) 
  {
    energie_ij[i,1]=grid[i,1]*(grid[i-1,1]+grid[i,2]+grid[i+1,1])
    energie_ij[i,N]=grid[i,N]*(grid[i-1,N]+grid[i,N-1]+grid[i+1,N])
  }
  
  #---- énergie pour colonnes gauche et droite (bordure) (3 voisins)
  for (j in 2:(N-1)) 
  {
    energie_ij[1,j]=grid[1,j]*(grid[1,j-1]+grid[2,j]+grid[1,j+1])
    energie_ij[N,j]=grid[N,j]*(grid[N,j-1]+grid[N-1,j]+grid[N,j+1])
  }
  
  #---- énergie pour le reste des sites (tout moins le contours)
  for (i in 2:(N-1)) 
  {
    for (j in 2:(N-1))
    {
      energie_ij[i,j]=grid[i,j]*(grid[i,j-1]+grid[i,j+1]+grid[i-1,j]+grid[i+1,j])
    }
    
  }
  
  energie_site=grid
  
  energie=-B*sum(energie_site)-beta*sum(energie_ij)
  
  
  
  #====================Monte-Carlo====================
  M=400 #nombre d'itération
  grid_mc = list()
  energie_mc = c()
  
  
  for (k in 1:M)  #fait tourner 400 fois le code
  {
    grid_mc[[k]]=matrix(2*rbinom(N*N,1,0.5)-1,nrow=N,ncol=N)
    energie_ij=matrix(0,nrow=N,ncol=N)
    
    #les quatres angles
    energie_ij[1,1]=grid_mc[[k]][1,1]*(grid_mc[[k]][1,2]+grid_mc[[k]][2,1])
    energie_ij[N,1]=grid_mc[[k]][N,1]*(grid_mc[[k]][N-1,1]+grid_mc[[k]][N,2])
    energie_ij[1,N]=grid_mc[[k]][1,N]*(grid_mc[[k]][1,N-1]+grid_mc[[k]][2,N])
    energie_ij[N,N]=grid_mc[[k]][N,N]*(grid_mc[[k]][N,N-1]+grid_mc[[k]][N-1,N])
    
    for (i in 2:(N-1)) #ligne du haut et du bas
    {
      energie_ij[i,1]=grid_mc[[k]][i,1]*(grid_mc[[k]][i-1,1]+grid_mc[[k]][i,2]+grid_mc[[k]][i+1,1])
      energie_ij[i,N]=grid_mc[[k]][i,N]*(grid_mc[[k]][i-1,N]+grid_mc[[k]][i,N-1]+grid_mc[[k]][i+1,N])
    }
    
    
    #colonnes de gauche et de droite
    
    for (j in 2:(N-1)) 
    {
      energie_ij[1,j]=grid_mc[[k]][1,j]*(grid_mc[[k]][1,j-1]+grid_mc[[k]][2,j]+grid_mc[[k]][1,j+1])
      energie_ij[N,j]=grid_mc[[k]][N,j]*(grid_mc[[k]][N,j-1]+grid_mc[[k]][N-1,j]+grid_mc[[k]][N,j+1])
    }
    
    
    #reste des sites (tous les sites moins les contours)
    
    for (i in 2:(N-1))
    {
      for (j in 2:(N-1))
      {
        energie_ij[i,j] = grid_mc[[k]][i,j]*(grid_mc[[k]][i,j-1]+grid_mc[[k]][i,j+1]+grid_mc[[k]][i-1,j]+grid_mc[[k]][i+1,j])
      }
      
    }
    energie_mc = c(energie_mc,-B*sum(grid_mc[[k]])-beta*sum(energie_ij))
  }
  
  
  #print(energie_mc)
  print(min(energie_mc))
  print(max(energie_mc))
  print(mean(energie_mc))
  minenergie=min(energie_mc)
  which.min(energie_mc)
  print(which.min(energie_mc))
  opti=min(which.min(energie_mc))
  print(grid_mc[[opti]])
  plot(grid_mc[[opti]])
  
}#function

globalU(B, beta, N)