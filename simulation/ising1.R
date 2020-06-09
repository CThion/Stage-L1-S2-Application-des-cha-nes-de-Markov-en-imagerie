library(stats)
library(ggplot2)
library(plot3D)
library('plot.matrix')
N=40
B=3
beta=1
grid=matrix(2*rbinom(N*N,1,0.5)-1,nrow=N,ncol=N)
grid
#print(grid)
#image2D(grid)
#plot(grid)

#for (i in 1:10)
#{grid=matrix(2*rbinom(N*N,1,0.1)-1,nrow=N,ncol=N)
#Sys.sleep(0.08)
#plot(grid)}

#Calcul de l'énergie

energie_ij=matrix(0,nrow=N,ncol=N)
energie_ij[1,1]=grid[1,1]*(grid[1,2]+grid[2,1])
energie_ij[N,1]=grid[N,1]*(grid[N-1,1]+grid[N,2])
energie_ij[1,N]=grid[1,N]*(grid[1,N-1]+grid[2,N])
energie_ij[N,N]=grid[N,N]*(grid[N,N-1]+grid[N-1,N])

for (i in 2:(N-1))
{
energie_ij[i,1]=grid[i,1]*(grid[i-1,1]+grid[i,2]+grid[i+1,1])
energie_ij[i,N]=grid[i,N]*(grid[i-1,N]+grid[i,N-1]+grid[i+1,N])
}

for (j in 2:(N-1))
{
  energie_ij[1,j]=grid[1,j]*(grid[1,j-1]+grid[2,j]+grid[1,j+1])
  energie_ij[N,j]=grid[N,j]*(grid[N,j-1]+grid[N-1,j]+grid[N,j+1])
}
for (i in 2:(N-1))
{
  for (j in 2:(N-1))
  {
    energie_ij[i,j]=grid[i,j]*(grid[i,j-1]+grid[i,j+1]+grid[i-1,j]+grid[i+1,j])
  }
  
}
#print(energie_ij)
energie_site=grid
energie=-B*sum(energie_site)-beta*sum(energie_ij)
#print(energie)


#Monte-Carlo

M=400
grid_mc=list()
energie_mc=c()
for (k in 1:M)
{
  grid_mc[[k]]=matrix(2*rbinom(N*N,1,0.5)-1,nrow=N,ncol=N)
  energie_ij=matrix(0,nrow=N,ncol=N)
  energie_ij[1,1]=grid_mc[[k]][1,1]*(grid_mc[[k]][1,2]+grid_mc[[k]][2,1])
  energie_ij[N,1]=grid_mc[[k]][N,1]*(grid_mc[[k]][N-1,1]+grid_mc[[k]][N,2])
  energie_ij[1,N]=grid_mc[[k]][1,N]*(grid_mc[[k]][1,N-1]+grid_mc[[k]][2,N])
  energie_ij[N,N]=grid_mc[[k]][N,N]*(grid_mc[[k]][N,N-1]+grid_mc[[k]][N-1,N])
  
  for (i in 2:(N-1))
  {
    energie_ij[i,1]=grid_mc[[k]][i,1]*(grid_mc[[k]][i-1,1]+grid_mc[[k]][i,2]+grid_mc[[k]][i+1,1])
    energie_ij[i,N]=grid_mc[[k]][i,N]*(grid_mc[[k]][i-1,N]+grid_mc[[k]][i,N-1]+grid_mc[[k]][i+1,N])
  }
  
  for (j in 2:(N-1))
  {
    energie_ij[1,j]=grid_mc[[k]][1,j]*(grid_mc[[k]][1,j-1]+grid_mc[[k]][2,j]+grid_mc[[k]][1,j+1])
    energie_ij[N,j]=grid_mc[[k]][N,j]*(grid_mc[[k]][N,j-1]+grid_mc[[k]][N-1,j]+grid_mc[[k]][N,j+1])
  }
  for (i in 2:(N-1))
  {
    for (j in 2:(N-1))
    {
      energie_ij[i,j]=grid_mc[[k]][i,j]*(grid_mc[[k]][i,j-1]+grid_mc[[k]][i,j+1]+grid_mc[[k]][i-1,j]+grid_mc[[k]][i+1,j])
    }
    
  }
  energie_mc=c(energie_mc,-B*sum(grid_mc[[k]])-beta*sum(energie_ij))
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
