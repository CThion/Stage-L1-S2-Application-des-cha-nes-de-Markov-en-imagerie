#!/usr/bin/env python
# coding: utf-8

# # Simulation of metropolis and Gibbs algorithms

# ## Manage pictures with python 

# In[3]:


import random as rd
from math import *
import PIL.Image as img  #https://effbot.org/imagingbook/image.htm


# In[4]:


#https://raymond-namyst.emi.u-bordeaux.fr/ens/lycee/TD2.html
#https://calque.pagesperso-orange.fr/langages/python/imagemanipy.html
# differents modes https://pillow.readthedocs.io/en/latest/handbook/concepts.html#modes

def image_generator(length):
    """create a picture of length*length pixels, of two color.
    """
    size = (length,length)
    pic = img.new('1', size)
    pix = pic.load()
    for i in range(size[0]): #0 to length in rows
        for j in range(size[1]): #0 to 255 in cols
            pix[i,j] = (rd.choice([0,1]))
    pic.save("picture.png")
    pic.show("picture.png")  #display the picture
    return


# In[88]:


image_generator(128)


# In[5]:


def image_modificator(pic, i, j, x2):
    """fonction to change the value x1 of pixel in row i and col j in an image (for an image in level of grey), by 
    the value x2. 
    """
    x1 = pic.img.getpixel(i,j) #gets the value of the pixel row i col j
    pic.img.putpixel((i,j), x2) #makes the change x1 -> x2
    return


# ## Compute the energy U and Us, for different markov random field (RMF)
# #### Ising RMF
# - <tex>$C_1:\ U(x_s) = -B x_s $</tex>  
# - <tex>$C_2:\ U_{c=(s,t)} = -\beta x_s x_t $</tex>
# <tex>$$ U(x) = -\sum_{c=(s,t)\in C}{\beta x_s x_t} - \sum_{s\in S}{B x_s} $$</tex>

# ### Global energy  <tex>$U$</tex>

# In[49]:


def ising_U_global(beta, B):
    """compute the global energy U of an image acording to the ising MRF, with a neighbourhood of 4
    """
    #VARIABLES
    pic = img.open("picture.png") #the picture, given by image generator
    size = pic.size #dimentions as tuple (width, height)
    print(size)
    C1 = C2 = 0 #futur sum results
    
    #FONCTION
    for row in range(size[0]-1):#rows -- skip last row because of C2 neighbourhood compute
        for col in range(size[1]-1): #cols -- skip last column because of C2 neighbourhood compute
            print(f"(row,col)=({row},{col})")
            xs = pic.getpixel((row,col)) #get pixel value
            xb = pic.getpixel((row+1,col)) #bottom neighbour
            xr = pic.getpixel((row,col-1)) #right neighbour
            print(xs,xb,xr)
            C1 += -B*xs #add value of the pixel at (i,j), multiplied by -B (clique d'ordre 1)
            C2 += (-beta*xs*xb -beta*xs*xr)
    for col in range(size[0]): #add energy of last row's c1
        C1 += -B*pic.getpixel((size[0]-1, col))
        print("col",(size[0]-1,col),pic.getpixel((size[0]-1, col))) #U_s of last row 
    for row in range(size[0]-1): #add energy of last column's c1
        C1 += -B*pic.getpixel((row, size[0]-1))
        print("row",(row, size[0]-1),pic.getpixel((row, size[0]-1))) #U_s of last column
    global_U = C1+C2
    print(f"global_U = {global_U}\n C1 = {C1}\n C2 = {C2}")
    return  global_U


# In[50]:


def ising_U_global_sumLoc(beta, B):
    """
    """
    #VARIABLES
    pic = img.open("picture.png") #the picture, given by image generator
    size = pic.size #dimentions as tuple (width, height)
    U = 0
    #FONCTION
    i=0
    for row in range(size[0]): 
        for col in range(size[1]):
            U += ising_U_local2(beta, B, row, col)
            i+=1
    print(f"global energy U = {U}, i={i}")
    return U


# ### Local energy  <tex>$U_s $</tex>
# p11 polymrf: l'énergie locale en un site est la somme des potentiels de toutes les cliques auquelles il appartient.

# In[55]:


def ising_U_local2(beta, B, row, col):
    """
    """
        #VARIABLES
    pic = img.open("picture.png") #the picture, given by image generator
    size = pic.size #dimentions as tuple (width, height)
    xs = pic.getpixel((row,col)) #energy of site s
    Us = -B*xs #energy initialised with site's energy
    l, r, u, d = (row,col-1), (row,col+1), (row-1,col), (row+1,col) #left, right, up, down
    
        #FONCTION
    for v in [l,r,u,d]:
        try: xv = pic.getpixel(v) #if s is in the image's outile, not every neighbour would work
        except: pass 
        else: Us += -beta*xs*xv #add clique energy to Us
    #print(f"local energy Us in {(row,col)} = {Us}")
    return Us


# ## Algorithms

# ### metropolis algorithm
# > 1. choix d'un site s  
# > 2. tirage aléatoire d'une descripteur  
# > 3. calcul de la variation d'énergie  
# > 4. décision


def metropolis(n):
    """recursive fonction to compute the metropolis algorithm  
    """
    #VARIABLES
    pic = img.open("picture.png") #the picture, given by image generator
    width, height = pic.size[0], pic.size[1]  #dimentions pic.size return a tuple (width, height)
    dU = 0 #initialise energy gap

        #FONCTION
    for _ in range(n): #run n times

        #CHOIX D'UN SITE
        s = (rd.randint(width), rd.randint(height))
        #TIRAGE DESCRIPTEUR
        x2 = rd.choice((0,255))
        #VARIATION Us
        

    


# In[19]:


rd.choice([0,50])


# ### Gibbs algoritm
# > 1. choix d'un site s
# > 2. calcul de la probabilité conditionnelle locale <tex>$P(X_s = x_s | V_s)$</tex> 
# > 3. mise à jour du site par tirage aléatoire selon la loi <tex>$P(X_s = x_s | V_s)$</tex>
# 





