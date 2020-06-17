#!/usr/bin/env python
# coding: utf-8

# # Simulation of metropolis and Gibbs algorithms

# https://www.markdownguide.org/basic-syntax

# ## Manage pictures with python 

# In[1]:


import random as rd  #https://www.w3schools.com/python/module_random.asp
from math import *

import PIL.Image as img  
#https://effbot.org/imagingbook/image.htm    
#https://pillow.readthedocs.io/en/4.2.x/reference/Image.html#PIL.Image.Image.load


# In[20]:


#https://raymond-namyst.emi.u-bordeaux.fr/ens/lycee/TD2.html
#https://calque.pagesperso-orange.fr/langages/python/imagemanipy.html
# differents modes https://pillow.readthedocs.io/en/latest/handbook/concepts.html#modes

def image_generator(length, name):
    """create a picture of length*length pixels, of two color.
    """
    size = (length,length)
    pic = img.new('1', size)
    pix = pic.load() #get pixel data
    for i in range(size[0]): #0 to length in rows
        for j in range(size[1]): #0 to 255 in cols
            pix[i,j] = (rd.choice([0,1]))
    pic.save(name+".png")
    #pic.show("picture.png")  #display the picture
    return


# In[24]:


#image_generator(8,"picture")


# In[43]:


def image_modificator(pic, i, j, x2):
    """
    """
    x2 = (x2+1)*255//2
    pic.putpixel((i,j), x2)


# In[3]:


def image_copy_modificator(pic, i, j, x2):
    """fonction to change the value x1 of pixel in row i and col j in an image (for an image in level of grey), by 
    the value x2. 
    """
    pic1 = img.open(pic)
    pic2 = pic1.copy() #create a copy of original image
    pic2.putpixel((i,j), x2) #makes the change x1 -> x2
    pic2.save("picture2.png") #save the copy


# In[25]:


#image_copy_modificator("picture.png",7,7,255)


# ## Compute the energy U and Us, for different markov random field (RMF)
# #### Ising RMF
# - <tex>$C_1:\ U(x_s) = -B x_s $</tex>  
# - <tex>$C_2:\ U_{c=(s,t)} = -\beta x_s x_t $</tex>
# <tex>$$ U(x) = -\sum_{c=(s,t)\in C}{\beta x_s x_t} - \sum_{s\in S}{B x_s} $$</tex>

# ### Global energy  <tex>$U$</tex>

# In[68]:


def ising_U_global(picture, beta, B):
    """Compute the global energy U of an image by totaling every local energy Us (one per pixel)
    """
    #VARIABLES
    pic = img.open(picture) #the picture, given by image generator
    size = pic.size #dimentions as tuple (width, height)
    U = 0
    #FONCTION
    for col in range(size[0]): 
        for row in range(size[1]):
            U += ising_U_local_tore(picture, beta, B, col, row)
            print("\n")
    return U


# ### Local energy  <tex>$U_s $</tex>
# p11 polymrf: l'énergie locale en un site est la somme des potentiels de toutes les cliques auquelles il appartient.

# In[75]:


def ising_U_local_tore(picture, beta, B, col, row):
    """
    """
    
        #VARIABLES
    pic = img.open(picture) #the picture, given by image generator
    size = pic.size #dimentions as tuple (width, height)
    
    assert((row,col)<size), "les coordonnées demandées sont en dehors du cadre des sites"
    
    xs = -1+pic.getpixel((row,col))*2/255 #energy of site s, with translation to (-1,1) value range
    Us = -B*xs #energy initialised with site's energy
    d, u, l, r = (col,row-1), (col,row+1), (col-1,row), (col+1,row) #down, up, left, right
    
        #FONCTION
    for v in [d,u,l,r]:
        print(v)
        xv = -1+pic.getpixel(v)*2/255 #if s is in the image's outile, not every neighbour would work 
        Us += -beta*xs*xv #add clique energy to Us 
        #print(f"xv={xv}") 
    #print(f"local energy Us in {(row,col)} = {Us}")
    return Us
            
    


# In[78]:


ising_U_local_tore("picture.png",1,1,7,7)


# In[67]:


ising_U_global("picture.png", 1, 1)


# In[74]:


imagex= img.open("picture.png")
imagex.getpixel((8,-1))


# In[13]:


def ising_U_local(picture, beta, B, col, row):
    """
    """
        #VARIABLES
    pic = img.open(picture) #the picture, given by image generator
    size = pic.size #dimentions as tuple (width, height)
    tst = 2//255
    xs = pic.getpixel((row,col))*tst-1 #energy of site s, with translation to (-1,1) value range
    Us = -B*xs #energy initialised with site's energy
    
    u, d, l, r = (col,row-1), (col,row+1), (col-1,row), (col+1,row) #down, up, left, right
    xu, xd, xl, xr, = pic.getpixel(u), pic.getpixel(d), pic.getpixel(l), pic.getpixel(r)
    u_rip, d_rip, l_rip, r_rip = (col, 0), (col, size[1]), (0, row), (size[0], row)
    V = [u,d,l,r] #list of neightbour coordinates
    xV = [xu, xd, xl, xr] #lsit of neighbour value
    RIP = [u_rip, d_rip, l_rip, r_rip] #list of try condition
    
    for (v, rip, xv) in zip(V, RIP, xV): #possibly four neightbours
        try: v>=rip #avoid u=(col,-1)
        except ExplicitException: del xV[xV.index(xv)] # then xv wont be used in Us compute
        try: v<=rip
        except ExplicitException: del xV[V.index(xv)] # then xv wont be used in Us compute
    for xv in xV:
        Us += -beta*xs*xv
    print(f"(xs, V, RIP, xV)={(xs, V, RIP, xV)}")
    return Us


# In[19]:


ising_U_local("picture.png",1,1,0,2)


# In[16]:


ising_U_local('picture.png', 1, 1, 0, 0)


# In[49]:


c.getpixel((0,0))


# In[37]:


c = img.open("picture.png")
print(c.getpixel((0,2)))
c.putpixel((0,2),1)
c.save("picture.png")
c.show("picture.png")


# In[ ]:





# ## Algorithms

# ### metropolis algorithm
# > 1. choix d'un site s  
# > 2. tirage aléatoire d'une descripteur  
# > 3. calcul de la variation d'énergie  
# > 4. décision

# In[7]:


def metropolis(picture1, beta, B, n):
    """fonction to compute the metropolis algorithm
    """
    #VARIABLES
    pic1 = img.open(picture1) #the picture, given by image generator
    width, height = pic1.size[0], pic1.size[1]  #dimentions pic.size return a tuple (width, height)
    coords = [] #list of coordinates
    for i in range(width): 
        for j in range(height): coords.append((i,j))
    dU = 0 #initialise energy gap
    not_change = 0 #stop condition

        #FONCTION
    for _ in range(n): #run n times     #while not_change <= 10 :

        #CHOIX D'UN SITE (aléatoire ici)
        s = rd.choice(coords) #(row, col) = (rd.randint(width), rd.randint(height))
        xs1 = pic1.getpixel(s) 
        
        #TIRAGE DESCRIPTEUR
        xs2 = rd.choice((0,255))
        
        #VARIATION Us
        image_copy_modificator(picture1, s[0], s[1], xs2) #create a copy "picture2" of image1 with value change of s 
        pic2 = img.open("picture2.png")
        U1 = ising_U_local(picture1, beta, B, s[0], s[1])
        U2 = ising_U_local("picture2.png", beta, B, s[0], s[1])
        dU = U2-U1
        
        #DECISION
        if dU < 0:  #global energy decreased
            pic1.putpixel(s, xs2) 
            not_change = 0
        else:  #global energy increased with p proba 
            continue
            p = exp(-beta*(dU))
            if rd.random() < p: 
                pic1.putpixel(s, xs2) 
                not_change = 0
            else:
                not_change += 1
        print(f"dU={dU}, U_global = {ising_U_global_sumLoc(picture1, beta, B)}")


# In[105]:


metropolis("picture1.png", 1, 1, 100)


# In[72]:


exp(13305)


# In[16]:


rd.choice([0,50])


# ### Gibbs algoritm
# > 1. choix d'un site s
# > 2. calcul de la probabilité conditionnelle locale <tex>$P(X_s = x_s | V_s)$</tex> 
# > 3. mise à jour du site par tirage aléatoire selon la loi <tex>$P(X_s = x_s | V_s)$</tex>
# 

# In[2]:


def gibbs():
    """recursive fonction to compute the gibbs' parser
    """
    return


# In[ ]:




