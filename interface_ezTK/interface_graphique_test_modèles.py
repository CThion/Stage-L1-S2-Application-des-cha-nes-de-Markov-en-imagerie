from ezTK import *
from math import *
import settings as st

#==================================================================================================
def refresh():
    """fonction to update all parameters in st.dic (initialised in settings)
    """
    #update values
    st.dic = {
        "nbSite":int(st.nbSite.state),
        "nbState":int(st.nbState.state)
        }
    del st.displayFrame #remove current displayed sites' sqarred
    display() #create a new one with updated values

def makeStates():
    """fonction to return a tuple containing all states
    """
    return ("blue")

def control():
    """fonction to create a control band at the left. Will contain all the parameters controlers.
    """
    controlFrame = Frame(st.win, width=700, height=500, fold=2, bg="red")
    st.nbSite = Scale(controlFrame, scale=(4,100), command=lambda:(refresh())) #sqarre from 2*2 to 10*10 sites. 
    st.nbSite.set(64) #initialised with a 8*8
    st.nbState = Entry(controlFrame)
    st.nbState.insert(0,st.dic["nbState"])

#==================================================================================================
def display():
    """fonction to display sites ("pixels")
    """
    dispH = int(sqrt(st.dic["nbSite"])) #side lengh of the sites' squarre
    st.displayFrame = Frame(st.win, width=400, height=400, fold=dispH, grow=False)
    st.siteList = [] #stock all frames simuling sites
    for _ in range(1, st.dic["nbSite"]+1): #make as many frame as specified in the scale st.nbSite
        st.siteList.append(Frame(st.displayFrame, width=400/dispH, height=400/dispH, bg=makeStates(), bd=2))
        
#==================================================================================================
def evolution():
    """fonction to follow the evolution of simulation wile it's running
    """
    evolutionFrame = Frame(st.win)
    
#==================================================================================================
def main():
    """main fonction of the programm, firstly called
    """
    st.win = Win(title="zone de test", flow='ES', fold=3)
    display()
    control()
    evolution()
    
    #----
    st.win.loop()




if __name__ == "__main__":
    main()
