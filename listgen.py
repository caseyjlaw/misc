import os
import glob
import pickle
import matplotlib
import matplotlib.pyplot as plt 

def noiseplot(directory,y):  #User gives dir to pull files from and column to pull
    biglist=[]           # creates an empty list to be filled with files

    os.chdir(directory)
    for filename in glob.glob("*.pkl"):  
            biglist.append(filename)       #appends all .pkl files to biglist

    a = biglist              #shortens biglist to a just to be lazy
    
    image = []              #creates an empty list to append column data to

    for file in a:          #loops through file names in 'a' 
        f = open(file,'r')
        while True:         #for each file it loops through to get 'y' column data
            try:
                value=pickle.load(f)
                image.append(value[y]) #appends each certain column onto image list
            except EOFError:           #when hits EOFError it breaks
                f.close()
                break

    if y == 0:
        b='Integration Number'       #Label's are per email
       
    if y ==1: 
        b ='Estimate of Noise in Raw Data'
        
    if y == 2:
        b ='Fraction of Data Flagged'
   
    if y== 3:
        b ='Estimate of Noise in Single 5ms Image'
                               #image contains all files 'y' column data now
    
    plt.hist(image, bins=200)    #plots a histogram of the list made in while loop
    plt.xlabel(b)                               #assigns correct x-label based off user input
    plt.ylabel('Counts')
    return plt.show()            


    
