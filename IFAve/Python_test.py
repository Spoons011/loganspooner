# Import the required libraries
import os
import numpy as np
import matplotlib.pyplot as plt
# Define the location of the directory
path = os.getcwd() +"/Data"
image_dir = os.getcwd()+"/Figures"
# Change the directory
os.chdir(path)


def make_graph(x,y,marker,color,fname):
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(x,y,'-',color=color,lw=2)
  fig.savefig(image_dir+'/'+fname,bbox_inches='tight',dpi=300)
  plt.close()

# Iterate over all the files in the directory
for file in os.listdir():
  if file.endswith('.txt'):
    x,y = np.genfromtxt(file, delimiter= '',skip_header=1,skip_footer=1, unpack=True)
    make_graph(x,y,'.','r',file.replace(".txt",""))
    
    
