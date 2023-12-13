import numpy as np
import matplotlib.pyplot as plt
import os

path = os.getcwd() 
# Change the directory

image_dir = os.getcwd()+'/SESAPS'
os.chdir(path)


x,y1,y2 = np.genfromtxt('Logan_Bandpass_filter_S21.txt', delimiter = '	',usecols = (0,1,2), unpack = True)
y=20*np.log10(y1**2+y2**2)

fig=plt.figure()
ax=fig.add_subplot(111)
ax.set_xticks([0,0.35e9,0.7e9,1.1e9,1.42e9,1.8e9,2.1e9,2.5e9,2.8e9])
ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel('S11 db')
ax.plot(x,y)
fig.savefig(image_dir+'Plot',bbox_inches='tight',dpi=300)
plt.close()
plt.show()
