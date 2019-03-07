import sys 
import numpy as np
import matplotlib.pyplot as plt


#for i in range(1, len(sys.argv)):
#    mas = np.loadtxt(sys.argv[i], skiprows=2)
#    print "File " + sys.argv[1] + " loaded!"
#    plt.plot(mas[:,0], mas[:,1])
#plt.show()

files = ["./LOCAL DENSITY HISTOGRAMS (Iter = 1)", "./LOCAL DENSITY HISTOGRAMS (Iter = 2)"]

for f in files:
    mas = np.loadtxt(f, skiprows=0)
    print "File " + f + " loaded!"
    plt.plot(mas[1:-1,0], mas[1-1:,1])
plt.show()
