#%%
import matplotlib.pyplot as plt
import numpy as np
import math as m

#variables
sigma_D = 500E6 # 
SF = 1.5  #facteur de securite
f_adm = 1.25E-3
E = 114E9   #young
contr_tech = 60
epsilon = 0.43#np.arange(0, 1, 0.01)

#calculs
sigma_adm = sigma_D/SF
l = (3*E*f_adm)/(epsilon*(3-3*epsilon + epsilon**2)*contr_tech*sigma_adm)
print(l)

#plot
axes = plt.gca()
axes.set_ylim([10E-3, 40E-3])
plt.plot(epsilon, l)
plt.xlabel("epsilon")
plt.ylabel("longueur lame[m]")
plt.grid()
plt.show()