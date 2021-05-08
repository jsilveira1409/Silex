#%%
import matplotlib.pyplot as plt
import numpy as np
import math


L = np.arange(1E-3, 5E-3 , 1E-9)   #longueur de la tige

#Parametres
E = 60E9           #module de young
F = 2.5               #force du moteur
D = 57.35E-3           #distance entre CR et application de F
d = 3E-3          #diametre de la tige

alpha = 2.618E-2       #angle dinclinaison en rad

I = (math.pi*(d**4))/64                #inertie
mu = (F*D)/(E*I)


deltaL = (1/6)*((mu**2)*(L**3) -3*mu*(L**2)*alpha + 3*L*(alpha**2))

plt.plot(L, deltaL)
plt.xlabel("Longueur de la tige [m]")
plt.ylabel("parasitisme deltaL [m]")
plt.show()
