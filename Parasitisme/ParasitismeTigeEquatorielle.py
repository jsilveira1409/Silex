#Calcul de la rotation parasite induit par la tige equatorielle
#une force positive ici veut dire que le miroir monte du cote des actionneurs
#d'ou le signe moins devant le deplacement selon z de la tige equatorielle

#%%
import matplotlib.pyplot as plt
import numpy as np
import math as m


# Constantes 

Fm_max = 2      #force moteur max
Dm = 12.7   #diametre miroir
a = 51      #distance entre miroir et application de la force, diametre du disque du support 
alpha = 1.5 #angle d'inclinaison du miroir en degree

# Variables

L = np.arange(1, 60, 0.01)

# Conversion

L = L/1000
Dm = Dm / 1000
a = a/1000
alpha = alpha*m.pi/180

# Calculs

Ft = m.sqrt(2)* Fm      #force sur la tige
zp = -m.sin(alpha)*(Dm/2 + a)   #deplacement de la tige equatorielle selon z 
deltaL = (3/5*L)*(m.sin(alpha)**2)*(Dm/2 + a)**2    #Raccourcissement de la tige dans le plan xy, indesire
Beta = 2*deltaL/(Dm + 2*a)      #rotation parasite selon z
Beta = Beta * 180/m.pi

# Plot


plt.plot (L, Beta)
plt.xlabel("Longueur de la tige equatorielle[m]")
plt.ylabel("Rotation parasite selon z [°]")
axes = plt.gca()
#axes.set_ylim([-0.7E-6, 0.7E-6])
#plt.figtext(.92, .8, "p = %sm"%(p) )
#plt.figtext(.92, .75, "l = %sm"%(l) )
#plt.figtext(.92, .7, "d = %sm"%(d) )
plt.savefig('line_plot.jpg')  
plt.show()
