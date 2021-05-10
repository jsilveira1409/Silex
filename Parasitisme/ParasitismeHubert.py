#%%
import matplotlib.pyplot as plt
import numpy as np
import math as m

#   Parametres
alpha = np.arange(-1.5, 1.5, 0.001)                      #en degrees, conversion apres, axe x
p = 5.1                                    #hauteur du centre de rotation en mm
l = 7                                   #longueur de la tige en mm
d = 1.5                         #diametre de la tige en mm
E = 210E9                        #module de young

#   Constantes

hmiroir = 6.4                    #hauteur du miroir en mm

#   Conversion

alpha = m.pi*alpha/180
p = p/1000
l = l/1000
d = d/1000
hmiroir = hmiroir/1000

#   Calcul

I = (m.pi*d**4)/64                 #inertie de la tige
sigma = E*I
V = (12*sigma/l**3)*((l*np.tan(alpha)/2) + p*np.sin(alpha))   #Force tranchante
M =  (-6*sigma/l**2 * (p*np.sin(alpha) + l*np.tan(alpha)/3))

#racourcissement de la tige
racTige = (1/(2 * sigma**2) * (( V**2 * l**5)/20 + (V * M * l**4)/4 + (M**2 * l**3)/3) ) 
racLevier = -p*(1 - np.cos(alpha))
deltaL = -racTige - racLevier

# deplacement du point d'incidence du faiseau
deplFaisceauInc = (hmiroir - p)*(1/np.cos(alpha) - 1)

#Mouvement parasite total

Zparasite = deplFaisceauInc + deltaL


#   Plot
limiteSup = [0.5E-6]* alpha.size
limiteInf = [-0.5E-6]* alpha.size
alpha = np.arange(-1.5, 1.5, 0.001)                      #en degrees, conversion apres, axe x

plt.plot(alpha, deltaL, 'y')
plt.plot(alpha, deplFaisceauInc, 'g')
plt.plot(alpha, Zparasite, 'b')
plt.plot(alpha, limiteSup, 'r--')
plt.plot(alpha, limiteInf, 'r--')
plt.legend(['deltaL', 'depl. faisceau incident', 'depl. z parasite total'])
plt.xlabel("Angle d'inclinaison du miroir [°]")
plt.ylabel("parasitisme deltaL [m]")
axes = plt.gca()
axes.set_ylim([-0.7E-6, 0.7E-6])
plt.figtext(.3, .9, "p = %sm"%(p) )
plt.figtext(.5, .9, "l = %sm"%(l) )
plt.figtext(.7, .9, "d = %sm"%(d) )
plt.savefig('line_plot.svg')  
plt.show()
