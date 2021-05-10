# _tp : tige polaire
# _m : miroir
# alpha : angle d'inclinaison du miroir
# _p : hauteur du centre de rotation virtuel des pivots rcc
# _mt : moteur, actionneur
# _eq : tige equatoriale


#%%
import matplotlib.pyplot as plt
import numpy as np
import math as m

#   Constantes

h_m = 6.4                    #hauteur du miroir en mm
F_mt = 2      #force moteur
D_m = 12.7   #diametre miroir
a = 51      #distance entre miroir et application de la force, diametre du disque du support 

#   Parametres

alphaVar = np.arange(-1.5, 1.5, 0.001)                      #en degrees, conversion apres, axe x
alphaFixe = 1.5                      #en degrees, conversion apres, axe x
ho_tp = 5.1                                    #hauteur du centre de rotation en mm
l_tp = 7                                   #longueur de la tige en mm
d_tp = 1.5                         #diametre de la tige en mm
E_tp = 210E9                        #module de young

l_eq = 10

#   Conversion

alphaVar = m.pi*alphaVar/180
alphaFixe = m.pi*alphaFixe/180
ho_tp = ho_tp/1000
l_tp = l_tp/1000
d_tp = d_tp/1000
h_m = h_m/1000
l_eq = l_eq/1000
D_m = D_m / 1000
a = a/1000


# Calcul d'inertie

I_tp = (m.pi*d_tp**4)/64                 #inertie de la tige

############################################################################  PARASITISME TIGE POLAIRE
#   Calcul

sigma_tp = E_tp*I_tp        #pour alleger les calculs
V_tp = (12*sigma_tp/l_tp**3)*((l_tp*np.tan(alphaVar)/2) + ho_tp*np.sin(alphaVar))   #Force tranchante sur la tige polaire
M_tp =  (-6*sigma_tp/l_tp**2 * (ho_tp*np.sin(alphaVar) + l_tp*np.tan(alphaVar)/3))       # moment de flexion sur la tige

#racourcissement de la tige polaire
racTigePolaire = (1/(2 * sigma_tp**2) * (( V_tp**2 * l_tp**5)/20 + (V_tp * M_tp * l_tp**4)/4 + (M_tp**2 * l_tp**3)/3) ) 
racLevierPolaire = -ho_tp*(1 - np.cos(alphaVar))
deltaL_tp = -racTigePolaire - racLevierPolaire

# deplacement du point d'incidence du faiseau
deplFaisceauInc = (h_m - ho_tp)*(1/np.cos(alphaVar) - 1)

#Mouvement parasite total

Zparasite = deplFaisceauInc + deltaL_tp


#   Plot du parasitisme tige polaire
limiteSup = [0.5E-6]* alphaVar.size
limiteInf = [-0.5E-6]* alphaVar.size

plt.plot(alphaVar, deltaL_tp, 'y')
plt.plot(alphaVar, deplFaisceauInc, 'g')
plt.plot(alphaVar, Zparasite, 'b')
plt.plot(alphaVar, limiteSup, 'r--')
plt.plot(alphaVar, limiteInf, 'r--')
plt.legend(['deltaL', 'depl. faisceau incident', 'depl. z parasite total'])
plt.xlabel("Angle d'inclinaison du miroir [°]")
plt.ylabel("parasitisme deltaL [m]")
axes = plt.gca()
axes.set_ylim([-0.7E-6, 0.7E-6])
plt.figtext(.92, .8, "p = %sm"%(ho_tp) )
plt.figtext(.92, .75, "l = %sm"%(l_tp) )
plt.figtext(.92, .7, "d = %sm"%(d_tp) )
plt.savefig('line_plot.svg')  
plt.show()

############################################################################  PARASITISME DE LA TIGE EQUATORIALE


F_eq = m.sqrt(2)* F_mt      #force sur la tige
zp_eq = -m.sin(alphaFixe)*(D_m/2 + a)   #deplacement de la tige equatorielle selon z 
deltaL = (3/5*l_eq)*(m.sin(alphaFixe)**2)*(D_m/2 + a)**2    #Raccourcissement de la tige dans le plan xy, indesire
Beta = 2*deltaL/(D_m + 2*a)      #rotation parasite selon z
Beta = Beta * 180/m.pi
print("rotation parasite selon z a cause de la tige equatoriel : ", Beta)


############################################################################  RIGIDITE EQUIVALENTE
