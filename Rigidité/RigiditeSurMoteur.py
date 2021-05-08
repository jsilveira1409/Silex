# toutes les longueurs sont en mm, la conversion a m se fait automatiquement
# _tlp = table a lames paralleles
# _rcc = pivot rcc
# _tp = tige polaire
# _t = tige normal
# _inv = inverseur
import matplotlib.pyplot as plt
import numpy as np
import math as m

#######################################################   Constantes

#######################################################   Parametres

# Table a lames paralleles x2
E_tlp = 190E9         #module de young
b_tlp = 90          #verifier
l_tlp =  4          #longueur des lames
h_tlp = 2           #

# Pivot RCC
E_rcc = 190E9         #module de young
b_rcc = 90            # diametre de la tige, egale à e et à b dans la formule
l_rcc =  4            #longueur des lames
p_rcc = 6              #hauteur du pivot virtuel

# Tige Polaire
E_tp = 190E9        # module de young
d_tp = 90           # verifier
l_tp =  4           # longueur des lames
v_tp = 0.3          # module de poisson
G_tp = E_tp/(2*(1+v_tp))

# Inverseur
E_inv = 190E9       # module de young
b_inv = 90          # largeur de la base
L_inv =  4          # longueur des lames
h_inv = 2           # epaisseur des lames



#######################################################   Conversion


# Differentes Inerties

I_tlp = (b_tlp*h_tlp**3)/12           
I_tp = m.pi*d_tp**4/64
I_rcc = b_rcc**4/12
I_inv = (b_inv*h_inv**3)/12

#Differentes Rigidites 
K_tlp = 24*E_tlp*I_tlp/(l_tlp**3)       
K_tp =  G_tp*I_tp/l_tp 
K_rcc = 8*E_rcc*I_rcc*(l_rcc**2 + e*p_rcc*l_rcc + 3*p_rcc**2)/(l_rcc**3)

K_eq = 2*K_tlp + K_tp + K_rcc



 #Plot
#plt.plot()
plt.legend(['deltaL', 'depl. faisceau incident', 'depl. z parasite total'])
plt.xlabel("Angle d'inclinaison du miroir [°]")
plt.ylabel("parasitisme deltaL [m]")
axes = plt.gca()
axes.set_ylim([-0.7E-6, 0.7E-6])

#plt.figtext(.92, .8, "p = %sm"%(p) )
#plt.figtext(.92, .75, "l = %sm"%(l) )
#plt.figtext(.92, .7, "d = %sm"%(d) )

plt.savefig('rigidite.svg')  
plt.show()


