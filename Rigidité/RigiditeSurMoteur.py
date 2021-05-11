# toutes les longueurs sont en mm, la conversion a m se fait automatiquement
# _tlp = table a lames paralleles
# _rcc = pivot rcc
# _tp = tige polaire
# _t = tige normal
# _inv = inverseur
# TOUT EST EN MM, la conversion se fait automatiquement

#%%
import matplotlib.pyplot as plt
import numpy as np
import math as m

#######################################################   Constantes

#######################################################   Parametres
# Constantes
a = 51  #distance support point d'application       
Dm = 12.7   #diametre miroir
gamma = Dm/2 + a 


    # Table a lames paralleles x2
    E_tlp = 190E9         #module de young
    b_tlp = 90          #largeur des lames
    h_tlp = 2           #epaisseur des lames
    l_tlp =  4          #longueur des lames

    # Pivot RCC
    E_rcc = 190E9         #module de young
    h_rcc = 90            # diametre de la tige, egale à e et à b dans la formule
    l_rcc =  4            #longueur des lames
    p_rcc = 6              #hauteur du pivot virtuel

    # Tige Polaire
    E_tp = 190E9        # module de young
    d_tp = 90           # diametre 
    l_tp =  4           # longueur 
    v_tp = 0.3          # module de poisson

    # Inverseur
    E_inv = 190E9       # module de young
    b_inv = 90          # largeur de la base
    L_inv =  4          # longueur des lames
    h_inv = 2           # epaisseur des lames
    l_inv = 5           # distance entre le pivot et la tige

    #Tige normal
    E_t = 190E9        # module de young
    d_t = 90           # diametre 
    l_t =  4           # longueur 
    v_t = 0.3          # module de poisson


# Conversion

b_tlp = b_tlp/1000 
h_tlp = h_tlp /1000
l_tlp =  l_tlp/1000
h_rcc = h_rcc/1000 
l_rcc =  l_rcc/1000         
p_rcc = p_rcc/1000         
d_tp = d_tp/1000          
l_tp =  l_tp/1000          
b_inv = b_inv/1000         
L_inv =  L_inv/1000         
h_inv = h_inv/1000          
d_t = d_t/1000           
l_t =  l_t/1000           

#######################################################   Conversion

#Differentes Rigidites Individuelles
K_tlp = (2*E_tlp*b_tlp*h_tlp**3)/l_tlp**3
K_rcc = (2*E_rcc*h_rcc**4*(l_rcc**2+3*p_rcc*l_rcc+3*p_rcc**2)/3*l_rcc**3)
K_t = (E_t*m.pi*d_t**4)/(128*l_t*(1+v_tp))
K_tp = (E_tp*m.pi*d_tp**4)/(128*l_tp*(1+v_tp))
K_inv = (2*E_inv*b_inv*h_inv**3)/3*L_inv

#rigidites equivalentes en paralleles
K_1 = K_rcc + K_t + K_tp
K_2 = 2* K_t

#Rigidite equivalente final

K_eq = 2*K_tlp + (K_2 + K_1 + K_t)/(gamma**2) + K_inv/(l_inv**2)

print(K_eq)


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