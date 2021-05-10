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
a = 51  #distance support point d'application       
Dm = 12.7   #diametre miroir
gamma = Dm/2 + a 

#######################################################   Parametres

M = 10          #contremasse/table a lames paralleles et bobine mobile

Mm = 10         #masse miroir
hm = 10         #hauteur miroir
Ms = 10         #masse support
hs = 10         #hauteur support
p = 6           #hauteur du pivot virtuel rcc
c = 5           #centre de masse miroir et support
Jinv = 10       #inertie inverseu
l_inv =10       #distance entre pivot de l'inverseur et tige attaché
M_inv =10       #masse inverseur
i = 10          #hauteur entre centre de masse de linverseur et centre de rotation

#######################################################   Calcul
Jms = Ms*a(a+Dm)/4 + (Ms*hs**2)/12 + Mm*((Dm**2)/2 + (hm**2)/3)/4         #inertie miroir et support ensemble
m_red = 2*M +Jms/(2*gamma**2) + 2*Mm*(p-c)**2 + Jinv/(l_inv**2) + 2*M_inv*i**2

#Plot
