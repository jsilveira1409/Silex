#%%
import matplotlib.pyplot as plt
import numpy as np
import math as m

#variables
sigma_D = 500E6 # 
SF = 1.5  #facteur de securite
l = 25E-3   #longueur total de la lame
l_c = 5.94E-3  #longueur cols
E = 114E9   #young
b = 10.5E-3 #largeur des lames
contr_tech = 60

#calculs
#course admissible
c = (2*l_c)/l
sigma_adm = sigma_D/SF
h = l/contr_tech
f_adm = (c*(3 - 3*c +c**3)*contr_tech*sigma_adm/(3*E))


#rigidité
Ko = (2*b*(h**3)*E)/(c*(3-3*c + c**2)*l**3)     #rigidité naturelle
No = l*Ko


#plot
print("epaisseur des lames fines: ", '{:.6f}'.format(h) )
print("rigidité naturelle : ", '{:.6f}'.format(Ko) )
print("force necessaire: ", '{:.6f}'.format(No) )
print("course admissible ", '{:.6f}'.format(f_adm) )