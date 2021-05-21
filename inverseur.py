#%%
import matplotlib.pyplot as plt
import numpy as np
import math as m

#pivot RCC a quatre cols
#parametres
sigma_adm = 500E6
h = 0.5
l = 10
E = 114E9
#calculs
theta_adm = 2*sigma_adm*l*180/(E*h*m.pi)
print(theta_adm)
