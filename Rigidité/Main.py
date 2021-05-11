# _tp : tige polaire
# _m : miroir
# alpha : angle d'inclinaison du miroir
# _p : hauteur du centre de rotation virtuel des pivots rcc
# _mt : moteur, actionneur
# _eq : tige equatoriale
#_s : support
#_ts : tige superieur des contrepoids
#_cp : contrepoids
#_ti : tige inferieur des contrepoids
#%%
import matplotlib.pyplot as plt
import numpy as np
import math as m

#   Constantes
F_mt = 2      #force moteur
z_mt = 0      # course du moteur, position


#miroir
h_m = 6.4                    #hauteur du miroir en mm
D_m = 12.7   #diametre miroir
a = 51      #distance entre miroir et application de la force, diametre du disque du support 
yo_m = h_m/2     #axe neutre
M_m = 7


#support
h_s = 0
l_s = 0
A_s = 0
yo_s = h_m - h_s + h_s/2      #axe neutre
pmat_s = 0      #densite materiau
M_s = pmat_s *m.pi ((l_s+D_m/2)**2 - (D_m/2)**2)*h_s


#   Parametres

alphaVar = np.arange(-1.5, 1.5, 0.001)                      #en degrees, conversion apres, axe x
alphaFixe = 1.5                      #en degrees, conversion apres, axe x


# Tige equatoriale
l_eq = 10

# Contrepoids
M_cp = 0



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
ho_tp = 5.1                                    #hauteur du centre de rotation en mm
l_tp = 7                                   #longueur de la tige en mm
leff_tp = 0.699*l_tp
d_tp = 1.5                         #diametre de la tige en mm
E_tp = 196E9                        #module de young, en acier
v_tp = 0.3          # module de poisson


# Inverseur
E_inv = 190E9       # module de young
b_inv = 90          # largeur de la base
L_inv =  4          # longueur des lames
h_inv = 2           # epaisseur des lames
g_inv = 5           # distance entre le pivot et la tige
i_inv = 0           # distance entre le centre de masse et le centre de rotation pivot
Theta_inv = m.atan(z_mt/g_inv)       # angle dinclinaison du miroir -----> a definir
M_inv = 0           # masse 


#Tige superieur
E_ts = 190E9        # module de young
d_ts = 90           # diametre 
l_ts =  4           # longueur 
v_ts = 0.3          # module de poisson



#Tige superieur
E_ti = 190E9        # module de young
d_ti = 90           # diametre 
l_ti =  4           # longueur 
v_ti = 0.3          # module de poisson

#   Conversion

M_m = M_m/1000
alphaVar = m.pi*alphaVar/180
alphaFixe = m.pi*alphaFixe/180
ho_tp = ho_tp/1000
l_tp = l_tp/1000
d_tp = d_tp/1000
h_m = h_m/1000
D_m = D_m / 1000
l_eq = l_eq/1000
a = a/1000
b_tlp = b_tlp/1000 
h_tlp = h_tlp /1000
l_tlp =  l_tlp/1000
h_rcc = h_rcc/1000 
l_rcc =  l_rcc/1000         
p_rcc = p_rcc/1000        
b_inv = b_inv/1000         
L_inv =  L_inv/1000         
h_inv = h_inv/1000          
d_t = d_t/1000           
l_t =  l_t/1000           


# Calcul d'inertie
I_tp = (m.pi*d_tp**4)/64                 #inertie de la tige


############################################################################   TIGE POLAIRE

sigma_tp = E_tp*I_tp        #pour alleger les calculs
V_tp = (12*sigma_tp/l_tp**3)*((l_tp*np.tan(alphaVar)/2) + ho_tp*np.sin(alphaVar))   #Force tranchante sur la tige polaire
M_tp =  (-6*sigma_tp/l_tp**2 * (ho_tp*np.sin(alphaVar) + l_tp*np.tan(alphaVar)/3))       # moment de flexion sur la tige

#Racourcissement de la tige polaire
racTigePolaire = (1/(2 * sigma_tp**2) * (( V_tp**2 * l_tp**5)/20 + (V_tp * M_tp * l_tp**4)/4 + (M_tp**2 * l_tp**3)/3) ) 
racLevierPolaire = -ho_tp*(1 - np.cos(alphaVar))
deltaL_tp = -racTigePolaire - racLevierPolaire

#deplacement du point d'incidence du faiseau
deplFaisceauInc = (h_m - ho_tp)*(1/np.cos(alphaVar) - 1)

#Mouvement parasite total
Zparasite = deplFaisceauInc + deltaL_tp

#Plot du parasitisme tige polaire
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
plt.figtext(.3, .92, "p = %sm"%(ho_tp) )
plt.figtext(.5, .92, "l = %sm"%(l_tp) )
plt.figtext(.7, .92, "d = %sm"%(d_tp) )
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

K_tlp = (2*E_tlp*b_tlp*h_tlp**3)/l_tlp**3
K_rcc = (2*E_rcc*h_rcc**4*(l_rcc**2+3*p_rcc*l_rcc+3*p_rcc**2)/3*l_rcc**3)
K_t = (E_t*m.pi*d_t**4)/(128*l_t*(1+v_tp))
K_tp = (E_tp*m.pi*d_tp**4)/(128*l_tp*(1+v_tp))
K_inv = (2*E_inv*b_inv*h_inv**3)/3*L_inv

K_1 = K_rcc + K_t + K_tp
K_2 = 2* K_t
K_eq = 0
############################################################################  INERTIE DU SYSTEM
#Miroir et support
A_s = 2*h_s*l_s
A_m = h_m*D_m
yo = (yo_s*A_s + yo_m*A_m)/( A_m + A_s)  #axe neutre du miroir et support
Irot_s = pmat_s*m.pi*h_s*((a + D_m/2)**4 - (D_m/2)**4 )/4 + pmat_s*m.pi*h_s**3*((a + D_m/2)**2 - (D_m/2)**2 )/12
Irot_m = M_m*( (D_m/2)**2) /4 + (M_m*h_m**2) /12
Irot_ms = Irot_m + Irot_s + A_m*(yo_m - yo)**2 +A_s*(yo_s - yo)**2

#Inverseur
I_inv = 0   #a definir la forme materiau etc



############################################################################  MASSE REDUITE DU SYSTEM
#pour alleger les calcules
c1 = (Irot_ms + (M_m + M_s)*(ho_tp - yo)**2 ) /(( a +D_m/2)**2)    
c2 = M_cp
c3 =( I_inv + M_inv*i_inv**2 )/(g_inv**2)

M_red = c1 + c2 + c3

############################################################################  LOIS DES MOUVEMENTS
a_max = (F_mt - K_eq*z_mt)/(M_red)*1000     #acceleration maximale
z_max = 1.5 #course max de l'actionneur
theta_s = 0.015                                    #course angluaire pour le balayage de precision, en degree
y_s = theta_s*(m.pi/180)*(a + D_m/2)
T = m.sqrt(4*m.pi*y_s/a_max)
omega = 2*m.pi/T 
t_base = np.arange(0, T, 0.001)



############################################################################  PARASITISME

#inverseur lames croise non separe
par_inv = m.sqrt(2) * L_inv * Theta_inv /15

#lames des tables
par_tlp = (3*z_mt**2)/5*l_tlp        #ce sera vers le meme cote des lames, Ox, Oy

#tiges superieurs des cotrepoids 
par_ts = 3*par_tlp/(5*l_ts)

#tiges inferieurs des cotrepoids 
par_ti = 3*par_tlp/(5*l_ti)
