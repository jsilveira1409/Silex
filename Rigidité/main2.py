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
va te pendre pelo
###############################################PROGRAM OPTIONS
graphs = True                   #print graphs ?

#Constantes
F_mt = 1.8      #force moteur
z_mt = 2.5      #course du moteur, position
z_max = 1       #course max du moteur, depuis la position nominal
alpha_max = 1.5 #inclinaison max
M_mt = 7

#Parametres

alphaVar = np.arange(-1.5, 1.5, 0.001)        #en degrees, conversion apres, axe x
alphaFixe = 1.5                               #en degrees, conversion apres, axe x

#miroir
D_m = 12.7                      #diametre miroir
h_m = 6.4                       #hauteur du miroir en mm
M_m = 2                         #masse du miroir
yo_m = h_m/2                    #axe neutre
a = (180*z_max/(alpha_max*m.pi)) - D_m/2                          #distance entre miroir et application de la force   
R_m = a + D_m/2                 #distance entre centre du miroir et application des rcc

#support
h_s = 2                         #epaisseur
d_s =  2*R_m                        #largeur    
A_s = 0.01                        #surface        ??
yo_s = h_m - h_s + h_s/2        #axe neutre
pmat_s = 0.01                     #densite materiau

# Tige Polaire
d_tp = 1                     #diametre de la tige en mm
l_tp = 11.5                        #longueur de la tige en mm
ho_tp = 1.5                    #hauteur du centre de rotation en mm, distance entre point de contacte de la tige et miroir et le centre de rotation
leff_tp = 0.699*l_tp            #longueur efficace, Pinned-fixed
E_tp = 196E9                    #module de young, en acier
v_tp = 0.3                      #module de poisson

# Masse inertielle
M_cp = 7                        #masse

# Table a lames paralleles x2
contr_adm_tlp = 800E6           #contrainte admissible, acier bohler 
b_tlp = 10                      #largeur des lames
e_tlp = 15                       #distance entre les lames
l_tlp = 30                      #longueur des lames
E_tlp = 196E9                       #module de young
h_tlp = 1                       #epaisseur des lames
N_tlp = 200                     #charge du ressort sur les tlb
N_tlp_c = 1             #charge critique, calculé apres
N_tlp_o = 1             #parametre pour rigidite avec charge, calculé apres
tlpCharge = True                #sont-elles chargées?

# Tige equatoriale
l_eq = 30                      #longueur
d_eq = 1                       #diametre
E_eq = 196E9                   #module de young

# Pivot RCC
d_rcc = 0.5                    #diametre de la tige, egale à e, à b et à h dans la formule
ht_rcc = 5                     #approximation des tiges rcc qui se trouvent sur le moteur actionné, projection sur z
l_rcc =  m.sqrt(2)*ht_rcc      #longueur des lames
h_s = h_m - ho_tp              #pour mettre hp_rcc a zero
E_rcc = 196E9                  #module de young

# Inverseur
L_inv =  4.79                  #longueur des lames
h_inv = 3                      #epaisseur des lames
b_inv = 2                      #largeur de la base
M_inv = 100                    #masse 
E_inv = 190E9                  #module de young
g_inv = 40                     #distance entre le pivot et la tige
i_inv = 2                      #distance entre le centre de masse et le centre de rotation pivot
Theta_inv = m.atan(z_mt/g_inv) #angle dinclinaison de l'inverseur -----> a definir

#Tige superieur
d_ts = 10                      #diametre 
l_ts =  20                       #longueur 
E_ts = 140E9                    #module de young
v_ts = 0.3                      #module de poisson

#Tige inferieur
d_ti = 1                     #diametre 
l_ti =  30                      #longueur 
E_ti = 190E9                    #module de young
v_ti = 0.3                      #module de poisson

###########################################################################   CONVERSION

z_mt = z_mt/1000
z_max = z_max/1000
M_mt = M_mt/1000     
alphaVar = alphaVar*m.pi/180
alphaFixe = alphaFixe*m.pi/180 
h_m = h_m/1000         
D_m = D_m/1000      
a = a/1000           
yo_m = yo_m/1000      
M_m = M_m/1000        
R_m = R_m/1000        
h_s = h_s/1000        
d_s = d_s/1000        
yo_s = yo_s/1000       
l_eq = l_eq/1000              
d_eq = d_eq/1000              
M_cp = M_cp/1000              
b_tlp = b_tlp/1000            
h_tlp = h_tlp/1000            
l_tlp =  l_tlp/1000           
d_rcc = d_rcc/1000            
ht_rcc = ht_rcc/1000
l_rcc =  l_rcc/1000           
d_tp = d_tp/1000  
l_tp = l_tp/1000  
ho_tp = ho_tp/1000 
leff_tp = leff_tp/1000          
b_inv = b_inv/1000                    
L_inv =  L_inv/1000                   
h_inv = h_inv/1000                    
g_inv = g_inv/1000                    
i_inv = i_inv/1000                    
M_inv = M_inv/1000                    
d_ts = d_ts/1000                      
l_ts =  l_ts/1000                     
d_ti = d_ti/1000                      
l_ti =  l_ti/1000                     


############################################################################  DIMENSIONS
M_s = pmat_s *m.pi*((d_s+D_m/2)**2 - (D_m/2)**2)*h_s    #masse
I_inv = 2E-3                           #inverseur, a definir en fonction de tout jai envie de mourir


#Table a lames
print("DIMENSIONS")

print("Support")
print("a :", "{:.3e}".format(a))
print("Table a lames ")
print("Longueur des lames minimales :", "{:.3e}".format((3*z_max**2)/(5*0.6E-3)))

print('\n')

############################################################################  INERTIE DU SYSTEM




#Miroir et support
A_s = 2*h_s*d_s
A_m = h_m*D_m
yo = (yo_s*A_s + yo_m*A_m)/( A_m + A_s)  #axe neutre du miroir et support
Irot_s = 4.4478E-5  #pmat_s*m.pi*h_s*((a + D_m/2)**4 - (D_m/2)**4 )/4 + pmat_s*m.pi*h_s**3*((a + D_m/2)**2 - (D_m/2)**2 )/12
Irot_m = 2.7E-8  #M_m*( (D_m/2)**2) /4 + (M_m*h_m**2) /12
Irot_ms = Irot_m + Irot_s + A_m*(yo_m - yo)**2 +A_s*(yo_s - yo)**2

############################################################################  PARASITISME
#inverseur lames croise non separe
par_inv = m.sqrt(2) * L_inv * Theta_inv /15

#lames des tables
par_tlp = (3*z_mt**2)/5*l_tlp        #ce sera vers le meme cote des lames, Ox, Oy

#tiges superieurs des cotrepoids 
par_ts = 3*par_tlp/(5*l_ts)

#tiges inferieurs des cotrepoids 
par_ti = 3*par_tlp/(5*l_ti)
print("Parasitisme" )
print("inverseur : ", "{:.3e}".format(par_inv))
print("tables à lames : ", "{:.3e}".format(par_tlp))
print("tiges supérieures : ", "{:.3e}".format(par_ts))
print("tiges inférieurs : ", "{:.3e}".format(par_ti ))
print('\n')

####################################################################### PIVOT RCC ANGULAIRE#
deltaL_rcc = -m.sqrt(2)*(alphaVar**2 * l_rcc)/15     #vers -z

####################################################################### TIGE POLAIRE
I_tp = (d_tp**4)/12                 #inertie de la tige polaire
sigma_tp = E_tp*I_tp        #pour alleger les calculs
V_tp = (12*sigma_tp/(l_tp**3))*((l_tp*np.tan(alphaVar)/2) + ho_tp*np.sin(alphaVar))       #force tranchante sur la tige polaire
M_tp =  ((-6*sigma_tp/(l_tp**2)) * (ho_tp*np.sin(alphaVar) + l_tp*np.tan(alphaVar)/3))      #moment de flexion sur la tige

#Racourcissement de la tige polaire
racTigePolaire = -(1/(2 * sigma_tp**2) * (( V_tp**2 * l_tp**5)/20 + (V_tp * M_tp * l_tp**4)/4 + (M_tp**2 * l_tp**3)/3) ) 
racLevierPolaire = -ho_tp*(1 - np.cos(alphaVar))
deltaL_tp = racTigePolaire + racLevierPolaire         #vers -z

#deplacement du point d'incidence du faiseau
deplFaisceauInc = (h_m - ho_tp)*(1/np.cos(alphaVar) - 1)        #vers +z

#Mouvement parasite total
Zparasite = deplFaisceauInc + deltaL_tp + deltaL_rcc

#Plot du parasitisme tige polaire
if(True):
    limiteSup = [0.5E-6]* alphaVar.size
    limiteInf = [-0.5E-6]* alphaVar.size

    plt.plot(alphaVar, deltaL_tp, 'y')
    plt.plot(alphaVar, deplFaisceauInc, 'g')
    plt.plot(alphaVar, deltaL_rcc, 'b')
    plt.plot(alphaVar, Zparasite, 'r')
    plt.plot(alphaVar, limiteSup, 'r--')
    plt.plot(alphaVar, limiteInf, 'r--')
    plt.legend(['$\Delta$L', 'depl. faisceau incident','tige RCC','depl. z parasite total'])
    plt.xlabel("Angle d'inclinaison du miroir [°]")
    plt.ylabel("parasitisme deltaL [m]")
    axes = plt.gca()
    axes.set_ylim([-0.7E-6, 0.7E-6])
    plt.figtext(.3, .92, "$h_o$ = %smm"%'{:.2f}'.format(ho_tp*1000) )
    plt.figtext(.5, .92, "l = %smm"%'{:.2f}'.format(l_tp*1000) )
    plt.figtext(.7, .92, "d = %smm"%'{:.2f}'.format(d_tp*1000) )
    plt.savefig('line_plot.svg')  
    plt.grid()
    plt.show()

############################################################################  PARASITISME DE LA TIGE EQUATORIALE
F_eq = m.sqrt(2)* F_mt                  #force sur la tige
zp_eq = -m.sin(alphaFixe)*(D_m/2 + a)   #deplacement de la tige equatorielle selon z 
deltaL = (3/5*l_eq)*(m.sin(alphaFixe)**2)*(D_m/2 + a)**2    #Raccourcissement de la tige dans le plan xy, indesire
Beta = 2*deltaL/(D_m + 2*a)             #rotation parasite selon z
Beta = Beta * 180/m.pi

############################################################################  RIGIDITE EQUIVALENTE
K_tlp = ( 2*E_tlp*b_tlp*h_tlp**3 )/( l_tlp**3 )
if(tlpCharge == True):
    N_tlp_c = (8*m.pi*E_tlp*b_tlp*h_tlp**3)/(12*l_tlp**2)
    if(N_tlp >= N_tlp_c):
        print("CHARGE CRITIQUE DEPASSEE")
    N_tlp_o = (2*m.pi*E_tlp*b_tlp*h_tlp**3) / (12*l_tlp**2)
    K_tlp = K_tlp - (K_tlp * N_tlp)/N_tlp_o
    print("Charge critique: ", '{:.2f}'.format(N_tlp_c))
K_rcc_angulaire = (2*E_rcc*d_rcc**4 *(l_rcc**2)/ (3*l_rcc**3))
K_rcc_translation = (E_rcc*d_rcc**4) /(12*ht_rcc)
K_ts = (E_ts*d_ts**4)/(12*l_ts)
K_ti = (E_ti*d_ti**4)/(12*l_ti)
K_tp = E_tp*d_tp**4 *(l_tp**2 + 3*ho_tp*l_tp + 3*ho_tp)
K_inv = (2*E_inv*b_inv*h_inv**3)/3*L_inv
K_eq = 2*K_tlp + 2*K_rcc_translation/(R_m**2) + K_tp/(R_m**2) + K_rcc_angulaire/(R_m**2) + (2*K_ts)/(R_m**2) + K_ti/(g_inv**2) + K_inv/(g_inv**2)

############################################################################  MASSE REDUITE DU SYSTEM
c1 = (Irot_m + Irot_s + (M_m + M_s)*(ho_tp - yo)**2 ) /(R_m**2)    #pour alleger les calcules
c2 = M_cp
c3 =( I_inv + M_inv*(i_inv**2) )/(g_inv**2)
M_red = c1 + c2 + c3

############################################################################  LOIS DES MOUVEMENTS
#CALCUL
def convDegToZ(angle):
    return 3*180*angle/m.pi + 1
z_max = 4E-3                              #course max de l'actionneur
z_eq = 2.5E-3
z_nominal = 2.5E-3
z_precision_fin = 4E-3
z_precision_debut = 3.97E-3          #position nominal où le ressort equivalent n'excerce pas de force
F_eq = K_eq * (z_max - z_eq)
print("Rigidité equivalente : ""{:.3e}".format(K_eq))
print("Force equivalente : ""{:.3e}".format(F_eq))
a_max = (F_mt - F_eq)/(M_red + M_m)         #acceleration maximale, conversion
print("accélération max : ""{:.3e}".format(a_max))
theta_precision = 0.015                         #course angluaire pour le balayage de precision, en degree
z_precision = theta_precision*(m.pi/180)*R_m
T = m.sqrt(4*m.pi*(z_precision_fin - z_precision_debut)/a_max)            #periode pour passer de -z_pre a + z_pre
omega = 2*m.pi/T 
t_base = np.arange(0, T, 0.00001)

z_jrk = a_max*omega*np.cos(omega*t_base)        #jerk
z_acc = a_max*np.sin(omega*t_base)              #acceleration
z_vit = a_max/omega -(a_max/omega)*np.cos(omega*t_base)               #vitesse
z_pos = (a_max/omega)*t_base - (a_max/(omega**2))*np.sin(omega*t_base) + z_precision_debut     #position

f_scan = 1/(2*T)        #frequence complete de scan
f_scan_min = 860        #contrainte impose par le cahier de charge
T_scan_max = 1/(2*f_scan_min)

#PLOT
if(graphs):
    plt.plot(t_base, z_pos)
    plt.xlabel("Temps [ms]")
    plt.ylabel("Angle du miroir [deg]")
    plt.grid()
    plt.show()

    plt.plot(t_base, z_vit)
    plt.xlabel("Temps [ms]")
    plt.ylabel("Vitesse du miroir [deg/s]")
    plt.grid()
    plt.show()

    plt.plot(t_base, z_acc)
    plt.xlabel("Temps [ms]")
    plt.ylabel("Accélération du miroir [$deg/s^2$]")
    plt.grid()
    plt.show()

    plt.plot(t_base,z_jrk)
    plt.xlabel("Temps [ms]")
    plt.ylabel("Jerk du miroir [deg/$s^3$]")
    plt.grid()
    plt.show()
 
################################################################## VALEURS NUMERIQUES
print("Rotation parasite selon z a cause de la tige equatoriel : %s"%"{:.3e}".format(Beta), '\n')
print("Rigidité équivalente :%s"%"{:.3e}".format(K_eq))
print("Masse réduite :%s"%"{:.3e}".format(M_red))
print("Frequence de scan :%s"%"{:.3e}".format(f_scan))

