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

###############################################PROGRAM OPTIONS
graphs = True                   #print graphs ?

#Constantes
F_mt = 1.8      #force moteur
z_mt = 2.5      #course du moteur, position
M_mt = 7

#Parametres

alphaVar = np.arange(-1.5, 1.5, 0.001)        #en degrees, conversion apres, axe x
alphaFixe = 1.5                               #en degrees, conversion apres, axe x

#contrepoids
M_cp = 20.64

#miroir
D_m = 12.7                      #diametre miroir
h_m = 6.4                       #hauteur du miroir en mm
M_m = 2   
a= 47.75                      #masse du miroir
yo_m = 6.175#h_m/2                    #axe neutre par rapport à l'axe z de CATIA
R_ms = 47.75 #distance entre centre du miroir et application de la force par l'actionneur  
R_m = a + D_m/2                 #distance entre centre du miroir et application des rcc

#roue inertiel
CR_ri = 18.75E-3

R_ri_1 = 12.5E-3 #rayon roue 1
p_ri_1 = 8.94   #densité roue 1
CG_ri_1 = 12.5E-3

L_ri_2 = 3E-3 #rayon carré 2
p_ri_2 = p_ri_1 ##densité carré 2
h_ri_2 = 4E-3  #largeur carré 2
V_ri_2 = (L_ri_2**2)*h_ri_2
CG_ri_2 = 18.75E-3

#pivot RCC roue inertiel
l_ri_rcc= m.sqrt(2)*CR_ri
h_ri_rcc= 0.4E-3
sigma_d_ri = 500E6
p_ri_rcc = L_ri_2/m.sqrt(2)
E_ri = 114E9
SF_ri = 1.5
sigma_adm = sigma_d_ri / SF_ri
theta_adm_ri = ((l_ri_rcc**2)*sigma_adm)/(E_ri*((2*h_ri_rcc*l_ri_rcc) + (3*h_ri_rcc*p_ri_rcc)))
theta_adm_ri = theta_adm_ri *(180/m.pi)
print('theta_adm_ri', theta_adm_ri)

theta_min = m.arcsin()

#support
h_s = 1.95                         #epaisseur
d_s =  2*R_m                       #largeur    
A_s = 168.05                    #surface calculée par CATIA     
yo_s = 2.702        #axe neutre par rapport à l'axe z de CATIA
pmat_s = 1                      #densite materiau

# Tige Polaire
d_tp = 1                     #diametre de la tige en mm
l_tp = 15                        #longueur de la tige en mm
ho_tp = 5                    #hauteur du centre de rotation en mm
leff_tp = 0.699*l_tp            #longueur efficace, Pinned-fixed
E_tp = 196E9                    #module de young, en acier
v_tp = 0.3                      #module de poisson

# Masse inertielle
M_cp = 10                        #masse

# Table a lames paralleles x2
M_blocC_tlp = 13.64E-3
b_tlp =  0                    #largeur des lames
e_tlp =   0                    #distance entre les lames
l_tlp =   0                     #longueur des lames
E_tlp = 196E9                       #module de young
h_tlp = 0.5                       #epaisseur des lames
N_tlp = 350                     #charge du ressort sur les tlb
N_tlp_c = 1             #charge critique, calculé apres
N_tlp_o = 1             #parametre pour rigidite avec charge, calculé apres
tlpCharge = True                #sont-elles chargées?

# Tige equatoriale
l_eq = 30                       #longueur
d_eq = 1                      #diametre
E_eq = 196E9                      #module de young

# Pivot RCC
d_rcc = 0.5                      #diametre de la tige, egale à e et à b dans la formule
l_rcc = 12                     #longueur des lames
p_rcc = 3                     #hauteur du pivot virtuel
E_rcc = 196E9                   #module de young
h_rcc = h_s/(3*m.sin(45*m.pi/180))   #distance entre point de contact de la tige et le pivot virtuel
l_rcc_eff = m.sqrt(2)*l_rcc/2   #approximation des tiges rcc qui se trouvent sur le moteur actionné

# Inverseur
L_inv =  4.79                      #longueur des lames
h_inv = 3                       #epaisseur des lames
b_inv = 2                      #largeur de la base
M_inv = 100                       #masse 
E_inv = 190E9                   #module de young
g_inv = 40                       #distance entre le pivot et la tige
i_inv = 2                       #distance entre le centre de masse et le centre de rotation pivot
Theta_inv = m.atan(z_mt/g_inv)  #angle dinclinaison de l'inverseur -----> a definir

#Tige superieur
d_ts = 1                      #diametre 
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
M_mt = M_mt/1000  
M_cp = 20.64 /1000   
alphaVar = alphaVar*m.pi/180
alphaFixe = alphaFixe*m.pi/180 
h_m = h_m/1000         
D_m = D_m/1000      
R_ms = R_ms /1000       
yo_m = yo_m/1000      
M_m = M_m/1000        
R_m = R_m/1000  
A_s = A_s/1000000      
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
l_rcc =  l_rcc/1000           
p_rcc = p_rcc/1000            
h_rcc = h_rcc/1000            
l_rcc_eff = l_rcc_eff/1000
ho_tp = ho_tp/1000             
l_tp = l_tp/1000               
leff_tp = leff_tp/1000         
d_tp = d_tp/1000                      
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
 


############################################################################  INERTIE DU SYSTEME
M_s = pmat_s *m.pi*((d_s+D_m/2)**2 - (D_m/2)**2)*h_s    #masse
I_inv = 2E-3                           #inverseur, a definir en fonction de tout jai envie de mourir
I_tp = (d_tp**4)/12                 #inertie de la tige polaire

#Miroir et support
CR_ms = 1.5E-3  #position centre de rotation miroir+support par rapport au bas du support
CG_ms = 2.968E-3  #position du centre de gravité miroir+support par rapport au bas du support
I_rot_ms_CG = 2.243E-5 #moment d'inertie du miroir+support par rapport au centre de gravité
I_rot_ms_CR = I_rot_ms_CG + (M_m + M_s)*(CR_ms - CG_ms)**2

#roue inertie 2
I_ri_2 = (1/6)*(p_ri_2*V_ri_2*L_ri_2**2) +(p_ri_2*V_ri_2*(CG_ri_2-CR_ri)**2)

############################################################################ EQUILIBRAGE DU MOMENT INERTIE
# l'équilibrage du moment d'inertie nous fixe l'inertie de la roue inertiel et donc ses dimensions

I_roue_in = (R_ri_1/R_ms)*(I_rot_ms_CR + (M_blocC_tlp*R_ms**2) + (M_cp*R_ms**2))
denom = (0.5*m.pi*p_ri_1*(R_ri_1**4))+(p_ri_1*m.pi*(R_ri_1**2)*(CG_ri_1-CR_ri)**2)
h_ri_1 = (I_roue_in - I_ri_2)/denom
print('h_ri_1=')
print(h_ri_1)

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

####################################################################### PIVOT RCC ANGULAIRE
I_rcc = d_rcc**4 /12
sigma_rcc = E_rcc*I_rcc        #pour alleger les calculs
V_rcc = (12*sigma_rcc/l_rcc**3)*((l_rcc*np.tan(alphaVar)/2) + h_rcc*np.sin(alphaVar))       #force tranchante sur la tige polaire
M_rcc =  (-6*sigma_rcc/l_rcc**2 * (h_rcc*np.sin(alphaVar) + l_rcc*np.tan(alphaVar)/3))      #moment de flexion sur la tige

#Racourcissement d'une seule tige rcc 
racTigeRCC = (1/(2 * sigma_rcc**2) * (( V_rcc**2 * l_rcc**5)/20 + (V_rcc * M_rcc * l_rcc**4)/4 + (M_rcc**2 * l_rcc**3)/3) ) 
racLevierRCC = -h_rcc*(1 - np.cos(alphaVar))
deltaL_rcc = m.sqrt(2)*(-racTigeRCC - racLevierRCC)

####################################################################### TIGE POLAIRE
sigma_tp = E_tp*I_tp        #pour alleger les calculs
V_tp = (12*sigma_tp/l_tp**3)*((l_tp*np.tan(alphaVar)/2) + ho_tp*np.sin(alphaVar))       #force tranchante sur la tige polaire
M_tp =  (-6*sigma_tp/l_tp**2 * (ho_tp*np.sin(alphaVar) + l_tp*np.tan(alphaVar)/3))      #moment de flexion sur la tige

#Racourcissement de la tige polaire
racTigePolaire = (1/(2 * sigma_tp**2) * (( V_tp**2 * l_tp**5)/20 + (V_tp * M_tp * l_tp**4)/4 + (M_tp**2 * l_tp**3)/3) ) 
racLevierPolaire = -ho_tp*(1 - np.cos(alphaVar))
deltaL_tp = -racTigePolaire - racLevierPolaire

#deplacement du point d'incidence du faiseau
deplFaisceauInc = (h_m - ho_tp)*(1/np.cos(alphaVar) - 1)

#Mouvement parasite total
Zparasite = deplFaisceauInc + deltaL_tp + deltaL_rcc

#Plot du parasitisme tige polaire
if(graphs):
    limiteSup = [0.5E-6]* alphaVar.size
    limiteInf = [-0.5E-6]* alphaVar.size

    plt.plot(alphaVar, deltaL_tp, 'y')
    plt.plot(alphaVar, deplFaisceauInc, 'g')
    plt.plot(alphaVar, deltaL_rcc, 'b')
    plt.plot(alphaVar, Zparasite, 'r')
    plt.plot(alphaVar, limiteSup, 'r--')
    plt.plot(alphaVar, limiteInf, 'r--')
    plt.legend(['$\Delta$L', 'depl. faisceau incident','depl. parasite RCC' ,'depl. z parasite total'])
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

############################################################################  MASSE REDUITE DU SYSTEM

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