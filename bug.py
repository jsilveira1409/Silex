K_tlp = ( 2*E_tlp*b_tlp*h_tlp**3 )/( l_tlp**3 )
if(tlpCharge == True):
    N_tlp_c = (8*m.pi*E_tlp*b_tlp*h_tlp**3)/(12*l_tlp**2)
    if(N_tlp >= N_tlp_c):
        print("CHARGE CRITIQUE DEPASSEE")
    N_tlp_o = (2*m.pi*E_tlp*b_tlp*h_tlp**3) / (12*l_tlp**2)
    K_tlp = K_tlp - (K_tlp * N_tlp)/N_tlp_o
    print("Charge critique: ", '{:.2f}'.format(N_tlp_c))
K_rcc_angulaire = (2*E_rcc*h_rcc**4 *(l_rcc**2 + 3*p_rcc*l_rcc + 3*p_rcc**2)/ (3*l_rcc**3))
K_rcc_translation = (E_rcc*h_rcc**4) /(12*l_rcc_eff)
K_ts = (E_ts*d_ts**4)/(12*l_ts)
K_ti = (E_ti*d_ti**4)/(12*l_ti)
K_tp = E_tp*d_tp**4 *(l_tp**2 + 3*ho_tp*l_tp + 3*ho_tp)
K_inv = (2*E_inv*b_inv*h_inv**3)/3*L_inv
K_eq = 2*K_tlp + 2*K_rcc_translation/(R_m**2) + K_tp/(R_m**2) + K_rcc_angulaire/(R_m**2) + (2*K_ts)/(R_m**2) + K_ti/(g_inv**2) + K_inv/(g_inv**2)


#masse reduite
c1 = (Irot_m + Irot_s + (M_m + M_s)*(ho_tp - yo)**2 ) /(R_m**2)    #pour alleger les calcules
c2 = M_cp
c3 =( I_inv + M_inv*(i_inv**2) )/(g_inv**2)
M_red = c1 + c2 + c3