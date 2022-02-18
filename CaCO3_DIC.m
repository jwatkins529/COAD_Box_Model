function [alpha_c,alpha_o,alpha63,alpha64,alpha65,R_c]=CaCO3_DIC(TC,Ks,pH, A, B0, B2, B1)
%%
%======Input conditions====================================================
TK = TC+273.15;
H = 10^-pH;
B = B1+B2;  
omega = A.*B1./Ks;
%======Input parameters====================================================
theta = 10^(8.6-pH);                                                                %Ratio of HCO3/CO3 on calcite surface
phi = B2/B1;                                                                        %Ratio of HCO3/CO3 in bulk solution
chi = phi/(1+phi);                                                                  %Fraction of HCO3 in EIC
boltz = 1.38065E-23;                                                                %Boltzmann constant
epsilon = 0.67E-20;                                                                 %Kink formation energy
gamma = 1.2E-19;                                                                    %Edge work
kA1 = 3E6;                                                                          %Wolthers et al. (2012)
kA2 = kA1;                                                                          %Wolthers et al. (2012)
kB1 = 2*kA1*(1+theta)/(1+phi);                                                      %Wolthers et al. (2012)
kB2 = kB1;                                                                          %Wolthers et al. (2012)
vA1 = 2E3;                                                                          %Wolthers et al. (2012)
vA2 = vA1;                                                                          %Wolthers et al. (2012)
%=====Dependent parameters=================================================
a = 3.199E-10;                                                                      %Closest spacing between A and B sites
d = 27100;                                                                          %Molar density of calcite
S = (A*B1/Ks)^0.5;                                                                  %Saturation ratio for calcite
k_bar_A = kA1+theta*kA2;                                                            %Rate coefficient for A attachment
k_bar_B = kB1+phi*kB2;                                                              %Rate coefficient for B1 and B2 attachment
v_bar_A = vA1+vA2;                                                                  %Rate coefficient for A detachment
vB1 = Ks*k_bar_A*k_bar_B/(v_bar_A*(1+theta));                                       %B1 detachment frequency
vB2 = vB1;                                                                          %Wolthers et al. (2011)
v_bar_B = vB1+theta*vB2;                                                            %Wolthers et al. (2011)
PB1 = (k_bar_B*B1+v_bar_A)/(k_bar_A*A+v_bar_B+(1+theta)*(k_bar_B*B1+v_bar_A));      %Probability that a given site is a B1 site
PA = 1-(1+theta)*PB1;                                                               %Probability that a given site is an A site
PB2 = 1-PA-PB1;                                                                     %Probability that a given site is a B2 site
uA = k_bar_A*A*PB1-v_bar_A*PA;                                                      %Kink propagation rate
uB = k_bar_B*B1*PA-v_bar_B*PB1;                                                     %Kink propagation rate
uC = uA+uB;                                                                         %Kink propagation rate
iA = 2*exp(-2*epsilon/(boltz*TK))*(S^2-1)*(v_bar_B*k_bar_A*A/(k_bar_A*A+v_bar_B));  %Rate of kink formation on B sites
iB = 2*exp(-2*epsilon/(boltz*TK))*(S^2-1)*(v_bar_A*k_bar_B*B/(k_bar_B*B+v_bar_A));  %Rate of kink formation on A sites
iC = (iA+iB)/2;                                                                     %Net rate of kink formation
rho_c = (2*iC/(uA+uB))^0.5;                                                         %Steady state kink density
y_o = 19*a*gamma/(boltz*TK*log(S));                                                 %Step spacing
R_c = rho_c*uC*a^2*d/y_o;                                                           %Calcite growth rate
%======Oxygen Isotopes=====================================================
d18Ow = 0;
rVSMOW = 0.0020052;
rw = (d18Ow/1000+1)*rVSMOW;
B2 = phi*B1;
alpha_dw = exp(2520/TK^2 + 0.01212);                                                %Beck et al. (2005)
alpha_bw = exp(2590/TK^2 + 0.00189);                                                %Beck et al. (2005)
alpha_cw = exp(2390/TK^2 - 0.00270);                                                %Beck et al. (2005)
alpha_xw = exp((17747/TK-29.777)/1000);
alpha_O_eq_1 = alpha_xw/alpha_cw;
alpha_O_eq_2 = alpha_xw/alpha_bw;
alpha_O_f_1 = 0.9995;                                                               %KFF
alpha_O_f_2 = 0.9966;                                                               %KFF
alpha_O_b_1 = alpha_O_f_1/alpha_O_eq_1;
alpha_O_b_2 = alpha_O_f_2/alpha_O_eq_2;
alpha_O_b2_b1 = alpha_bw/alpha_cw;
alpha_eq_c_w = alpha_xw;
%=========Carbon isotopes =================================================
d13C_DIC = -17;
rDIC = 0.01118;                                                                     %This is PDB, but the value is arbitrary
epsilon_gb = -9483/TK+23.89;                                                        %Mook (1986)
epsilon_cb = -867/TK+2.52;                                                          %Mook (1986) - small offset between HCO3-CO3
epsilon_db = -9866/TK+24.12;                                                        %Mook (1986)                                                  
delta_B2 = (d13C_DIC*(B0+B1+B2)-(epsilon_db*B0+epsilon_cb*B1))/((1+epsilon_db/1000)*B0+B2+(1+epsilon_cb/1000)*B2);
delta_B1 = delta_B2*(1+epsilon_cb/1000)+epsilon_cb;
alpha_gb = epsilon_gb/1000+1;
alpha_gx = exp((-2.4612+(7.6663*1000/TK)-(2.988*1000000/(TK^2)))/1000);
alpha_xb = alpha_gb/alpha_gx;
alpha_cb = epsilon_cb/1000+1;
alpha_bc = 1/alpha_cb;
alpha_xc = alpha_xb/alpha_cb;
alpha_C_eq_1 = alpha_xc ;
alpha_C_eq_2 = alpha_xb;
alpha_C_f_1 = 1.0000;
alpha_C_f_2 = 1.0000;
alpha_C_b_1 = alpha_C_f_1/alpha_C_eq_1;
alpha_C_b_2 = alpha_C_f_2/alpha_C_eq_2;
% =====Isotopologue k's and nu's ==========================================
d1 = 3*alpha_eq_c_w*rw;
r_B1 = (delta_B1/1000+1)*rDIC;
r_B2 = (delta_B2/1000+1)*rDIC;
d2 = alpha_C_eq_2*r_B2;
k1 = kB1;
k2 = kB2;
k3 = kB1*alpha_C_f_1;                                                               %13kB1
k4 = kB2*alpha_C_f_2;                                                               %13kB2
k5 = kB1*alpha_O_f_1;                                                               %18kB1
k6 = kB2*alpha_O_f_2;                                                               %18kB2
v1 = vB1;
v2 = vB2;
v3 = vB1*alpha_C_b_1;                                                               %13vB1
v4 = vB2*alpha_C_b_2;                                                               %13vB2
v5 = vB1*alpha_O_b_1;                                                               %18vB1
v6 = vB2*alpha_O_b_2;                                                               %18vB2
%=================Isotopologue concentrations =============================
B3 = B1*(delta_B1/1000+1)*rDIC;                                                     %Mass 61
B4 = B2*(delta_B2/1000+1)*rDIC;                                                     %Mass 61
B5 = 3*alpha_cw*rw*B1;                                                              %Mass 62
B6 = 3*alpha_bw*rw*B2;                                                              %Mass 62
%===========Carbon isotope composition of calcite==========================
D = ((k3*B3+k4*B4)/(k1*B1+k2*B2)) * (v1*PB1+v2*PB2) / (d2*(v3-v4))-v4*(PB1+PB2)/(v3-v4);    
uB = k1*B1*PA+k2*B2*PA-v1*PB1-v2*PB2;                                                       
r_c = (k3*B3*PA+k4*B4*PA)/(uB+v3*D+v4*PB1+v4*PB2-v4*D);
alpha_c = r_c/r_B1;
delta_c = (r_c/rDIC-1)*1000;
%============Oxygen isotope composition of calcite=========================
D = ((k5*B5+k6*B6)/(k1*B1+k2*B2)) * (v1*PB1+v2*PB2)/(3*alpha_eq_c_w*rw*(v5-v6))-v6*(PB1+PB2)/(v5-v6);   %Eq. 14 of Watkins et al. (2014)
uB = k1*B1*PA+k2*B2*PA-v1*PB1-v2*PB2;                                                                   %Eq. 6 of Watkins et al. (2014)
r_o = (k5*B5*PA+k6*B6*PA)/(uB+(v5-v6)*D+v6*(PB1+PB2));                                                  %Eq. 15b of Watkins et al. (2014)
alpha_o = r_o/(B5/B1);
%alpha_o = r_o/3/rw;
%============Clumped isotope composition of calcite========================
D47_B0 = 26447/(TK^2)+285.51/TK-0.3004;                                             %Wang et al. (2004) - assume CO2(aq) = CO2(g), as done by Guo (2020)
D48_B0 = 29306/(TK^2)+93.885/TK-0.2914;                                             %Wang et al. (2004) - assume CO2(aq) = CO2(g), as done by Guo (2020)
D49_B0 = 108776/(TK^2)+477.14/TK-0.5954;                                            %Wang et al. (2004) - assume CO2(aq) = CO2(g), as done by Guo (2020)
%D63_B1 = 9.03e-6*TC^2-(3.13e-3)*TC+7.23e-1 - 0.28;                                 %Uchikawa et al. (2021) shifted from Hill et al. (2020) - negative AFF to convert to D63
%D63_B2 = 9.39e-6*TC^2-(3.31e-3)*TC+7.91e-1 - 0.28;                                 %Uchikawa et al. (2021) shifted from Hill et al. (2020) - negative AFF to convert to D63
D63_B1 = 43187/(TK^2)-34.833/TK+0.0007;                                             %Hill et al. 2020
D63_B2 = 43655/(TK^2)-23.643/TK-0.0088;                                             %Hill et al. 2020
D64_B1 = 23492/(TK^2)-52.842/TK+0.0304;                                             %Hill et al. 2020
D64_B2 = 21842/(TK^2)-50.457/TK+0.0291;                                             %Hill et al. 2020
D65_B1 = 112667/(TK^2)-123.11/TK+0.0304;                                            %Hill et al. 2020
D65_B2 = 112026/(TK^2)-97.208/TK+0.009;                                             %Hill et al. 2020
D63eq = 43159/(TK^2)-25.095/TK-0.0078;                                              %Calcite; Hill et al. 2020
D64eq = 23566/(TK^2)-52.319/TK+0.0297;                                              %Calcite; Hill et al. 2020
D65eq = 112667/(TK^2)-102.28/TK+0.012;                                              %Calcite; Hill et al. 2020
%==================13C-18O clumped isotopes=================================
epsilon = 0;
alpha_63_f_1 = alpha_O_f_1*alpha_C_f_1+epsilon;
alpha_63_f_2 = alpha_O_f_2*alpha_C_f_2+epsilon;
AttRate = k1*B1*PA+k2*B2*PA;
DetRate = v1*PB1+v2*PB2;
AR = k1*B1*PA+k2*B2*PA;
DR = v1*PB1+v2*PB2;
carbon_isotopes = alpha_xc*(AR-DR)/((1-chi)*alpha_C_f_1+chi*alpha_C_f_2*alpha_bc)+DR;
oxygen_isotopes = alpha_xw/alpha_cw*(AR-DR)/((1-chi)*alpha_O_f_1+chi*alpha_O_f_2*alpha_bw/alpha_cw)+DR;
RR63 = carbon_isotopes*oxygen_isotopes*(D63eq/1000+1)/(AR*(DR+alpha_xc*alpha_xw/alpha_cw*(D63eq/1000+1)*(AR-DR)/((1-chi)*alpha_63_f_1*(D63_B1/1000+1)+chi*alpha_63_f_2*alpha_bc*alpha_bw/alpha_cw*(D63_B2/1000+1))));
D63_EIC = B2/(B1+B2)*D63_B2 + B1/(B1+B2)*D63_B1;
alpha63 = RR63/(D63_EIC/1000+1);
%===============18O-18O clumped isotopes====================================
epsilon = 0;
alpha_64_f_1 = alpha_O_f_1*alpha_O_f_1+epsilon;
alpha_64_f_2 = alpha_O_f_2*alpha_O_f_2+epsilon;
D64_EIC = B2/(B1+B2)*D64_B2 + B1/(B1+B2)*D64_B1;
RR64 = oxygen_isotopes*oxygen_isotopes*(D64eq/1000+1)/(AR*(DR+(alpha_xw/alpha_cw)^2*(D64eq/1000+1)*(AR-DR)/((1-chi)*alpha_64_f_1*(D64_B1/1000+1)+chi*alpha_64_f_2*(alpha_bw/alpha_cw)^2*(D64_B2/1000+1))));
alpha64 = RR64/(D64_EIC/1000+1);
%===============13C-18O-18O clumped isotopes====================================
epsilon = 0;
alpha_65_f_1 = alpha_O_f_1*alpha_O_f_1*alpha_C_f_1+epsilon;
alpha_65_f_2 = alpha_O_f_2*alpha_O_f_2*alpha_C_f_2+epsilon;
D65_EIC = B2/(B1+B2)*D64_B2 + B1/(B1+B2)*D65_B1;
RR65 = carbon_isotopes*oxygen_isotopes^2*(D65eq/1000+1)/(AR^2*(DR+alpha_xc*(alpha_xw/alpha_cw)^2*(D65eq/1000+1)*(AR-DR)/((1-chi)*alpha_65_f_1*(D65_B1/1000+1)+chi*alpha_64_f_2*alpha_bc*(alpha_bw/alpha_cw)^2*(D65_B2/1000+1))));
alpha65 = RR65/(D65_EIC/1000+1);%now refers to r/r* of CaCO3 relative to r/r* of EIC
%========================================================================
end