%% File: Ion_by_ion_D61_D62_D63_D64_D65.m
%Created: 10/14/2021
%Last modified: 1/10/2022
%Author: James M. Watkins
%Department of Earth Sciences, University of Oregon
%watkins4@uoregon.edu
%Description: This script was used to produce Fig. 2 in Watkins and Devriendt
%(2022). If you use it, please cite: Watkins, J. and Devriendt, L., 2022, A
%combined model for kinetic clumped isotope effects in the
%CaCO3-DIC-H2O system, Geochemistry, Geophysics, Geosystems, v. X, p. Y-Z. 

%% Description
% Ion-by-ion model updated from Watkins and Hunt (2015)
% Update #1: This code uses Beck et al. (2005) instead of Wang et
% al. (2013) for the equilibrium d18O of CO2, HCO3, and CO3. This leads to
% different KFFs for the HCO3 --> CaCO3 and CO3 --> CaCO3 reactions
% Update #2: This code produces curves for D64 and D65. 
clear all
clc
%%
%=========================Input conditions=================================
TC = 5;
TK = TC+273.15;
salinity = 0;                                                                %permil
mNaCl = 0.015;                                                               %(seawater mNaCl = 0.537) (15 mM from Morrill et al., 2013)
A = 0.015;                                                                   %[Ca2+] concentration (M)
Rate_spec = 10^-6;                                                           %Use 10^-5.5 for Romanek data and KFFs = 1.  Use 10-6 for Watkins data. 
logOmega = 0.0005:0.05:8;                                                    %Range of omegas to produce a wide range of growth rates
Omega = exp(logOmega);
niter=1;
pH_stored = [0 0.1 14]';
%%
%==========================Initializing====================================
%------------------------ion-by-ion parameters-----------------------------
S = zeros(length(Omega),1);                                                 %Saturation ratio of Wolthers et al. (2012)
B = zeros(length(Omega),1);                                                 %HCO3+CO3
B0 = zeros(length(Omega),1);                                                %CO2
B1 = zeros(length(Omega),1);                                                %CO3
B2 = zeros(length(Omega),1);                                                %HCO3
DIC = zeros(length(Omega),1);                                               %DIC    
PA = zeros(length(Omega),1);                                                %Probability that a given site is an A site
PB1 = zeros(length(Omega),1);                                               %Probability that a given site is a B1 site
PB2 = zeros(length(Omega),1);                                               %Probability that a given site is a B2 site
uA = zeros(length(Omega),1);                                                %Kink propagation rate
uB = zeros(length(Omega),1);                                                %Kink propagation rate
uB_prime = zeros(length(Omega),1);                                          %Kink propagation rate for C18O16O16O2-
uC = zeros(length(Omega),1);                                                %Kink propagation rate
iA = zeros(length(Omega),1);                                                %Rate of kink formation on B sites
iB = zeros(length(Omega),1);                                                %Rate of kink formation on A sites
iC = zeros(length(Omega),1);                                                %Net rate of kink formation
rho_c = zeros(length(Omega),1);                                             %Steady state kink density
y_o = zeros(length(Omega),1);                                               %Step spacing
AR = zeros(length(Omega),1);                                                %Net attachment rate
DR = zeros(length(Omega),1);                                                %Net detachment rate
R_c = zeros(length(Omega),1);                                               %Calcite growth rate in moles/m2/s
Rate = zeros(length(Omega),1);                                              %log10R_c
%------------------------oxygen isotopes-----------------------------------
oxygen_isotopes = zeros(length(Omega),1);
B1_prime = zeros(length(Omega),1);
B2_prime = zeros(length(Omega),1);
PB1_prime = zeros(length(Omega),1);
PB2_prime = zeros(length(Omega),1);
Delta_dw = zeros(length(Omega),1);
Delta_bw = zeros(length(Omega),1);
Delta_cw = zeros(length(Omega),1);
R_o = zeros(length(Omega),1);                                               %2866/2666 of calcite
alpha_o = zeros(length(Omega),1);                                           %alpha-cw
D18O = zeros(length(Omega),1);                                              %1000ln18alpha_cw
%------------------------carbon isotopes-----------------------------------
carbon_isotopes = zeros(length(Omega),1);
Delta_cb = zeros(length(Omega),1);
Delta_bb = zeros(length(Omega),1);
delta_B0 = zeros(length(Omega),1);
delta_B1 = zeros(length(Omega),1);
delta_B2 = zeros(length(Omega),1);
D = zeros(length(Omega),1);
r_c = zeros(length(Omega),1);                                               %3666/2666 of calcite
delta_c = zeros(length(Omega),1);                                           %d13C of calcite (assumes d13C of DIC is 0)
alpha_cB2 = zeros(length(Omega),1);                                         %13alpha_c-B2
alpha_c_B2=zeros(length(Omega),1);
%------------------------Clumped isotopes-----------------------------------
d2 = zeros(length(Omega),1);                                                %Used for calculating isotopologue k's and nu's
r_B1 = zeros(length(Omega),1);                                              %Used for calculating isotopologue k's and nu's
r_B2 = zeros(length(Omega),1);                                              %Used for calculating isotopologue k's and nu's
B3 = zeros(length(Omega),1);                                                %Mass 61 CO3
B4 = zeros(length(Omega),1);                                                %Mass 61 HCO3
B5 = zeros(length(Omega),1);                                                %Mass 62 CO3
B6 = zeros(length(Omega),1);                                                %Mass 62 HCO3
term1 = zeros(length(Omega),1);
term2 = zeros(length(Omega),1);
RR63 = zeros(length(Omega),1);                                              %63R/63R*
Delta63 = zeros(length(Omega),1);
RR64 = zeros(length(Omega),1);                                              %64R/64R*
Delta64 = zeros(length(Omega),1);
RR65 = zeros(length(Omega),1);                                              %65R/65R*
Delta65 = zeros(length(Omega),1);
%================Arrays for plotting vs pH=================================
Rate_spec_stored = zeros(length(pH_stored),1);                              %Specified growth rate at top of code
colorcode_stored = zeros(length(pH_stored),1);                              %log10R for coloring the curves
%----------------------oxygen isotopes-------------------------------------
D_CO2_water = zeros(length(pH_stored),1);
D_B1_water = zeros(length(pH_stored),1);
D_B2_water = zeros(length(pH_stored),1);
D_DIC_water = zeros(length(pH_stored),1);
alpha_stored_slow = zeros(length(pH_stored),1);                             %18alpha_c-w slow growth
alpha_stored_fast = zeros(length(pH_stored),1);                             %18alpha_c-w fast growth limit
D18O_stored_slow = zeros(length(pH_stored),1);
D18O_stored_fast = zeros(length(pH_stored),1);
alpha_DIC_stored = zeros(length(pH_stored),1);                              %18alpha_DIC-w
D18O_spec_stored = zeros(length(pH_stored),1);                              %1000lna_c-w at specified growth rate
%----------------------carbon isotopes-------------------------------------
D_B0 = zeros(length(pH_stored),1);                                          %1000lna_CO2-DIC
D_B2 = zeros(length(pH_stored),1);                                          %1000lna_HCO3-DIC
D_B1 = zeros(length(pH_stored),1);                                          %1000lna_CO3-DIC
D_DIC = zeros(length(pH_stored),1);                                         %1000lna_DIC-DIC
D13C_stored_slow = zeros(length(pH_stored),1);                              %1000lna_c-DIC at slow growth  
D13C_stored_fast = zeros(length(pH_stored),1);                              %1000lna_c-DIC at fast growth  
D13C_spec_stored = zeros(length(pH_stored),1);                              %1000lna_c-DIC at specified growth rate  
%----------------------clumped isotopes-------------------------------------
D63_B1_stored  = zeros(length(pH_stored),1);                                %D63 of CO3(eq)
D63_B2_stored  = zeros(length(pH_stored),1);                                %D63 of HCO3(eq)
D63_calcite_equilibrium = zeros(length(pH_stored),1);                       %D63 of calcite(eq)
D63_calcite_kinetic = zeros(length(pH_stored),1);                           %D63 of calcite at fast growth
D63_spec_stored = zeros(length(pH_stored),1);                               %D63 of calcite at specified growth rate
D64_B1_stored  = zeros(length(pH_stored),1);                                %D64 of CO3(eq)
D64_B2_stored  = zeros(length(pH_stored),1);                                %D64 of HCO3(eq)
D64_calcite_equilibrium = zeros(length(pH_stored),1);                       %D64 of calcite(eq)
D64_spec_stored = zeros(length(pH_stored),1);                               %D64 of calcite at specified growth rate
D65_B1_stored  = zeros(length(pH_stored),1);                                %D65 of CO3(eq)
D65_B2_stored  = zeros(length(pH_stored),1);                                %D65 of HCO3(eq)
D65_calcite_equilibrium = zeros(length(pH_stored),1);                       %D65 of calcite(eq)
D65_spec_stored = zeros(length(pH_stored),1);                               %D65 of calcite at specified growth rate
%==========Arrays for plotting isotope comp of DIC=========================
xB1 = zeros(length(pH_stored),1);
xDIC = zeros(length(pH_stored),1);
xB2 = zeros(length(pH_stored),1); 
xB0 = zeros(length(pH_stored),1);
xDelta_EIC = zeros(length(pH_stored),1);
xD63_EIC = zeros(length(pH_stored),1);
xD64_EIC = zeros(length(pH_stored),1);
xD65_EIC = zeros(length(pH_stored),1);
%=========Matrices for plotting rate dependence at different pH values=====
% Rate_stored = zeros(length(Omega),niter);
% calcite_stored = zeros(length(Omega),niter);
% carbon_stored = zeros(length(Omega),niter);
% Delta63_stored = zeros(length(Omega),niter);
% Delta_B1_stored(i) = zeros(length(Omega),niter);
% Delta_B2_stored(i) = zeros(length(Omega),niter);
% Delta_B1_carbon_stored(i) = zeros(length(Omega),niter);
% Delta_B2_carbon_stored(i) = zeros(length(Omega),niter);
%%
iter=1;
for pH = 0:0.1:14
for i = 1: length(Omega)                                                    %This inner loop calculates the growth rate dependence at each pH
%%
%-------------Millero et al. (2007) - NaCl solutions-----------------------
AA1 = 35.2911*mNaCl^0.5+0.8491*mNaCl-0.32*mNaCl^1.5+0.055*mNaCl^2;
BB1 = -1583.09*mNaCl^0.5;
CC1 = -5.4366*(mNaCl^0.5);
pK1o = -114.3106+5773.67/TK+17.779524*log(TK);
pK1 = AA1+BB1/TK+CC1*log(TK)+pK1o;                                           
AA2 = 38.2746*mNaCl^0.5+1.6057*mNaCl-0.647*mNaCl^1.5+0.113*mNaCl^2;
BB2 = -1738.16*mNaCl^0.5;
CC2 = -6.0346*mNaCl^0.5;
pK2o = -83.2997+4821.38/TK+13.5962*log(TK);
pK2 = AA2+BB2/TK+CC2*log(TK)+pK2o;                                         
Ks = 10^-8.48;                                                              %Jacobsen (1974) calculated by Ellen using PHREEQC for freshwater
K1 = 10^(-pK1);
K2 = 10^(-pK2);
%%
%===========================Speciation=====================================
H = 10^-pH;
B1(i) = Omega(i)*Ks/A;% [CO32-] (M)
DIC(i) = B1(i)*(1+H/K2+(H^2)/(K1*K2));
B2(i) = DIC(i)/(1+H/K1+K2/H);% [HCO3-] (M)
B0(i) = DIC(i)/(1+K1/H+K1*K2/(H^2)); %[CO2 + H2CO3] (M)
B(i) = B1(i)+B2(i);
theta = 10^(8.6-pH);                                                        %Ratio of HCO3/CO3 on calcite surface
phi = (1+H/K2+H*H/K1/K2)/(1+H/K1+K2/H);                                     %Ratio of HCO3/CO3 in bulk solution
chi = phi/(1+phi);
%========================== Zeebe 2007 ====================================
%-----Note: THIS IS ONLY USED FOR PLOTTING THE DIC d18O CURVE--------------
CO2s = DIC(1)/(1+K1/H+K1*K2/(H^2));
xCO2 = CO2s/(10^(-2.98)+1);
H2CO3 = xCO2*10^(-2.98);
d1=CO2s;
a1=H2CO3;
b=DIC(1)/(1+H/K1+K2/H);
cr=DIC(1)/(1+H/K2+H*H/K1/K2);
V = 2*d1+3*(a1+b+cr);
zd = 2*d1/V;
za = 3*a1/V;
zb = 3*b/V;
zc = 3*cr/V;
lnad1000 = 1000*log(exp(2520/TK^2 + 0.01212));                               %Beck et al. (2005)
lnab1000 = 1000*log(exp(2590/TK^2 + 0.00189));                               %Beck et al. (2005)
lnac1000 = 1000*log(exp(2390/TK^2 - 0.00270));                               %Beck et al. (2005)
aa = 1.04;
ad = exp(lnad1000/1000);
ab = exp(lnab1000/1000);
ac = exp(lnac1000/1000);
RH2O = 1;
RRd = ad*RH2O;
RRa = aa*RH2O;
RRb = ab*RH2O;
RRc = ac*RH2O;
rd = RRd/(1+RRd);
ra = RRa/(1+RRa);
rb = RRb/(1+RRb);
xrc = RRc/(1+RRc);
RT =(rd*zd+ra*za+rb*zb+xrc*zc)/(1-(rd*zd+ra*za+rb*zb+xrc*zc));
AT = RT/RH2O;
lnAT1000 = 1000*log(AT);
%%
%===========================Input parameters===============================
boltz=1.38065E-23;                                                          %Boltzmann constant
epsilon = 0.67E-20;                                                         %kink formation energy
gamma = 1.2E-19;                                                            %edge work
kA1 = 3E6;                                                                  %Wolthers et al. (2012)
kA2 = kA1;                                                                  %Wolthers et al. (2012)
kB1 = 2*kA1*(1+theta)/(1+phi);                                              %Wolthers et al. (2012)
kB2 = kB1;                                                                  %Wolthers et al. (2012)
vA1 = 2E3;                                                                  %Wolthers et al. (2012)
vA2 = vA1;                                                                  %Wolthers et al. (2012)
%%
%========================Dependent parameters==============================
a = 3.199E-10;                                                              %Closest spacing between A and B sites
d = 27100;                                                                  %Molar density of calcite
S(i) = (A*B1(i)/Ks)^0.5;                                                    %Saturation ratio for calcite
k_bar_A = kA1+theta*kA2;                                                    %Rate coefficient for A attachment
k_bar_B = kB1+phi*kB2;                                                      %Rate coefficient for B1 and B2 attachment
v_bar_A = vA1+vA2;                                                          %Rate coefficient for A detachment
vB1 = Ks*k_bar_A*k_bar_B/(v_bar_A*(1+theta));                               %B1 detachment frequency
vB2 = vB1;                                                                  %Wolthers et al. (2012)
v_bar_B = vB1+theta*vB2;                                                    %Wolthers et al. (2012)
PB1(i) = (k_bar_B*B1(i)+v_bar_A)/(k_bar_A*A+v_bar_B+(1+theta)*(k_bar_B*B1(i)+v_bar_A));% Probability that a given site is a B1 site
PA(i) = 1-(1+theta)*PB1(i);                                                 %Probability that a given site is an A site
PB2(i) = 1-PA(i)-PB1(i);                                                    %Probability that a given site is a B2 site
uA(i) = k_bar_A*A*PB1(i)-v_bar_A*PA(i);                                     %Kink propagation rate
uB(i) = k_bar_B*B1(i)*PA(i)-v_bar_B*PB1(i);                                 %Kink propagation rate
uC(i) = uA(i)+uB(i);                                                        %Kink propagation rate
iA(i) = 2*exp(-2*epsilon/(boltz*TK))*(S(i)^2-1)*(v_bar_B*k_bar_A*A/(k_bar_A*A+v_bar_B));%Rate of kink formation on B sites
iB(i) = 2*exp(-2*epsilon/(boltz*TK))*(S(i)^2-1)*(v_bar_A*k_bar_B*B(i)/(k_bar_B*B(i)+v_bar_A));%Rate of kink formation on A sites
iC(i) = (iA(i)+iB(i))/2;                                                    %Net rate of kink formation
rho_c(i) = (2*iC(i)/(uA(i)+uB(i)))^0.5;                                     %Steady state kink density
y_o(i) = 19*a*gamma/(boltz*TK*log(S(i)));                                   %Step spacing
R_c(i) = rho_c(i)*uC(i)*a^2*d/y_o(i);                                       %Calcite growth rate
%%
%=========================Oxygen Isotopes==================================
rw = 0.00208835 / 1.01025;
alpha_dw = exp(2520/TK^2 + 0.01212);                                        %Beck et al. (2005)
alpha_bw = exp(2590/TK^2 + 0.00189);                                        %Beck et al. (2005)
alpha_cw = exp(2390/TK^2 - 0.00270);                                        %Beck et al. (2005)
alpha_xw = exp((17747/TK-29.777)/1000);                                     %Watkins et al. (2013)
alpha_O_eq_1 = alpha_xw/alpha_cw;                                           %EFF
alpha_O_eq_2 = alpha_xw/alpha_bw;                                           %EFF
alpha_O_f_1 = 0.9995;                                                       %KFF
alpha_O_f_2 = 0.9966;                                                       %KFF
alpha_O_b_1 = alpha_O_f_1/alpha_O_eq_1;                                 
alpha_O_b_2 = alpha_O_f_2/alpha_O_eq_2; 
alpha_O_b2_b1 = alpha_bw/alpha_cw;
alpha_eq_c_w = exp((1000*log(alpha_O_eq_1)+1000*log(ac))/1000);
Delta_cw(i) = 1000*log(alpha_cw);
Delta_bw(i) = 1000*log(alpha_bw);
%%
%==========================Carbon isotopes ================================
d13C_DIC = 0;
rDIC = 0.01118;                                                             %This is PDB, but the value is arbitrary
epsilon_bg = -0.1141*TC + 10.78;                                            %Zhang et al. (1995) HCO3 - CO2(g)
epsilon_dg = 0.0049*TC - 1.31;                                              %Zhang et al. (1995) CO2(aq) - CO2(g)
epsilon_cb = -867/TK+2.52;                                                  %Mook (1986) - 0.4 permil difference between CO3 and HCO3
epsilon_db = -9866/TK+24.12;                                                %Mook (1986) - CO2(aq)-HCO3 (is nearly identical to that in Zhang et al. 1995)
delta_B2(i) = (d13C_DIC*(B0(i)+B1(i)+B2(i))-(epsilon_db*B0(i)+epsilon_cb*B1(i)))/((1+epsilon_db/1000)*B0(i)+B2(i)+(1+epsilon_cb/1000)*B1(i));
delta_B1(i) = delta_B2(i)*(1+epsilon_cb/1000)+epsilon_cb;
delta_B0(i) = delta_B2(i)*(1+epsilon_db/1000)+epsilon_db;%CO2(aq)
alpha_gb = 1/(epsilon_bg/1000+1);
alpha_gx = exp((-2.4612+(7.6663*1000/TK)-(2.988*1000000/(TK^2)))/1000);
alpha_xb = alpha_gb/alpha_gx;
alpha_cb = epsilon_cb/1000+1;
alpha_bc = 1/alpha_cb;
alpha_xc = alpha_xb/alpha_cb;
Delta_cb(i) = 1000*log(alpha_cb);
Delta_bb(i) = 0;
alpha_C_eq_1 = alpha_xc;                                                    %EFF
alpha_C_eq_2 = alpha_xb;                                                    %EFF
alpha_C_f_1 = 1.000;                                                        %KFF
alpha_C_f_2 = 1.000;                                                        %KFF
alpha_C_b_1 = alpha_C_f_1/alpha_C_eq_1;
alpha_C_b_2 = alpha_C_f_2/alpha_C_eq_2;
%=====================Clumped isotopes=====================================
D47_B0 = 26447/(TK^2)+285.51/TK-0.3004;                                     %Wang et al. (2004) - assume CO2(aq) = CO2(g), as done by Guo (2020)
D48_B0 = 29306/(TK^2)+93.885/TK-0.2914;                                     %Wang et al. (2004) - assume CO2(aq) = CO2(g), as done by Guo (2020)
D49_B0 = 108776/(TK^2)+477.14/TK-0.5954;                                    %Wang et al. (2004) - assume CO2(aq) = CO2(g), as done by Guo (2020)
%D63_B1 = 9.03e-6*TC^2-(3.13e-3)*TC+7.23e-1 - 0.28;                         %Uchikawa et al. (2021) shifted from Hill et al. (2020) - negative AFF to convert to D63
%D63_B2 = 9.39e-6*TC^2-(3.31e-3)*TC+7.91e-1 - 0.28;                         %Uchikawa et al. (2021) shifted from Hill et al. (2020) - negative AFF to convert to D63
D63_B1 = 43187/(TK^2)-34.833/TK+0.0007;                                     %Hill et al. 2020
D63_B2 = 43655/(TK^2)-23.643/TK-0.0088;                                     %Hill et al. 2020
D64_B1 = 23492/(TK^2)-52.842/TK+0.0304;                                     %Hill et al. 2020
D64_B2 = 21842/(TK^2)-50.457/TK+0.0291;                                     %Hill et al. 2020
D65_B1 = 112667/(TK^2)-123.11/TK+0.0304;                                    %Hill et al. 2020
D65_B2 = 112026/(TK^2)-97.208/TK+0.009;                                     %Hill et al. 2020
D63eq = 43159/(TK^2)-25.095/TK-0.0078;                                      %Calcite; Hill et al. 2020
D64eq = 23566/(TK^2)-52.319/TK+0.0297;                                      %Calcite; Hill et al. 2020
D65eq = 112667/(TK^2)-102.28/TK+0.012;                                      %Calcite; Hill et al. 2020
%=====================Isotopologue k's and nu's ===========================
d1 = 3*alpha_eq_c_w*rw;
r_B2(i) = (delta_B2(i)/1000+1)*rDIC;
r_B1(i) = (delta_B1(i)/1000+1)*rDIC;
d2(i) = alpha_C_eq_2*r_B2(1);                                               %Old note: =r_c^eq; replaces alpha_eq_c_w for O isotopes
k1 = kB1;                                                                   %Mass 60    
k2 = kB2;                                                                   %Mass 60
k3 = kB1*alpha_C_f_1;                                                       %Mass 61 
k4 = kB2*alpha_C_f_2;                                                       %Mass 61 
k5 = kB1*alpha_O_f_1;                                                       %Mass 62 
k6 = kB2*alpha_O_f_2;                                                       %Mass 62 
v1 = vB1;                                                                   %Mass 60 
v2 = vB2;                                                                   %Mass 60 
v3 = vB1*alpha_C_b_1;                                                       %Mass 61 
v4 = vB2*alpha_C_b_2;                                                       %Mass 61 
v5 = vB1*alpha_O_b_1;                                                       %Mass 62 
v6 = vB2*alpha_O_b_2;                                                       %Mass 62 
%%
%===================Isotopologue concentrations ===========================
B3(i) = B1(i)*(delta_B1(i)/1000+1)*rDIC;                                    %Mass 61
B4(i) = B2(i)*(delta_B2(i)/1000+1)*rDIC;                                    %Mass 61
B5(i) = 3*alpha_cw*rw*B1(i);                                                %Mass 62
B6(i) = 3*alpha_bw*rw*B2(i);                                                %Mass 62
%===============Carbon isotope composition of calcite======================
AR(i) =  k1*B1(i)*PA(i)+k2*B2(i)*PA(i);
DR(i) = v1*PB1(i)+v2*PB2(i);
r_c(i) = d2(i)*AR(i)/(((d2(i)*(AR(i)-DR(i))*(1+phi))/(alpha_C_f_1*r_B1(i)+phi*alpha_C_f_2*r_B2(i)))+DR(i));
alpha_cB2(i) = r_c(i)/r_B2(i);
delta_c(i) = (r_c(i)/rDIC-1)*1000;
%===============Oxygen isotope composition of calcite======================
R_o(i) = d1*AR(i)/(((d1*(AR(i)-DR(i))*(1+phi))/(alpha_O_f_1*3*alpha_cw*rw+phi*alpha_O_f_2*3*alpha_bw*rw))+DR(i));
alpha_o(i) = R_o(i)/3/rw;
%===============13-18 clumped isotope composition of calcite=====================
epsilon = 0;
alpha_63_f_1 = alpha_O_f_1*alpha_C_f_1+epsilon;
alpha_63_f_2 = alpha_O_f_2*alpha_C_f_2+epsilon;
carbon_isotopes(i) = alpha_xc*(AR(i)-DR(i))/((1-chi)*alpha_C_f_1+chi*alpha_C_f_2*alpha_bc)+DR(i);
oxygen_isotopes(i) = alpha_xw/alpha_cw*(AR(i)-DR(i))/((1-chi)*alpha_O_f_1+chi*alpha_O_f_2*alpha_bw/alpha_cw)+DR(i);
RR63(i) = carbon_isotopes(i)*oxygen_isotopes(i)*(D63eq/1000+1)/(AR(i)*(DR(i)+alpha_xc*alpha_xw/alpha_cw*(D63eq/1000+1)*(AR(i)-DR(i))/((1-chi)*alpha_63_f_1*(D63_B1/1000+1)+chi*alpha_63_f_2*alpha_bc*alpha_bw/alpha_cw*(D63_B2/1000+1))));
Delta63(i) = (RR63(i)-1)*1000;
%===============18-18 clumped composition of calcite======================
epsilon = 0;
alpha_64_f_1 = alpha_O_f_1*alpha_O_f_1+epsilon;
alpha_64_f_2 = alpha_O_f_2*alpha_O_f_2+epsilon;
RR64(i) = oxygen_isotopes(i)*oxygen_isotopes(i)*(D64eq/1000+1)/(AR(i)*(DR(i)+(alpha_xw/alpha_cw)^2*(D64eq/1000+1)*(AR(i)-DR(i))/((1-chi)*alpha_64_f_1*(D64_B1/1000+1)+chi*alpha_64_f_2*(alpha_bw/alpha_cw)^2*(D64_B2/1000+1))));
Delta64(i) = (RR64(i)-1)*1000;
%===============13-18-18 clumped composition of calcite======================
epsilon = 0;
alpha_65_f_1 = alpha_O_f_1*alpha_O_f_1*alpha_C_f_1+epsilon;
alpha_65_f_2 = alpha_O_f_2*alpha_O_f_2*alpha_C_f_2+epsilon;
RR65(i) = carbon_isotopes(i)*oxygen_isotopes(i)^2*(D65eq/1000+1)/(AR(i)^2*(DR(i)+alpha_xc*(alpha_xw/alpha_cw)^2*(D65eq/1000+1)*(AR(i)-DR(i))/((1-chi)*alpha_65_f_1*(D65_B1/1000+1)+chi*alpha_64_f_2*alpha_bc*(alpha_bw/alpha_cw)^2*(D65_B2/1000+1))));
Delta65(i) = (RR65(i)-1)*1000;
%==========================================================================
end
%%
disp(pH);
%For pH plot at specified rate
alpha18_spec = interp1(R_c, 1000*log(alpha_o), Rate_spec);
delta13_spec = interp1(R_c, delta_c, Rate_spec);
D63_spec = interp1(R_c, Delta63, Rate_spec);
D64_spec = interp1(R_c, Delta64, Rate_spec);
D65_spec = interp1(R_c, Delta65, Rate_spec);
%---------------------------Carbon-----------------------------------------
pH_stored(iter)=pH;
Rate_spec_stored(iter) = log10(Rate_spec);
colorcode_stored(iter) = log10(Rate_spec);
D_B2(iter)=1000*log((delta_B2(1)+1000)/(d13C_DIC+1000));
D_B0(iter)=1000*log((delta_B0(1)+1000)/(d13C_DIC+1000));
D_B1(iter)=1000*log((delta_B1(1)+1000)/(d13C_DIC+1000));
D_DIC(iter)=1000*log((d13C_DIC+1000)/(d13C_DIC+1000));
D13C_spec_stored(iter) = 1000*log((delta13_spec+1000)/(d13C_DIC+1000));
D13C_stored_slow(iter) = 1000*log((delta_c(1)+1000)/(d13C_DIC+1000));   
D13C_stored_fast(iter) = 1000*log((delta_c(length(Omega))+1000)/(d13C_DIC+1000));
%---------------------------Oxygen-----------------------------------------
D_CO2_water(iter)=1000*log(alpha_dw);
D_B1_water(iter)=1000*log(alpha_cw);
D_B2_water(iter)=1000*log(alpha_bw);
D_DIC_water(iter)=1000*log(AT);
alpha_stored_fast(iter)=alpha_o(length(Omega));
alpha_stored_slow(iter)=alpha_o(1);
alpha_DIC_stored(iter) = AT;
D18O_stored_slow(iter) = 1000*log(alpha_stored_slow(iter));
D18O_stored_fast(iter) = 1000*log(alpha_stored_fast(iter));
D18O_spec_stored(iter) = (alpha18_spec);
%---------------------------Clumped----------------------------------------
D63_B1_stored(iter) = D63_B1;
D63_B2_stored(iter) = D63_B2;
D64_B1_stored(iter) = D64_B1;
D64_B2_stored(iter) = D64_B2;
D65_B1_stored(iter) = D65_B1;
D65_B2_stored(iter) = D65_B2;
D63_calcite_kinetic(iter) = Delta63(length(Omega));
D63_calcite_equilibrium(iter) = Delta63(1);
D63_spec_stored(iter) = D63_spec;
D64_calcite_equilibrium(iter) = Delta64(1);
D64_spec_stored(iter) = D64_spec;
D65_calcite_equilibrium(iter) = Delta65(1);
D65_spec_stored(iter) = D65_spec;
%===========Subsection for plotting rate dep at different pH values========
% if rem(pH,1) == 0 % then pH is an integer 
% for i = 1:length(Omega)
% Rate_stored(i,niter)=log10(R_c(i));
% oxygen_stored(i,niter) = alpha_o(i);
% carbon_stored(i,niter) = alpha_c(i);
% Delta63_stored(i,niter) = Delta63(i);
% Delta_B1_stored(i) = 1000*log(alpha_cw);
% Delta_B2_stored(i) = 1000*log(alpha_bw);
% Delta_B1_carbon_stored(i) = Delta_cb(1);
% Delta_B2_carbon_stored(i) = Delta_bb(1);
% end
% niter=niter+1;
% else
% end
%==========================================================================
iter=iter+1;
end
%%
%For plotting D63,D64, and D65 of EIC as function of pH'
xpH = (0:0.1:14)';
xH = 10.^-xpH;
for i = 1:length(pH_stored)
xB1(i) = 1e-4;
xDIC(i) = xB1(i)*(1+xH(i)/K2+(xH(i)^2)/(K1*K2));
xB2(i) = xDIC(i)/(1+xH(i)/K1+K2/xH(i));% 
xB0(i) = xDIC(i)/(1+K1/xH(i)+K1*K2/(xH(i)^2)); 
xDelta_EIC(i) = xB1(i)/(xB2(i)+xB1(i))*D_B1_water(i) + xB2(i)/(xB2(i)+xB1(i))*D_B2_water(i); 
xD63_EIC(i) = xB2(i)/(xB1(i)+xB2(i))*D63_B2 + xB1(i)/(xB1(i)+xB2(i))*D63_B1;
xD64_EIC(i) = xB2(i)/(xB1(i)+xB2(i))*D64_B2 + xB1(i)/(xB1(i)+xB2(i))*D64_B1;
xD65_EIC(i) = xB2(i)/(xB1(i)+xB2(i))*D65_B2 + xB1(i)/(xB1(i)+xB2(i))*D65_B1;
end
%%
%==========================Import data=====================================
Romanek_data = load('Romanek_to_5C.txt'); 
    pH_Romanek = Romanek_data(:,1);
    D13C_Romanek = Romanek_data(:,2);
    rate_Romanek = Romanek_data(:,3); 
Watkins_data = load('Watkins_to_5C.txt');
    pH_Watkins = Watkins_data(:,1);
    D18O_Watkins = Watkins_data(:,2);  
    rate_Watkins = Watkins_data(:,3);
Baker_data = load('Baker_to_5C.txt');
    pH_Baker = Baker_data(:,1);
    D18O_Baker = Baker_data(:,2);
    rate_Baker = Baker_data(:,2);
Levitt_data = load('Levitt_to_5C.txt');
    pH_Levitt = Levitt_data(:,1);
    D18O_Levitt = Levitt_data(:,2);
Levitt_data_2 = load('Levitt_to_5C_carbon.txt');
    pH_Levitt_carbon = Levitt_data_2(:,1);
    D13C_Levitt_carbon = Levitt_data_2(:,3);
    rate_Levitt_carbon = Levitt_data_2(:,4);
    
%=============================pH Plot======================================
figure
set(gcf,'DefaultAxesLineWidth',1.2)
set(gcf,'DefaultLineLineWidth',1.2)
set(gcf,'Position',[1500, 1000, 600, 1000])

subplot(5,1,1)
plot(pH_stored, D_B1_water, pH_stored, D_B2_water, xpH, xDelta_EIC, pH_stored, D18O_spec_stored, pH_stored, D18O_stored_slow,pH_Watkins, D18O_Watkins, 'ro',pH_Baker, D18O_Baker, 'bs', pH_Levitt,D18O_Levitt,'kd')
xlabel('pH')
ylabel('1000ln\alpha^{18}O')
legend('HCO_{3}^{-}','CO_{3}^{2-}','EIC','calcite (fast)','calcite(slow)','Watkins et al. (2014)','Baker et al. (2015)','Levitt et al. (2018)','Location','west')
subplot(5,1,2)
plot(pH_stored, D_B1, pH_stored, D_B2, pH_stored, D_DIC, pH_stored, D13C_spec_stored, pH_stored, D13C_stored_slow, pH_Romanek, D13C_Romanek, 'o',pH_Levitt_carbon,D13C_Levitt_carbon,'kd')
xlabel('pH')
ylabel('1000ln\alpha^{13}C')
legend('HCO_{3}^{-}','CO_{3}^{2-}','DIC','calcite(fast)','calcite (slow)','Romanek et al. (1992)','Levitt et al. (2018)','Location','west')
subplot(5,1,3)
hold on
plot(pH_stored,D63_B2_stored, pH_stored,D63_B1_stored, xpH,xD63_EIC,pH_stored,D63_spec_stored, pH_stored,D63_calcite_equilibrium)
xlabel('pH')
ylabel('\Delta_{63}')
legend('HCO_{3}^{-}','CO_{3}^{2-}','EIC','calcite (fast)','calcite(slow)','Location','west')
subplot(5,1,4)
plot(pH_stored,D64_B2_stored, pH_stored,D64_B1_stored, xpH',xD64_EIC, pH_stored,D64_spec_stored, pH_stored,D64_calcite_equilibrium)
xlabel('pH')
ylabel('\Delta_{64}')
legend('HCO_{3}^{-}','CO_{3}^{2-}','EIC','calcite (fast)','calcite(slow)','Location','west')
subplot(5,1,5)
plot(pH_stored,D65_B2_stored, pH_stored,D65_B1_stored,xpH,xD65_EIC,pH_stored,D65_spec_stored, pH_stored,D65_calcite_equilibrium)
xlabel('pH')
ylabel('\Delta_{65}')
legend('HCO_{3}^{-}','CO_{3}^{2-}','EIC','calcite (fast)','calcite(slow)','Location','west')

%% Outputting the isotopic composition of DIC species and calcite as a function of pH
%  file=horzcat(pH_stored,D_B1_water,D_B2_water,D_DIC_water,D18O_stored_slow);
%  save Ion_by_ion_updated_oxygen.txt file -ascii
%  file=horzcat(pH_stored, D_B0, D_B1, D_B2, D_DIC,D13C_stored_slow);
%  save Ion_by_ion_updated_carbon.txt file -ascii
%  file=horzcat(pH_stored,D63_B1_stored, D63_B2_stored, D63_calcite_equilibrium,D64_B1_stored, D64_B2_stored, D64_calcite_equilibrium,D65_B1_stored,D65_B2_stored, D65_calcite_equilibrium);
%  save Ion_by_ion_updated_clumped.txt file -ascii
%  file=horzcat(xpH,xD63_EIC,xD64_EIC,xD65_EIC);
%  save Ion_by_ion_updated_clumped_EIC.txt file -ascii
%  file=horzcat(pH_stored,D18O_spec_stored,D13C_spec_stored,D63_spec_stored,D64_spec_stored,D65_spec_stored,Rate_spec_stored);
%  save Ion_by_ion_updated_CaCO3.txt file -ascii
