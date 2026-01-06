%% %File: COAD_Behavior.m
%Created: 10/14/2021
%Last modified: 1/10/2022
%Author: James M. Watkins
%Department of Earth Sciences, University of Oregon
%watkins4@uoregon.edu
%Description: This script was used to produce Fig. 5 in Watkins and Devriendt
%(2022). If you use it, please cite: Watkins, J. and Devriendt, L., 2022, A
%combined model for kinetic clumped isotope effects in the
%CaCO3-DIC-H2O system, Geochemistry, Geophysics, Geosystems, v. X, p. Y-Z. 

%% %Instructions
%Box model with two fluxes: (1) FCO2 in and (2) FCaCO3 out. FCaCO3 is calculated from 'CaCO3_DIC.m'
%The post-processing script is 'Run_COAD_Behavior.m'
%To run, put all files in same directory and execute 'Run_COAD_Behavior.m'
%Input variables and parameters can be changed in this file. 

%%
function varargout = COAD_Behavior(varargin)
   [varargout{1:nargout}] = feval(varargin{:});
end

function o=setopts()
    %----------Set parameters----------------------------------------------
    o.TC = 5;
    o.TK = o.TC+273.15;
    o.S = 0;
    o.pH = 10.5;
    o.Ca = 10e-3;
    o.rVSMOW = 2005.2/1e6;
    o.V = 5;
    o.kcatKM = 2.7e7;                                                           %Uchikawa and Zeebe (2012)
    o.CA = 0;                                                                   %Units of moles/L    
    o.Fluxin = 0.01;                                                             %0.025 is close to the data (mmol/h)
    o.FCO2 = o.Fluxin/1000/60/60;                                               %moles/s
    o.FCaCO3 = 0;
    o.totalcarb = 0;
    o.SAcarb = 0.1;                                                             %Representative value. Ranges between 0.05-0.35
    %-----------Millero et al. (2007) - NaCl-------------------------------
    o.mNaCl = 0.0;                                                              %(seawater mNaCl = 0.537) (15 mM from Morrill et al., 2013)
    o.A1 = 35.2911*o.mNaCl^0.5+0.8491*o.mNaCl-0.32*o.mNaCl^1.5+0.055*o.mNaCl^2;
    o.B1 = -1583.09*o.mNaCl^0.5;
    o.C1 = -5.4366*(o.mNaCl^0.5);
    o.pK1o = -114.3106+5773.67/o.TK+17.779524*log(o.TK);
    o.pK1 = o.A1+o.B1/o.TK+o.C1*log(o.TK)+o.pK1o;                               %Millero et al. (2007) - NaCl  
    o.A2 = 38.2746*o.mNaCl^0.5+1.6057*o.mNaCl-0.647*o.mNaCl^1.5+0.113*o.mNaCl^2;
    o.B2 = -1738.16*o.mNaCl^0.5;
    o.C2 = -6.0346*o.mNaCl^0.5;
    o.pK2o = -83.2997+4821.38/o.TK+13.5962*log(o.TK);
    o.pK2 = o.A2+o.B2/o.TK+o.C2*log(o.TK)+o.pK2o;                               %Millero et al. (2007) - NaCl
    o.Ksp = 10^-8.38;                                                           %Jacobsen (1974) calculated by Ellen Olsen
    o.lnKw = 148.96502-13847.26/o.TK-23.6521*log(o.TK)+(118.67/o.TK-5.977+1.0495*log(o.TK))*o.S^0.5-0.01615*o.S;%p258
    o.lnK0 = 9345.17/o.TK-60.2409+23.3585*log(o.TK/100)+o.S*(0.023517-0.00023656*o.TK+0.0047036*(o.TK/100)^2);%p.257
    o.K1 = 10^(-o.pK1);
    o.K2 = 10^(-o.pK2);
    o.Kw = exp(o.lnKw);
    o.K0 = exp(o.lnK0);
    %----------Speciate----------------------------------------------------
    o.H = 10^-o.pH;
    o.OH = o.Kw/o.H;
    o.omega = 7.0;                                                              %Estimated from SIcrit (Dietzel et al. 2009)
    o.CO3 = o.omega*o.Ksp/o.Ca;	
    o.DIC = o.CO3*(1+o.H/o.K2+o.H^2/(o.K1*o.K2));                               %Appendix B of Zeebe and Wolf-Gladrow (2001)
    o.CO2 = o.DIC/(1+o.K1/o.H+o.K1*o.K2/(o.H^2));                               %Eq. 1.1.9 of Zeebe and Wolf-Gladrow (2001)
    o.HCO3 = o.DIC/(1+o.H/o.K1+o.K2/o.H);                                       %Eq. 1.1.10 of Zeebe and Wolf-Gladrow (2001)
    o.EIC = o.HCO3+o.CO3;
    o.chi = 1/(1+o.K2/o.H);
    o.omega = o.Ca*o.CO3/o.Ksp;
    %----------Oxygen isotopes---------------------------------------------
    o.d18Ow = 0;                                                                %Arbitrary
    o.rw = (o.d18Ow/1000+1)*o.rVSMOW;                                           %Definition
    o.alpha_dw = exp(2520/o.TK^2 + 0.01212);                                    %Beck et al. 2005
    o.alpha_bw = exp(2590/o.TK^2 + 0.00189);                                    %Beck et al. 2005
    o.alpha_cw = exp(2390/o.TK^2 - 0.00270);                                    %Beck et al. 2005
    o.alpha_xw = exp((17747/o.TK-29.777)/1000);                                 %Watkins et al. (2014); very similar to Coplen (2007)
    o.eOH = -4.4573+10.3255e3/o.TK -0.5976e6/(o.TK^2);                          %Zeebe (2020)
    o.alpha_OH = 1/(o.eOH/1000+1);                                              %Zeebe (2020)
    o.rOH = o.alpha_OH*o.rw;                                                    %Definition
    o.chi18 = 1/(1+o.K2*o.alpha_cw/o.alpha_bw/o.H);                             %Chen et al. (2018)
    %----------Carbon isotopes---------------------------------------------
    o.d13C_DIC = -10;                                                           %Estimated from Tang experiments
    o.rDIC = (o.d13C_DIC/1000+1)*0.01118;                                       %Definition
    o.epsilon_gb = -9483/o.TK+23.89;                                            %Mook (1986)
    o.epsilon_cb = -867/o.TK+2.52;                                              %Mook (1986)
    o.epsilon_db = -9866/o.TK+24.12;                                            %Mook (1986)
    o.delta_HCO3 = (o.d13C_DIC*o.DIC-(o.epsilon_db*o.CO2+o.epsilon_cb*o.CO3))/((1+o.epsilon_db/1000)*o.CO2+o.HCO3+(1+o.epsilon_cb/1000)*o.CO3);
    o.delta_CO3 = o.delta_HCO3*(1+o.epsilon_cb/1000)+o.epsilon_cb;
    o.delta_CO2 = o.delta_HCO3*(1+o.epsilon_db/1000)+o.epsilon_db;
    o.alpha_gb = o.epsilon_gb/1000+1;
    o.alpha_cb = o.epsilon_cb/1000+1;
    o.alpha_db = o.epsilon_db/1000+1;
    o.alpha_gx = exp((-2.4612+(7.6663*1000/o.TK)-(2.988*1000000/(o.TK^2)))/1000);
    o.alpha_xb = o.alpha_gb/o.alpha_gx+0.002;
    o.alpha_xc = o.alpha_xb/o.alpha_cb+0.002;
    %----------Clumped isotopes--------------------------------------------
    o.D47_CO2eq = 26447/(o.TK^2)+285.51/o.TK-0.3004;                            %Wang et al. (2004) - assume CO2(aq) = CO2(g), as done by Guo (2020)
    %o.D63_HCO3eq = 9.39e-6*o.TC^2-(3.31e-3)*o.TC+7.91e-1 - 0.28;               %Uchikawa et al. (2021) shifted from Hill et al. (2020) - negative AFF to convert to D63
    %o.D63_CO3eq = 9.03e-6*o.TC^2-(3.13e-3)*o.TC+7.23e-1 - 0.28;                %Uchikawa et al. (2021) shifted from Hill et al. (2020) - negative AFF to convert to D63
    o.D63_HCO3eq = 43655/(o.TK^2)-23.643/o.TK-0.0088;                           %Hill et al. 2020
    o.D63_CO3eq = 43187/(o.TK^2)-34.833/o.TK+0.0007;                            %Hill et al. 2020
    o.RR47_CO2 = 1+o.D47_CO2eq/1000;                                            %A58 in Uchikawa (2021)
    o.RR63_HCO3 = 1+o.D63_HCO3eq/1000;                                          %A61 in Uchikawa (2021)
    o.RR63_CO3 = 1+o.D63_CO3eq/1000;                                            %A62 in Uchikawa (2021)
    o.alpha63_cb = (1+o.D63_CO3eq/1000)/(1+o.D63_HCO3eq/1000);                  %A63 in Uchikawa (2021)
    o.K2_63 = o.K2*o.alpha63_cb*o.alpha_cb*o.alpha_cw/o.alpha_bw;               %A80 in Uchikawa (2021)
    %----------Double clumped isotopes-------------------------------------
    o.D48_CO2eq = 29306/(o.TK^2)+93.885/o.TK-0.2914;                            %Wang et al. (2004) - assume CO2(aq) = CO2(g), as done by Guo (2020)
    o.D64_HCO3eq = 21842/(o.TK^2)-50.457/o.TK+0.0291;                           %Hill et al. 2020
    o.D64_CO3eq = 23492/(o.TK^2)-52.842/o.TK+0.0304;                            %Hill et al. 2020
    o.RR48_CO2 = 1+o.D48_CO2eq/1000;                               	
    o.RR64_HCO3 = 1+o.D64_HCO3eq/1000;                             	
    o.RR64_CO3 = 1+o.D64_CO3eq/1000;                               	
    o.alpha64_cb = (1+o.D64_CO3eq/1000)/(1+o.D64_HCO3eq/1000);        	
    o.K2_64 = o.K2*o.alpha64_cb*(o.alpha_cw/o.alpha_bw)^2;            	
    %----------Forward k's-------------------------------------------------
    o.kf1o = 10^(329.85-110.541*log10(o.TK)-(17265.4/o.TK));                    %Uchikawa and Zeebe (2012)
    o.kf1 = o.kf1o + o.kcatKM*o.CA;                                             %Michaelis-Menten kinetics
    o.af1 = o.kf1*1.0000;                                                       %Yumol et al. (2020)
    o.bf1 = o.kf1*0.9812;                                                       %Yumol et al. (2020)
    o.cf1 = o.kf1*0.9824;                                                       %Yumol et al. (2020)
    %o.kf4 = 10^(13.635-2895/o.TK);                                             %Uchikawa and Zeebe (2012)
    o.A4 = 499002.24*exp(4.2986e-4*o.S^2+5.75499e-5*o.S);                       %Schulz et al. (2006)
    o.kf4 = o.A4*exp(-90166.83/(8.3145*o.TK))/o.Kw;                             %Schulz et al. (2006)
    o.af4 = o.kf4*0.978;                                                        %Christensen et al. (2021) - corrected
    o.bf4 = o.kf4*1.000;                                                        %Christensen et al. (2021)
    o.cf4 = o.kf4/1.019;                                                        %Christensen et al. (2021)
    %----------Backward k's------------------------------------------------
    o.kb1 = o.kf1/o.K1;
    o.ab1 = o.af1/(o.alpha_bw*o.K1);
    o.bb1 = o.bf1/(o.alpha_bw/o.alpha_dw*o.K1);
    o.cb1 = o.cf1*o.alpha_db/o.K1;
    o.kb4 = o.kf4*o.Kw/o.K1;
    o.ab4 = o.af4*o.Kw/(o.alpha_bw/o.alpha_OH*o.K1);
    o.bb4 = o.bf4*o.Kw/(o.alpha_bw/o.alpha_dw*o.K1);
    o.cb4 = o.cf4*o.Kw*o.alpha_db/o.K1;
    %-------Single clumped-------------------------------------------------
    o.KIE_p1 = (4613.8393/o.TK^2-4.0389/o.TK-0.185)/1000+1;                     %Guo (2020) 
    o.KIE_p4 = (-902.7635/o.TK^2+157.1718/o.TK-0.533)/1000+1;                   %Guo (2020)
    o.KIE_s1 = (-5705.688/o.TK^2-41.5925/o.TK-0.015)/1000+1;                    %Guo (2020)
    o.KIE_s4 = (-11771.2832/o.TK^2-62.7060/o.TK+.168)/1000+1;                   %Guo (2020)
    o.pf1 = o.KIE_p1*o.cf1*o.af1/o.kf1;                                         %Table 4
    o.sf1 = o.KIE_s1*o.cf1*o.bf1/o.kf1;                                         %Table 4
    o.pf4 = o.KIE_p4*o.cf4*o.af4/o.kf4;                                         %Table 4
    o.sf4 = o.KIE_s4*o.cf4*o.bf4/o.kf4;                                         %Table 4
    o.pb1 = o.pf1*o.alpha_db/(o.K1*o.RR63_HCO3*o.alpha_bw);                     %Table 4                   
    o.sb1 = o.sf1*o.RR47_CO2*o.alpha_db/(o.K1*o.RR63_HCO3*o.alpha_bw/o.alpha_dw);%Table 4     
    o.pb4 = o.pf4*o.Kw*o.alpha_db*o.alpha_OH/(o.K1*o.RR63_HCO3*o.alpha_bw);     %Table 4         
    o.sb4 = o.sf4*o.Kw*o.RR47_CO2*o.alpha_db/(o.K1*o.RR63_HCO3*o.alpha_bw/o.alpha_dw);%Table 4 
    %======Double clumped======================================================
    o.KIE_p1p = (13249.5324/o.TK^2-37.8964/o.TK+0.027)/1000+1;                  %Guo (2020) 
    o.KIE_p4p = (5859.1625/o.TK^2-3.7964/o.TK-0.197)/1000+1;                    %Guo (2020)
    o.KIE_s1p = (-18411.4121/o.TK^2-3.7575/o.TK+0.074)/1000+1;                  %Guo (2020)
    o.KIE_s4p = (-12333.8137/o.TK^2+8.6005/o.TK+0.024)/1000+1;                  %Guo (2020)
    o.pf1p = o.KIE_p1p*o.bf1*o.af1/o.kf1;                                       %Table 5
    o.sf1p = o.KIE_s1p*o.bf1*o.bf1/o.kf1;                                       %Table 5
    o.pf4p = o.KIE_p4p*o.bf4*o.af4/o.kf4;                                       %Table 5
    o.sf4p = o.KIE_s4p*o.bf4*o.bf4/o.kf4;                                       %Table 5
    o.pb1p = o.pf1p/(o.K1*o.RR64_HCO3*(o.alpha_bw/o.alpha_dw)*o.alpha_bw);      %Table 5               
    o.sb1p = o.sf1p*o.RR48_CO2/(o.K1*o.RR64_HCO3*(o.alpha_bw/o.alpha_dw)^2);    %Table 5             
    o.pb4p = o.pf4p/(o.K1/o.Kw*o.RR64_HCO3*(o.alpha_bw/o.alpha_dw)*o.alpha_bw/o.alpha_OH);%Table 5  
    o.sb4p = o.sf4p*o.RR48_CO2/(o.K1/o.Kw*o.RR64_HCO3*(o.alpha_bw/o.alpha_dw)^2);%Table 5             
    %----------chi's-------------------------------------------------------
    o.chi18 = 1/(1+o.K2*o.alpha_cw/o.alpha_bw/o.H);                             %Table 3
    o.chi13 = 1/(1+o.K2*o.alpha_cb/o.H);                                        %Table 3
    o.chi63 = 1/(1+o.K2_63/o.H);                                                %Table 4
    o.chi64 = 1/(1+o.K2_64/o.H);                                                %Table 5
    %----------Isotope ratios----------------------------------------------
    o.oCO2 = o.alpha_dw*o.rw;
    o.oHCO3 = o.alpha_bw*o.rw;
    o.oCO3 = o.alpha_cw*o.rw;
    o.cCO2 = (o.delta_CO2/1000+1)*0.01118;
    o.cHCO3 = (o.delta_HCO3/1000+1)*0.01118;
    o.cCO3 = (o.delta_CO3/1000+1)*0.01118;
    %----------Isotopologue concentrations---------------------------------
    o.C6OH = o.OH;
    o.C8OH = o.rOH*o.C6OH;
    o.C266 = o.CO2;
    o.C268 = 2*o.oCO2*o.C266;
    o.C366 = o.cCO2*o.C266;
    o.B2666 = o.HCO3;
    o.B2668 = 3*o.oHCO3*o.B2666;
    o.B3666 = o.cHCO3*o.B2666;
    o.C2666 = o.CO3;
    o.C2668 = 3*o.oCO3*o.C2666;
    o.C3666 = o.cCO3*o.C2666;  
    o.E2666 = o.B2666+o.C2666;
    o.E2668 = o.B2668+o.C2668;
    o.E3666 = o.B3666+o.C3666;
    o.C386 = o.RR47_CO2*o.C366*o.C268/o.C266;                
    o.B3866 = o.RR63_HCO3*o.B3666*o.B2668/o.B2666;  
    o.C3866 = o.RR63_CO3*o.C3666*o.C2668/o.C2666;          
    o.E3866 = o.B3866/o.chi63;     
    o.C288 = 1/4*o.RR48_CO2*o.C268*o.C268/o.C266;
    o.B2886 = 1/3*o.RR64_HCO3*o.B2668*o.B2668/o.B2666;
    o.C2886 = 1/3*o.RR64_CO3*o.C2668*o.C2668/o.C2666;
    o.E2886 = o.B2886/o.chi64;
    %----------CO2 composition---------------------------------------------
    o.d13CCO2in = -17;                                                          %Approx from Tang data
    o.D47CO2in = 26447/(o.TK^2)+285.51/o.TK-0.3004;                             %Assume CO2 in is in clumped isotope equilibrium
    o.D48CO2in = 29306/(o.TK^2)+93.885/o.TK-0.2914;                             %Assume CO2 in is in clumped isotope equilibrium
    o.r18CO2in = o.oCO2;%(o.d18OCO2in/1000+1)*o.rVSMOW;%                        %Case A
%     o.d18OCO2in = 32;                                                         %Case C (45.7 is the equilibrium value) - lower d18O of CO2 due to diffusion?
%     o.r18CO2in = (o.d18OCO2in/1000+1)*o.rVSMOW;%                              %Case C
%     o.D47CO2in = 26447/(o.TK^2)+285.51/o.TK-0.3004+0.10;                      %Case C - the '+0.10' represents higher D47 of CO2 due to diffusion
    o.r13CO2in = (o.d13CCO2in/1000+1)*0.01118;
    o.R47CO2in = (o.D47CO2in/1000+1)*o.r13CO2in*2*o.r18CO2in;
    o.R48CO2in = (o.D48CO2in/1000+1)*o.r18CO2in^2;
end

function [o, time, out]=doit(~)
    close all
    o=setopts();
    [time, out]=driver(o);
end
function [time, out]=driver(o)
    t=[1e-8 10*60*60];     % time (s)
    C2660 = o.C266;         %1 = 12C16O16O
    C2680 = o.C268;         %2 = 12C16O18O
    E26660 = o.E2666;       %3 = E12C16O16O16O
    E26680 = o.E2668;       %4 = E12C16O16O18O
    Ca0 = o.Ca;             %5 = Ca2+
    Carb0 = 0;              %6 = total carbonate    
    C3660 = o.C366;         %7 = 13C16O16O
    E36660 = o.E3666;       %8 = E13C16O16O16O
    C3860 = o.C386;         %9 = clumped CO2
    E38660 = o.E3866;       %10 = clumped EIC
    C2880 = o.C288;         %11 = double clumped CO2
    E28860 = o.E2886;       %12 = double clumped EIC
    
    y0=[C2660,C2680,E26660,E26680,Ca0,Carb0,C3660,E36660,C3860,E38660,C2880,E28860]; %Initial conditions vector
    options = odeset('MaxOrder', 4, 'MaxStep', 1e6, 'RelTol', 10^-9, 'AbsTol',10^-9);
    [time, out]=ode23s(@(t,y0) carbonate(t, y0, o.TC,o.TK,o.pH,o.FCO2,o.kf1,o.af1,o.bf1,o.cf1,o.pf1,o.sf1,o.pf1p,o.sf1p,o.kf4,o.af4,o.bf4,o.cf4,o.pf4,o.sf4,o.pf4p,o.sf4p,o.kb1,o.ab1,o.bb1,o.cb1,o.pb1,o.sb1,o.pb1p,o.sb1p,o.kb4,o.ab4,o.bb4,o.cb4,o.pb4,o.sb4,o.pb4p,o.sb4p,o.H,o.C6OH,o.C8OH,o.rw,o.Ksp,o.chi,o.chi13,o.chi18,o.chi63,o.chi64,o.r18CO2in,o.r13CO2in,o.R47CO2in,o.R48CO2in,o.V,o.SAcarb), t, y0, options);
end

function dy=carbonate(t,y,TC,TK,pH,FCO2,kf1,af1,bf1,cf1,pf1,sf1,pf1p,sf1p,kf4,af4,bf4,cf4,pf4,sf4,pf4p,sf4p,kb1,ab1,bb1,cb1,pb1,sb1,pb1p,sb1p,kb4,ab4,bb4,cb4,pb4,sb4,pb4p,sb4p,H,C6OH,C8OH,rw,Ksp,chi,chi13,chi18,chi63,chi64,r18CO2in,r13CO2in,R47CO2in,R48CO2in,V,SAcarb);
    dy=zeros(size(y));
    DIC = y(1)+y(3);
    omega = y(3)*(1-chi)*y(5)/Ksp;                                              %CO3 = y(3)(1-chi)
        [alpha_c,alpha_o,alpha63,alpha64,alpha65,R_c]=CaCO3_DIC(TC, Ksp, pH, y(5), y(1), y(3)*chi, y(3)*(1-chi));
        alpha_63EIC = alpha63;
        alpha_64EIC = alpha64;
        alpha_65EIC = alpha65;
        alpha_cEIC = alpha_c*((1-chi13)/(1-chi));
        alpha_oEIC = alpha_o*((1-chi18)/(1-chi));
    if omega < 1
          FCaCO3 = 0;
    else
          FCaCO3 = SAcarb*R_c;
    end
%Validation notes: (1) if FCO2 = 0 and omega = 1, then the output is the
%equilibrium value with no time-dependence. 
    dy(1) = -kf1*y(1)+kb1*y(3)*chi*H - kf4*y(1)*C6OH+kb4*y(3)*chi + FCO2/V;   
    dy(2) = -bf1*y(2)+2/3*bb1*y(4)*chi18*H - bf4*y(2)*C6OH+2/3*bb4*y(4)*chi18 + FCO2/V*2*r18CO2in;    
    dy(3) =  kf1*y(1)-kb1*y(3)*chi*H + kf4*y(1)*C6OH-kb4*y(3)*chi -FCaCO3/V;        
    dy(4) =  af1*y(1)*rw-1/3*ab1*y(4)*chi18*H + bf1*y(2)-2/3*bb1*y(4)*chi18*H + af4*y(1)*C8OH-1/3*ab4*y(4)*chi18 + bf4*y(2)*C6OH-2/3*bb4*y(4)*chi18 - FCaCO3/V*y(4)/y(3)*alpha_oEIC;
    dy(5) = 0;%-FCaCO3/V; This is necessary for long runs to avoid [Ca2+] --> 0
    dy(6) = FCaCO3;
    dy(7) = -cf1*y(7)+cb1*y(8)*chi13*H - cf4*y(7)*C6OH+cb4*y(8)*chi13 + FCO2/V*r13CO2in;
    dy(8) = cf1*y(7)-cb1*y(8)*chi13*H + cf4*y(7)*C6OH-cb4*y(8)*chi13 - FCaCO3/V*y(8)/y(3)*alpha_cEIC;
    dy(9) = -sf1*y(9)+2/3*sb1*y(10)*chi63*H - sf4*y(9)*C6OH+2/3*sb4*y(10)*chi63 + FCO2/V*R47CO2in;
    dy(10) = sf1*y(9)-2/3*sb1*y(10)*chi63*H + sf4*y(9)*C6OH-2/3*sb4*y(10)*chi63 + pf1*y(7)*rw-1/3*pb1*y(10)*chi63*H + pf4*y(7)*C8OH-1/3*pb4*y(10)*chi63 - FCaCO3/V*y(10)/y(3)*alpha_63EIC*alpha_cEIC*alpha_oEIC; 
    dy(11) = -sf1p*y(11)+1/3*sb1p*y(12)*chi64*H - sf4p*y(11)*C6OH+1/3*sb4p*y(12)*chi64 + FCO2/V*R48CO2in;
    dy(12) =  sf1p*y(11)-1/3*sb1p*y(12)*chi64*H + sf4p*y(11)*C6OH-1/3*sb4p*y(12)*chi64 + pf1p*y(2)*rw-2/3*pb1p*y(12)*chi64*H + pf4p*y(2)*C8OH-2/3*pb4p*y(12)*chi64 - FCaCO3/V*y(12)/y(3)*alpha_64EIC*alpha_oEIC*alpha_oEIC;  
end

