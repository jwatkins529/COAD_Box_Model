clear all

pH=(8:0.05:12.5)'; 
xpH = zeros(length(pH),1);
xd13C_HCO3 = zeros(length(pH),1);
xd13C_CO3 = zeros(length(pH),1);
xd13C_CaCO3 = zeros(length(pH),1);
xD18O_HCO3 = zeros(length(pH),1);
xD18O_CO3 = zeros(length(pH),1);
xD18O_CaCO3 = zeros(length(pH),1);
xD63_HCO3 = zeros(length(pH),1);
xD63_CO3 = zeros(length(pH),1);
xD63_EIC = zeros(length(pH),1);
xD63_CaCO3 = zeros(length(pH),1);
xD64_CaCO3 = zeros(length(pH),1);
xomega = zeros(length(pH),1);
xRc = zeros(length(pH),1);
xalpha_c = zeros(length(pH),1);

for ii = 1:length(pH)
[o,time,out]=COAD_Box_Model('doit',pH(ii));

dt = zeros(length(out(:,1)),1);
alpha_o = zeros(length(out(:,1)),1);
alpha_c = zeros(length(out(:,1)),1);
alpha63 = zeros(length(out(:,1)),1);
alpha64 = zeros(length(out(:,1)),1);
R_c = zeros(length(out(:,1)),1);
SAcarb = zeros(length(out(:,1)),1);
FCO2 = zeros(length(out(:,1)),1);
Fluxin = zeros(length(out(:,1)),1);
DIC = zeros(length(out(:,1)),1);
H = zeros(length(out(:,1)),1);
chi = zeros(length(out(:,1)),1);
omega = zeros(length(out(:,1)),1);
Ca = zeros(length(out(:,1)),1);
FCaCO3 = zeros(length(out(:,1)),1);
alpha_cEIC = zeros(length(out(:,1)),1);
alpha_oEIC = zeros(length(out(:,1)),1);
R_oCaCO3 = zeros(length(out(:,1)),1);
R_cCaCO3 = zeros(length(out(:,1)),1);
C266 = zeros(length(out(:,1)),1);
C268 = zeros(length(out(:,1)),1);
C366 = zeros(length(out(:,1)),1);
C386 = zeros(length(out(:,1)),1);
C288 = zeros(length(out(:,1)),1);
B2666 = zeros(length(out(:,1)),1); 
B2668 = zeros(length(out(:,1)),1);
B3666 = zeros(length(out(:,1)),1);
B3866 = zeros(length(out(:,1)),1);
B2886 = zeros(length(out(:,1)),1);
C2666 = zeros(length(out(:,1)),1); 
C2668 = zeros(length(out(:,1)),1);
C3666 = zeros(length(out(:,1)),1);
C3866 = zeros(length(out(:,1)),1);
C2886 = zeros(length(out(:,1)),1);
E2666 = zeros(length(out(:,1)),1); 
E2668 = zeros(length(out(:,1)),1);
E3666 = zeros(length(out(:,1)),1);
E3866 = zeros(length(out(:,1)),1);
E2886 = zeros(length(out(:,1)),1);
d18O_CO2 = zeros(length(out(:,1)),1);
d18O_HCO3 = zeros(length(out(:,1)),1);
d18O_CO3 = zeros(length(out(:,1)),1);
d18O_EIC = zeros(length(out(:,1)),1);
d18O_CaCO3 = zeros(length(out(:,1)),1);
D18O_HCO3 = zeros(length(out(:,1)),1);
D18O_CO3 = zeros(length(out(:,1)),1);
D18O_EIC = zeros(length(out(:,1)),1);
D18OCaCO3 = zeros(length(out(:,1)),1);
d13C_CO2 = zeros(length(out(:,1)),1);
d13C_HCO3 = zeros(length(out(:,1)),1);
d13C_CO3 = zeros(length(out(:,1)),1);
d13C_EIC = zeros(length(out(:,1)),1);
d13C_CaCO3 = zeros(length(out(:,1)),1);
D63_HCO3 = zeros(length(out(:,1)),1);
D63_CO3 = zeros(length(out(:,1)),1);
D63_EIC = zeros(length(out(:,1)),1);
D63_CaCO3 = zeros(length(out(:,1)),1);
D64_HCO3 = zeros(length(out(:,1)),1);
D64_CO3 = zeros(length(out(:,1)),1);
D64_EIC = zeros(length(out(:,1)),1);
D64_CaCO3 = zeros(length(out(:,1)),1);

for i=1:length(out(:,1))
    [alpha_c(i),alpha_o(i),alpha63(i),alpha64(i),~,R_c(i)]=CaCO3_DIC(o.TC,o.Ksp,o.pH, out(i,5), out(i,1), out(i,3)*o.chi, out(i,3)*(1-o.chi));
end
for i = 1:length(out(:,1))
    totalcarb=out(:,6);
    SAcarb(i) = o.SAcarb;     
    DIC(i) = out(i,1)+out(i,3);
    Ca(i) = out(i,5);
    H(i) = 10^-o.pH;
    C266(i) = out(i,1);
    C268(i) = out(i,2);
    C366(i) = out(i,7);
    C386(i) = out(i,9); 
    C288(i) = out(i,11);    
    E2666(i) = out(i,3); 
    E2668(i) = out(i,4);
    E3666(i) = out(i,8);
    E3866(i) = out(i,10);
    E2886(i) = out(i,12);
    B2666(i) = out(i,3)*o.chi;
    B2668(i) = out(i,4)*o.chi18;
    B3666(i) = out(i,8)*o.chi13;
    B3866(i) = out(i,10)/(1+o.K2_63/H(i));
    B2886(i) = out(i,12)/(1+o.K2_64/H(i));
    C2666(i) = out(i,3)*(1-o.chi); 
    C2668(i) = out(i,4)*(1-o.chi18);
    C3666(i) = out(i,8)*(1-o.chi13);
    C3866(i) = out(i,10)-B3866(i);
    C2886(i) = out(i,12)-B2886(i);
    d13C_HCO3(i) =(B3666(i)/B2666(i)/0.01118-1)*1000; 
    d13C_CO3(i) = (C3666(i)/C2666(i)/0.01118-1)*1000; 
    d13C_EIC(i) = (E3666(i)/E2666(i)/0.01118-1)*1000;
    d18O_CO2(i) = ((0.5*out(i,2)/out(i,1))/o.rVSMOW-1)*1000;                  
    d18O_HCO3(i) = ((1/3*B2668(i)/B2666(i))/o.rVSMOW-1)*1000;               
    d18O_CO3(i) =  ((1/3*C2668(i)/C2666(i))/o.rVSMOW-1)*1000;               
    d18O_EIC(i) = ((1/3*E2668(i)/E2666(i))/o.rVSMOW-1)*1000;                
    D18O_HCO3(i) = 1000*log((d18O_HCO3(i)+1000)/(o.d18Ow+1000));
    D18O_CO3(i) = 1000*log((d18O_CO3(i)+1000)/(o.d18Ow+1000));
    D63_HCO3(i) =(B3866(i)*B2666(i)/(B3666(i)*B2668(i))-1)*1000;            %63K = R/R*
    D63_CO3(i) = (C3866(i)*C2666(i)/(C3666(i)*C2668(i))-1)*1000;            %63K = R/R*
    D63_EIC(i) = (E3866(i)*E2666(i)/(E3666(i)*E2668(i))-1)*1000;            %63K = R/R*
    D64_HCO3(i) =(3*B2886(i)*B2666(i)/(B2668(i)*B2668(i))-1)*1000;          %64K = 3*R/R*
    D64_CO3(i) = (3*C2886(i)*C2666(i)/(C2668(i)*C2668(i))-1)*1000;          %64K = 3*R/R*
    D64_EIC(i) = (3*E2886(i)*E2666(i)/(E2668(i)*E2668(i))-1)*1000;          %64K = 3*R/R*
    omega(i) = Ca(i)*C2666(i)/o.Ksp; 
    FCO2(i) = o.FCO2;
    Fluxin(i) = o.Fluxin;

    if omega(i) < 1
        FCaCO3(i) = 0;
        else
        FCaCO3(i) = SAcarb(i)*R_c(i)*1000*60*60;                            %mmoles/h
    end
    alpha_oEIC(i) = alpha_o(i)*((1-o.chi18)/(1-o.chi));
    R_oCaCO3(i) = alpha_oEIC(i)*(1/3*E2668(i)/E2666(i));
    R_cCaCO3(i) = alpha_c(i)*C3666(i)/C2666(i);
    d18O_CaCO3(i) = ((R_oCaCO3(i)/o.rVSMOW-1)*1000);                 
    d13C_CaCO3(i) = (R_cCaCO3(i)/0.01118-1)*1000;   
    D63_CaCO3(i) = (alpha63(i)*(D63_EIC(i)/1000+1)-1)*1000;
    D64_CaCO3(i) = (alpha64(i)*(D64_EIC(i)/1000+1)-1)*1000;
end
for i = 2:length(out(:,1))
     dt(i) = time(i)-time(i-1);
end

%------------Store the steady state values---------------------------------
xpH(ii) = pH(ii);
xalpha_o(ii) = alpha_o(end);
xd13C_HCO3(ii) = d13C_HCO3(end);
xd13C_CO3(ii) = d13C_CO3(end);
xd13C_CaCO3(ii) = d13C_CaCO3(end);
xD18O_HCO3(ii) = D18O_HCO3(end);
xD18O_CO3(ii) = D18O_CO3(end);
xD18O_CaCO3(ii) = 1000*log((d18O_CaCO3(end)+1000)/(o.d18Ow+1000));
xD63_HCO3(ii) = D63_HCO3(end);
xD63_CO3(ii) = D63_CO3(end);
xD63_EIC(ii) = D63_EIC(end);
xD63_CaCO3(ii) = D63_CaCO3(end);
xD64_CaCO3(ii) = D64_CaCO3(end);
xomega(ii) = omega(end);
xRc(ii) = R_c(end);	
disp(pH(ii));
end
%% Load data
data = load('Tang_5C.txt'); 
      pH_Tang = data(:,1);
      lna_Tang = data(:,2);                                                 %1000ln(18)alpha
      D47_Tang = data(:,3);                                                 %Absolute reference frame (AFF = 0.28; Tripati et al., 2015)
      rate_Tang = data(:,4);                                                %log10R (moles/m2/s)
%% Make plots
AFF = 0.280;                                                                %Absolute reference frame
AFF2 = 0.131;                                                               %Lucarelli pers. comm.

figure(1)
set(gcf,'DefaultAxesLineWidth',2)
set(gcf,'DefaultLineLineWidth',2)
set(gcf,'Position',[500, 500, 500, 1000])                                   %[a, b, L, W] (a,b) is the lower left corner. L is width. W is height
 
% subplot(4,1,1)
% plot(xpH,xd13C_CaCO3,xpH,xd13C_HCO3,xpH,xd13C_CO3,'k-')
% xlabel('pH')
% ylabel('d13C')
% legend('CaCO3 (steady state)','HCO3','CO3','Location','southwest')

subplot(2,1,1)
plot(xpH,xD18O_CaCO3,'k-',xpH,xD18O_HCO3,'b-',xpH,xD18O_CO3,'r-',pH_Tang,lna_Tang,'bo')
xlabel('pH')
ylabel('1000ln\alpha')
legend('CaCO3 (steady state)','HCO3','CO3','Tang et al. (2014)','Location','southwest')
axis([8 12.5 5 35])

subplot(2,1,2)
plot(xpH,xD63_CaCO3,'k-',pH_Tang,D47_Tang-AFF,'bo')
xlabel('pH')
ylabel('\Delta_{63}')
legend('Model (steady state)','Tang et al. (2014)', 'Location','northwest')

% subplot(4,1,4)
% plot(xpH,xD64_CaCO3+AFF2,'k-')
% xlabel('pH')
% ylabel('\Delta_{48}')
% legend('Model (steady state)', 'Location','east')

%%  
% file=horzcat(xpH,xD18O_HCO3,xD18O_CO3,xD18O_CaCO3,xD63_HCO3+AFF,xD63_CO3+AFF,xD63_CaCO3+AFF,log10(xRc),xd13C_HCO3,xd13C_CO3,xd13C_CaCO3,xD64_CaCO3+AFF2);
% save Guo_Tang_0p003_lowSA_bottcher.txt file -ascii
% copyfile ('Guo_Tang_0p003_lowSA_bottcher.txt','/Users/James/Desktop/GMT/GMT_Guo')