clc
clear

%**************************************************************************
% Goals of this code: 
% (1) test leaf photosynthesis response to light
% (2) test leaf photosynthesis response to temperature
% (3) test leaf photosynthesis response to PSII(maximum intrinsic quantum yield) and Phi (curvature factor)
%**************************************************************************

%% Step 1: default model parameters
V25=40; % Kattge and Knorr, 2007 for tropical evergreen forests
J25=1.67*V25; % Melyn et al. 2002; Bonan et al. 2012
PAR=2000; % umol/m2/s
Tl=28; % Leaf temperature
Topt=35; % Optimal temperature, based on Lloyd and Farquhar, 2008  
AmbCO2=380; % umol/mol
Ci=AmbCO2*0.7; % Inter-celluar CO2 concentration
Press=10.^5; % atmosphere pressure Pa
PII_in=0.7; % Maximum quantumn yield: 0.7--Bernacchi et al. 2003; 0.85--von Cammera 2000
Phi_in=0.9; % Theta (curvature factor): 0.7--dePury & Farquhar, 1997; Bonan et al., 2012; 0.9--Bernacchi et al. 2003 and Mdelyn et al. 2002
FvCB=Func_Leaf_FvCB_Photosynthesis_Model(V25, J25, Tl, Topt, PAR, Ci, Press, PII_in, Phi_in);

%% Step 2: test 1--photosynthesis vs. light
PAR0=[10:10:2000]';
Vc=[100:-10:10]';
for i=1:length(PAR0)
    PAR=PAR0(i,1);
    for j=1:length(Vc)
        V25=Vc(j,1);
        J25=1.67*V25; % Melyn et al. 2002; Bonan et al. 2012
        FvCB=Func_Leaf_FvCB_Photosynthesis_Model(V25, J25, Tl, Topt, PAR, Ci, Press, PII_in, Phi_in);
        An(i,j)=FvCB.An;
        clear FvCB
    end
end

figure('color','white');
color_index=[1 0 0
    0.9 0.2 0.2
    0.7 0.5 0.2
    0.5 0.7 0.2
    0.2 0.9 0.2
    0.2 0.7 0.5
    0.2 0.5 0.7
    0.2 0.2 0.9
    0    0   1
    0.7 0.2 0.5
    0.5 0.2 0.7
    0.2 0.2 0.9
    0   1    0
    ];
for i=1:length(Vc)
    plot(PAR0,An(:,i),'-','LineWidth',1,'color',color_index(i,:));
    hold on
end
xlabel('PAR(umol/m2/s)','fontsize',14);
ylabel('An(umol CO2/m2/s)','fontsize',14);
set(gca,'fontsize',12);

%% Step 3: test 2--photosynthesis vs. temperature
PAR=2000;
V25=60; %Bonan et al. 2012
J25=1.67*V25; % Melyn et al. 2002; Bonan et al. 2012
T=[10:1:40];

clear An PAR0
for i=1:length(T)
    Tl=T(i);
    FvCB=Func_Leaf_FvCB_Photosynthesis_Model(V25, J25, Tl, Topt, PAR, Ci, Press, PII_in, Phi_in);
    An(i,1)=FvCB.An;
    Rd(i,1)=FvCB.Rd;
    Vcmax(i,1)=FvCB.Vcmax;
    Jmax(i,1)=FvCB.Jmax;
    Gama_star(i,1)=FvCB.Tau_star;
    Kc(i,1)=FvCB.Kc;
    Ko(i,1)=FvCB.Ko;
    clear FvCB
end

figure('color','white');
subplot(2,4,1);
plot(T,An,'k-','LineWidth',2);
xlabel('Leaf Temperature','fontsize',14);
ylabel('An (umol CO2/m2/s)','fontsize',14);
set(gca,'fontsize',12);

subplot(2,4,2);
plot(T,Rd,'k-','LineWidth',2);
xlabel('Leaf Temperature','fontsize',14);
ylabel('Rd (umol CO2/m2/s)','fontsize',14);
set(gca,'fontsize',12);

subplot(2,4,3);
plot(T,Vcmax,'k-','LineWidth',2);
xlabel('Leaf Temperature','fontsize',14);
ylabel('Vcmax (umol CO2/m2/s)','fontsize',14);
set(gca,'fontsize',12);

subplot(2,4,4);
plot(T,Jmax,'k-','LineWidth',2);
xlabel('Leaf Temperature','fontsize',14);
ylabel('Jmax (umol CO2/m2/s)','fontsize',14);
set(gca,'fontsize',12);

subplot(2,4,5);
plot(T,Gama_star,'k-','LineWidth',2);
xlabel('Leaf Temperature','fontsize',14);
ylabel('Gama-star','fontsize',14);
set(gca,'fontsize',12);

subplot(2,4,6);
plot(T,Kc,'k-','LineWidth',2);
xlabel('Leaf Temperature','fontsize',14);
ylabel('Kc','fontsize',14);
set(gca,'fontsize',12);

subplot(2,4,7);
plot(T,Ko,'k-','LineWidth',2);
xlabel('Leaf Temperature','fontsize',14);
ylabel('Ko','fontsize',14);
set(gca,'fontsize',12);


% %% Step 4: test 3--photosynthesis vs. Phi/PSII
% clear Kc Ko Gamma_star Vcmax Jmax Rd An
% clear T
% 
% PAR=2000;
% V25=60; %Bonan et al. 2012
% J25=1.64*V25+29.1; % Bonan et al. 2012
% Tl=30;
% 
% Phi=[0.1:0.1:1]';
% PSII=[0.1:0.1:1]';
% 
% PSII_in=0.7;
% for i=1:length(Phi)
%     Phi_in=Phi(i,1);
%     FvCB=Func_Leaf_FvCB_Photosynthesis_Model(V25, J25, Tl, Topt, PAR, Ci, Press, PII_in, Phi_in);
%     An(i,1)=FvCB.An;
%     clear FvCB
% end
% 
% PSII_in=0.3;
% for i=1:length(Phi)
%     Phi_in=Phi(i,1);
%     FvCB=Func_Leaf_FvCB_Photosynthesis_Model(V25, J25, Tl, Topt, PAR, Ci, Press, PII_in, Phi_in);
%     An(i,2)=FvCB.An;
%     clear FvCB
% end
% 
% Phi_in=0.7;
% for i=1:length(PSII)
%     PSII_in=PSII(i,1);
%     FvCB=Func_Leaf_FvCB_Photosynthesis_Model(V25, J25, Tl, Topt, PAR, Ci, Press, PII_in, Phi_in);
%     An(i,3)=FvCB.An;
%     clear FvCB
% end
% 
% Phi_in=0.3;
% for i=1:length(PSII)
%     PSII_in=PSII(i,1);
%     FvCB=Func_Leaf_FvCB_Photosynthesis_Model(V25, J25, Tl, Topt, PAR, Ci, Press, PII_in, Phi_in);
%     An(i,4)=FvCB.An;
%     clear FvCB
% end
