clc
clear

%**************************************************************************
% Goals:
% (1) test Multi-Layer Photosyntehsis model
% (2) test Depury and Farquhar Photosynthesis model
%**************************************************************************

%% Step 1--default model parameters
FLAG=1; % model version control; 1--Lloyd et al. 2010 Vc-LAI relationship; 2--Mercado et al. 2006 Vc-LAI relationship;
SZA=30; % Solar Zenith Angle, in degree
Pres=10.^5; % Atmosphere Pressure, in pa
LAI=6; % Leaf Area Index

Tl=28; % leaf temperature for sunlit leaf
Tldiff=0; % leaf temperature difference between sunlit and shade leaves

ambCO2=380; % Ambient CO2 in ppm
Vcmax0_25=40; % Bonan et al., 2012 for the tropcis
CI=0.63; % Clumping index, from Chen etal., 2005 for tropical evergreen forests
Topt=35; % optimal leaf temperature for the tropics, from Lloyd and Farquhar, 2008

N=20; % number of layers for Multi-Layer Canopy Photosynthesis Modeling

% variables related to sun/shade leaf maximum intrinsic quantum yield
Phi_sun=0.7; 
PSII_sun=0.7; 
Phi_shade=0.7;
PSII_shade=0.7; 

% varaibles related to leaf age effect; scale factor 
sf_sun=1; 
sf_shade=1;
sf=1;

%% Step 2--Call Light Partitioning Function
PAR0=1320; % top canop irradiance, in umol/m2/s
LQ=Func_Light_Partitioning(SZA,Pres, PAR0);

%% Step3--Call Multi-Layer Photosynthesis Model
[DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vcmax0_25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf);

%% Step4--Canopy Photosynthesis vs. ML
NL=[1:1:30]'; % number of Layer
for i=1:length(NL)
    N=NL(i);
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vcmax0_25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf);
    An_ML(i,1)=ML.An_tot;
    clear DF ML
end

figure('color','white');
plot(NL,An_ML,'k-o','LineWidth',2);
xlabel('Number of Layer','fontsize',14);
ylabel('Canopy Photosynthesis','fontsize',14);
set(gca,'fontsize',12);


%% Step5--Canopy Photosynthesis vs. PAR
clear DF ML
PARobs=[10:10:2150]'*cos(SZA./180.*pi);
N=20;
for i=1:length(PARobs)
    PAR0=PARobs(i,1);
    LQ=Func_Light_Partitioning(SZA,Pres, PAR0);
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vcmax0_25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf);
    
    An_DF(i,1)=DF.An_tot;
    An_ML(i,1)=ML.An_tot;
        
    clear DF ML LQ PAR0    
end

figure('color','white');
plot(PARobs,An_DF(:,1),'r-','LineWidth',2);
hold on
plot(PARobs,An_ML(:,1),'b-','LineWidth',2);
xlabel('PAR(umol/m2/s)','fontsize',14);
ylabel('Canopy Photosynthesis','fontsize',14);
set(gca,'fontsize',12);


%% Step6--Canopy Photosynthesis vs. PAR on different LAI
clear An_ML
Vcmax0_25=60;
PARobs=[10:10:2150]'*cos(SZA./180.*pi);
N=20;
CI=0.8;
LAI0=[1 2 3 4 5 6]';
for i=1:length(LAI0)
    for j=1:length(PARobs)
        PAR0=PARobs(j,1);
        LQ=Func_Light_Partitioning(SZA,Pres, PAR0);
        LAI=LAI0(i,1);
        [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vcmax0_25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf);
        An_ML(j,i)=ML.An_tot;
        clear DF ML
    end
end

figure('color','white');
plot(PARobs,An_ML(:,1),'-','LineWidth',2,'color',[0 0 0]);
hold on
plot(PARobs,An_ML(:,2),'-','LineWidth',2,'color',[0 0 1]);
plot(PARobs,An_ML(:,3),'-','LineWidth',2,'color',[0 1 0]);
plot(PARobs,An_ML(:,4),'-','LineWidth',2,'color',[0 1 1]);
plot(PARobs,An_ML(:,5),'-','LineWidth',2,'color',[1 0 1]);
plot(PARobs,An_ML(:,6),'-','LineWidth',2,'color',[1 0 0]);
xlabel('PAR(umol/m2/s)','fontsize',14);
ylabel('Canopy Photosynthesis','fontsize',14);
set(gca,'fontsize',12);
