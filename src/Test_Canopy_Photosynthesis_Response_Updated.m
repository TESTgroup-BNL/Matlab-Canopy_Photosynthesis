clc
clear

%**************************************************************************
% Goals:
% ? (1) Modeled GPP against Number of Layers (for multi-layer model), under sunlit, medium diffuse and high diffuse light
% ? (2) Modeled GPP-PAR response
% ? (3) Modeled GPP against temperature response under high, medium, low light 
% Jin Wu, BNL, jinwu@bnl.gov, Jul-2016
%**************************************************************************

%% Step 1--Model parameters
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

N=15; % number of layers for Multi-Layer Canopy Photosynthesis Modeling

% variables related to sun/shade leaf maximum intrinsic quantum yield
Phi_sun=0.7; 
PSII_sun=0.7; 
Phi_shade=0.7;
PSII_shade=0.7; 

% varaibles related to leaf age effect; scale factor 
sf_sun=1; 
sf_shade=1;
sf=1;

LAI_cut=2;
%% Step 2--Test the number of layer on Modeled GPP (for Multi-layer purpose)
% Low light level
PAR0_1=250; % top canop irradiance, in umol/m2/s
LQ_1=Func_Light_Partitioning(SZA,Pres, PAR0_1);
NL=[1:1:30];
for i=1:length(NL)
    N=NL(i);
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ_1, LAI, Tl, Tldiff, ambCO2, Vcmax0_25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);
    An_ML_1(i,1)=ML.An_tot;
    clear DF ML
end

% Medium light level
PAR0_2=1000; % top canop irradiance, in umol/m2/s
LQ_2=Func_Light_Partitioning(SZA,Pres, PAR0_2);
for i=1:length(NL)
    N=NL(i);
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ_2, LAI, Tl, Tldiff, ambCO2, Vcmax0_25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);
    An_ML_2(i,1)=ML.An_tot;
    clear DF ML
end

% High light level
PAR0_3=2000; % top canop irradiance, in umol/m2/s
LQ_3=Func_Light_Partitioning(SZA,Pres, PAR0_3);
for i=1:length(NL)
    N=NL(i);
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ_3, LAI, Tl, Tldiff, ambCO2, Vcmax0_25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);
    An_ML_3(i,1)=ML.An_tot;
    clear DF ML
end

figure('color','white');
plot(NL,An_ML_1,'k-o','LineWidth',1, 'MarkerSize',6);
hold on
plot([0 NL(1, end)+1],[An_ML_1(end,1) An_ML_1(end,1)],'--','LineWidth',1, 'color', [0.5 0.5 0.5]);
xlabel('Layer N','fontsize',14);
ylabel('GPP (umol/m2/s)','fontsize',14);
set(gca,'fontsize',12);
axis([0 31 5 7]);

figure('color','white');
plot(NL,An_ML_2,'k-o','LineWidth',1, 'MarkerSize',6);
hold on
plot([0 NL(1, end)+1],[An_ML_2(end,1) An_ML_2(end,1)],'--','LineWidth',1, 'color', [0.5 0.5 0.5]);
xlabel('Layer N','fontsize',14);
ylabel('GPP (umol/m2/s)','fontsize',14);
set(gca,'fontsize',12);
axis([0 31 18 20]);

figure('color','white');
plot(NL,An_ML_3,'k-o','LineWidth',1, 'MarkerSize',6);
hold on
plot([0 NL(1, end)+1],[An_ML_3(end,1) An_ML_3(end,1)],'--','LineWidth',1, 'color', [0.5 0.5 0.5]);
xlabel('Layer N','fontsize',14);
ylabel('GPP (umol/m2/s)','fontsize',14);
set(gca,'fontsize',12);
axis([0 31 13 16]);

%% Step 3---Modeled GPP-PAR by using both ML and DF1997
clear DF ML
PARobs=[100:10:2150]'*cos(SZA./180.*pi);
N=30;
for i=1:length(PARobs)
    PAR0=PARobs(i,1);
    LQ=Func_Light_Partitioning(SZA,Pres, PAR0);
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vcmax0_25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);
    
    An_DF(i,1)=DF.An_tot;
    An_ML(i,1)=ML.An_tot;
        
    clear DF ML LQ PAR0    
end

figure('color','white');
plot(PARobs,An_DF(:,1),'r-','LineWidth',2);
hold on
plot(PARobs,An_ML(:,1),'b-','LineWidth',2);
xlabel('PAR(umol/m2/s)','fontsize',14);
ylabel('Photosynthesis (umol/m2/s)','fontsize',14);
set(gca,'fontsize',12);


%% Step 4--Modeled GPP-Temperature response
clear An_DF An_ML
T=[10:1:40];
N=30;

% Low PAR
PAR0=300;
for i=1:length(T)
    Tl=T(i);
    LQ=Func_Light_Partitioning(SZA,Pres, PAR0);
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vcmax0_25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);   
    An_DF_lo(i,1)=DF.An_tot;
    An_ML_lo(i,1)=ML.An_tot;
        
    clear DF ML LQ   
end

figure('color','white');
plot(T,An_DF_lo(:,1),'r-','LineWidth',1);
hold on
plot(T,An_ML_lo(:,1),'b-','LineWidth',1);
xlabel('Leaf temperature','fontsize',14);
ylabel('Photosynthesis (umol/m2/s)','fontsize',14);
set(gca,'fontsize',12);

% Medium PAR
PAR0=1000;
for i=1:length(T)
    Tl=T(i);
    LQ=Func_Light_Partitioning(SZA,Pres, PAR0);
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vcmax0_25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);   
    An_DF_med(i,1)=DF.An_tot;
    An_ML_med(i,1)=ML.An_tot;
        
    clear DF ML LQ   
end

figure('color','white');
plot(T,An_DF_med(:,1),'r-','LineWidth',1);
hold on
plot(T,An_ML_med(:,1),'b-','LineWidth',1);
xlabel('Leaf temperature','fontsize',14);
ylabel('Photosynthesis (umol/m2/s)','fontsize',14);
set(gca,'fontsize',12);


% High PAR
PAR0=2000;
for i=1:length(T)
    Tl=T(i);
    LQ=Func_Light_Partitioning(SZA,Pres, PAR0);
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vcmax0_25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);   
    An_DF_hi(i,1)=DF.An_tot;
    An_ML_hi(i,1)=ML.An_tot;
        
    clear DF ML LQ   
end

figure('color','white');
plot(T,An_DF_hi(:,1),'r-','LineWidth',1);
hold on
plot(T,An_ML_hi(:,1),'b-','LineWidth',1);
xlabel('Leaf temperature','fontsize',14);
ylabel('Photosynthesis (umol/m2/s)','fontsize',14);
set(gca,'fontsize',12);


%% Step 5--Modeled GPP-LAI
clear An_DF_lo An_ML_lo An_DF_med An_ML_med An_DF_hi An_med_hi
Tl=28;
LAIx=[1:0.1:7];
% Low PAR
PAR0=300;
for i=1:length(LAIx)
    LAI=LAIx(i);
    LQ=Func_Light_Partitioning(SZA,Pres, PAR0);
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vcmax0_25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);   
    An_DF_lo(i,1)=DF.An_tot;
    An_ML_lo(i,1)=ML.An_tot;
        
    clear DF ML LQ   
end

figure('color','white');
plot(LAIx,An_DF_lo(:,1),'r-','LineWidth',1);
hold on
plot(LAIx,An_ML_lo(:,1),'b-','LineWidth',1);
xlabel('Leaf Area Index','fontsize',14);
ylabel('Photosynthesis (umol/m2/s)','fontsize',14);
set(gca,'fontsize',12);


% Mid PAR
PAR0=1000;
for i=1:length(LAIx)
    LAI=LAIx(i);
    LQ=Func_Light_Partitioning(SZA,Pres, PAR0);
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vcmax0_25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);   
    An_DF_med(i,1)=DF.An_tot;
    An_ML_med(i,1)=ML.An_tot;
        
    clear DF ML LQ   
end

figure('color','white');
plot(LAIx,An_DF_med(:,1),'r-','LineWidth',1);
hold on
plot(LAIx,An_ML_med(:,1),'b-','LineWidth',1);
xlabel('Leaf Area Index','fontsize',14);
ylabel('Photosynthesis (umol/m2/s)','fontsize',14);
set(gca,'fontsize',12);


% Hi PAR
PAR0=2000;
for i=1:length(LAIx)
    LAI=LAIx(i);
    LQ=Func_Light_Partitioning(SZA,Pres, PAR0);
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vcmax0_25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);   
    An_DF_hi(i,1)=DF.An_tot;
    An_ML_hi(i,1)=ML.An_tot;
        
    clear DF ML LQ   
end

figure('color','white');
plot(LAIx,An_DF_hi(:,1),'r-','LineWidth',1);
hold on
plot(LAIx,An_ML_hi(:,1),'b-','LineWidth',1);
xlabel('Leaf Area Index','fontsize',14);
ylabel('Photosynthesis (umol/m2/s)','fontsize',14);
set(gca,'fontsize',12);



% 
% %% Step 2--Call Light Partitioning Function
% PAR0=1320; % top canop irradiance, in umol/m2/s
% LQ=Func_Light_Partitioning(SZA,Pres, PAR0);
% 
% %% Step3--Call Multi-Layer Photosynthesis Model
% [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vcmax0_25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf);
% 
% 
% 
% %% Step5--Canopy Photosynthesis vs. PAR
% clear DF ML
% PARobs=[10:10:2150]'*cos(SZA./180.*pi);
% N=20;
% for i=1:length(PARobs)
%     PAR0=PARobs(i,1);
%     LQ=Func_Light_Partitioning(SZA,Pres, PAR0);
%     [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vcmax0_25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf);
%     
%     An_DF(i,1)=DF.An_tot;
%     An_ML(i,1)=ML.An_tot;
%         
%     clear DF ML LQ PAR0    
% end
% 
% figure('color','white');
% plot(PARobs,An_DF(:,1),'r-','LineWidth',2);
% hold on
% plot(PARobs,An_ML(:,1),'b-','LineWidth',2);
% xlabel('PAR(umol/m2/s)','fontsize',14);
% ylabel('Canopy Photosynthesis','fontsize',14);
% set(gca,'fontsize',12);
% 
% 
% %% Step6--Canopy Photosynthesis vs. PAR on different LAI
% clear An_ML
% Vcmax0_25=60;
% PARobs=[10:10:2150]'*cos(SZA./180.*pi);
% N=20;
% CI=0.8;
% LAI0=[1 2 3 4 5 6]';
% for i=1:length(LAI0)
%     for j=1:length(PARobs)
%         PAR0=PARobs(j,1);
%         LQ=Func_Light_Partitioning(SZA,Pres, PAR0);
%         LAI=LAI0(i,1);
%         [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vcmax0_25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf);
%         An_ML(j,i)=ML.An_tot;
%         clear DF ML
%     end
% end
% 
% figure('color','white');
% plot(PARobs,An_ML(:,1),'-','LineWidth',2,'color',[0 0 0]);
% hold on
% plot(PARobs,An_ML(:,2),'-','LineWidth',2,'color',[0 0 1]);
% plot(PARobs,An_ML(:,3),'-','LineWidth',2,'color',[0 1 0]);
% plot(PARobs,An_ML(:,4),'-','LineWidth',2,'color',[0 1 1]);
% plot(PARobs,An_ML(:,5),'-','LineWidth',2,'color',[1 0 1]);
% plot(PARobs,An_ML(:,6),'-','LineWidth',2,'color',[1 0 0]);
% xlabel('PAR(umol/m2/s)','fontsize',14);
% ylabel('Canopy Photosynthesis','fontsize',14);
% set(gca,'fontsize',12);
