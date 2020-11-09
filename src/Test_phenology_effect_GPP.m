clc
clear

% *************************************************************************
% Goals of this code:
% Implement phenology analysis
% Jin Wu, BNL, 2016-July
% *************************************************************************

%% Step 1: default model parameters
FLAG=1; % model version control; 1--Lloyd et al. 2010 Vc-LAI relationship; 2--Mercado et al. 2006 Vc-LAI relationship; 3--Ryu et al. version
SZA=30; % Solar Zenith Angle, in degree
Pres=10.^5; % Atmosphere Pressure, in pa
LAI=6; % Leaf Area Index

Tl=28; % leaf temperature for sunlit leaf

ambCO2=380; % Ambient CO2 in ppm
Vcmax0_25=40; % Bonan et al., 2012 for the tropcis
Flag_Scale=0; % 1--apply scale factor; 0--do not apply; to scale Vcmax into a given value
Topt=35; % optimal leaf temperature for the tropics, from Lloyd and Farquhar, 2008

N=30; % number of layers for Multi-Layer Canopy Photosynthesis Modeling

% variables related to sun/shade leaf maximum intrinsic quantum yield
Phi_sun=0.7; 
PSII_sun=0.7; 
Phi_shade=0.7;
PSII_shade=0.7; 

% varaibles related to leaf age effect; scale factor 
sf_sun=1; 
sf_shade=1;
sf=1;
LAI_cut=2.5;

%% Uncertainty variables:
PC_age=[7 34 22];
Tldiff=0; % leaf temperature difference between sunlit and shade leaves
CI=0.66; % Clumping index, 0.63 from Chen etal., 2005 for tropical evergreen forests; 0.66 from He et al. 2012


%% Step 2: prescribed leaf phenology
LAI_t=[6.009206308
5.983849077
5.901438077
5.717598154
5.734714285
5.387954154
5.626312123
5.882420154
5.958491846
6.049777877
6.072599385
6.039001054];

Lage_t=[0.092439089	0.268547395	0.639013516
0.076009841	0.195546513	0.728443646
0.050167844	0.157419017	0.792413139
0.045700477	0.128920456	0.825379067
0.047578381	0.099266329	0.85315529
0.077997627	0.082206855	0.839795518
0.189053017	0.068833906	0.742113077
0.27597517	0.096407936	0.627616894
0.265153338	0.199657126	0.535189536
0.205006481	0.316981136	0.478012383
0.146203762	0.386872503	0.466923735
0.108593617	0.346785476	0.544620907];

PC=[0.020264256
0.020563209
0.020025085
0.020077366
0.018945277
0.018313793
0.016658999
0.016831411
0.019045785
0.019786155
0.021301711
0.021359915];


V25_t(:,1)=PC_age(1,1)*Lage_t(:,1)+PC_age(1,2)*Lage_t(:,2)+PC_age(1,3)*Lage_t(:,3);

if Flag_Scale==1
    V25_t(:,1)=V25_t(:,1)/34*Vcmax0_25;
end


%% Step 3: PC-calculation based Wu2016 and DG2008
% calculate Wu2016 model results
cos_sza=cos(SZA./180*pi);
FAPAR=0.95-exp(-0.5*LAI_t./cos_sza);
PC_Wu(:,1)=V25_t.*FAPAR;
PC_Wu(:,2)=V25_t;
PC_Wu(:,3)=FAPAR;

% calculate DG2008 model results
Lsun=(1-exp(-0.5*CI*LAI_t./cos_sza))./(0.5/cos_sza);
Lshade=LAI_t-Lsun;
PC_DG(:,1)=10*V25_t.*Lsun+3*V25_t.*Lshade;
PC_DG(:,2)=10*V25_t+3*V25_t;
PC_DG(:,3)=10*Lsun+3*Lshade;


%% Step 4: test phenology effect
PAR0=1320; % 1320 top canop irradiance, in umol/m2/s
LQ=Func_Light_Partitioning(SZA,Pres, PAR0);

% DF1997 and ML
% 1. allow for varying leaf quantity and quality
for i=1:12
    LAI=LAI_t(i,1);
    Vc25=V25_t(i,1);
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vc25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);   
    An_DF(i,1)=DF.An_tot;
    An_ML(i,1)=ML.An_tot;
    clear DF ML
end

% 2. allow for varying leaf quality only
for i=1:12
    LAI=mean(LAI_t(:,1));
    Vc25=V25_t(i,1);
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vc25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);   
    An_DF(i,2)=DF.An_tot;
    An_ML(i,2)=ML.An_tot;
    clear DF ML
end

% 3. allow for varying leaf quantity only
for i=1:12
    LAI=LAI_t(i,1);
    Vc25=mean(V25_t(:,1));
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vc25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);   
    An_DF(i,3)=DF.An_tot;
    An_ML(i,3)=ML.An_tot;
    clear DF ML
end

%% Both quantity and quality
figure('color','white');
plot(1:12,An_DF(:,1)./max(An_DF(:,1)),'m-+','LineWidth',1,'MarkerSize',6);
hold on
plot(1:12,An_ML(:,1)./max(An_ML(:,1)),'b-o','LineWidth',1,'MarkerSize',5);
plot(1:12,PC./max(PC),'k-s','LineWidth',1,'MarkerSize',5);
plot(1:12,PC_Wu(:,1)./max(PC_Wu(:,1)),'r-*','LineWidth',1,'MarkerSize',6);
plot(1:12,PC_DG(:,1)./max(PC_DG(:,1)),'->','LineWidth',1,'color',[0 0.5 0],'MarkerSize',5);

xlabel('Month','fontsize',14);
ylabel('GPP (umol/m2/s)','fontsize',14);
set(gca,'fontsize',12);
axis([0 13 0.7 1.05]);

%% Step 5: calculate R2 for each scenario
Stat(1,:)=Estimate_Rsqure(PC, An_DF); % An_DF
Stat(2,:)=Estimate_Rsqure(PC, An_ML); % An_ML
Stat(3,:)=Estimate_Rsqure(PC, PC_Wu); % PC_Wu
Stat(4,:)=Estimate_Rsqure(PC, PC_DG); % PC_DG

%% quantity only
figure('color','white');
plot(1:12,An_DF(:,3)./max(An_DF(:,3)),'m-+','LineWidth',1,'MarkerSize',6);
hold on
plot(1:12,An_ML(:,3)./max(An_ML(:,3)),'b-o','LineWidth',1,'MarkerSize',5);
plot(1:12,PC./max(PC),'k-s','LineWidth',1,'MarkerSize',5);
plot(1:12,PC_Wu(:,3)./max(PC_Wu(:,3)),'r-*','LineWidth',1,'MarkerSize',6);
plot(1:12,PC_DG(:,3)./max(PC_DG(:,3)),'->','LineWidth',1,'color',[0 0.5 0],'MarkerSize',5);

xlabel('Month','fontsize',14);
ylabel('GPP (umol/m2/s)','fontsize',14);
set(gca,'fontsize',12);
axis([0 13 0.7 1.05]);

%% quality only
figure('color','white');
plot(1:12,An_DF(:,2)./max(An_DF(:,2)),'m-+','LineWidth',1,'MarkerSize',6);
hold on
plot(1:12,An_ML(:,2)./max(An_ML(:,2)),'b-o','LineWidth',1,'MarkerSize',5);
plot(1:12,PC./max(PC),'k-s','LineWidth',1,'MarkerSize',5);
plot(1:12,PC_Wu(:,2)./max(PC_Wu(:,2)),'r-*','LineWidth',1,'MarkerSize',6);
plot(1:12,PC_DG(:,2)./max(PC_DG(:,2)),'->','LineWidth',1,'color',[0 0.5 0],'MarkerSize',5);

xlabel('Month','fontsize',14);
ylabel('GPP (umol/m2/s)','fontsize',14);
set(gca,'fontsize',12);
axis([0 13 0.7 1.05]);


%% Step 6: diagnosis for the underlying process
for i=1:12
    LAI=LAI_t(i,1);
    Vc25=V25_t(i,1);
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vc25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);   
    An_DF(i,1)=DF.An_tot;
    An_ML(i,1)=ML.An_tot;
    
    Diag_LAI(i,1)=ML.LAI_sun;
    Diag_LAI(i,2)=ML.LAI_shade;
    Diag_LAI(i,3)=ML.LAI_sun+ML.LAI_shade;
    
    Diag_I(i,1)=ML.I_sun;
    Diag_I(i,2)=ML.I_shade;
    Diag_I(i,3)=ML.I_sun+ML.I_shade;
    
    Diag_A(i,1)=ML.An_sun;
    Diag_A(i,2)=ML.An_shade;
    Diag_A(i,3)=ML.An_sun+ML.An_shade;
    
    Diag_Vc(i,1)=DF.Vcmax_sun;
    Diag_Vc(i,2)=DF.Vcmax_shade;
    Diag_Vc(i,3)=DF.Vcmax_sun+DF.Vcmax_shade;
    
    clear DF ML
end

% display leaf quantity seasonality
figure('color','white');
plot(1:12,Diag_LAI(:,3),'k-o','LineWidth',1,'MarkerSize',5);
hold on
plot(1:12,Diag_LAI(:,1),'k-->','LineWidth',1,'MarkerSize',5);
plot(1:12,Diag_LAI(:,2),'-->','LineWidth',1,'MarkerSize',5,'color',[0.5 0.5 0.5]);
xlabel('Month','fontsize',14);
ylabel('LAI','fontsize',14);
set(gca,'fontsize',12);
axis([0 13 0 7]);

% display leaf quality seasonality
figure('color','white');
plot(1:12,Diag_Vc(:,3),'k-o','LineWidth',1,'MarkerSize',5);
hold on
plot(1:12,Diag_Vc(:,1),'k-->','LineWidth',1,'MarkerSize',5);
plot(1:12,Diag_Vc(:,2),'-->','LineWidth',1,'MarkerSize',5,'color',[0.5 0.5 0.5]);
xlabel('Month','fontsize',14);
ylabel('Vcmax','fontsize',14);
set(gca,'fontsize',12);
axis([0 13 0 200]);

% display PARabsorbed seasonality
figure('color','white');
plot(1:12,Diag_I(:,3),'k-o','LineWidth',1,'MarkerSize',5);
hold on
plot(1:12,Diag_I(:,1),'k-->','LineWidth',1,'MarkerSize',5);
plot(1:12,Diag_I(:,2),'-->','LineWidth',1,'MarkerSize',5,'color',[0.5 0.5 0.5]);
xlabel('Month','fontsize',14);
ylabel('PARabsorbed','fontsize',14);
set(gca,'fontsize',12);
axis([0 13 0 1300]);

% display An seasonality
figure('color','white');
plot(1:12,Diag_A(:,3),'k-o','LineWidth',1,'MarkerSize',5);
hold on
plot(1:12,Diag_A(:,1),'k-->','LineWidth',1,'MarkerSize',5);
plot(1:12,Diag_A(:,2),'-->','LineWidth',1,'MarkerSize',5,'color',[0.5 0.5 0.5]);
xlabel('Month','fontsize',14);
ylabel('An','fontsize',14);
set(gca,'fontsize',12);
axis([0 13 0 20]);


% display An seasonality-relative
figure('color','white');
plot(1:12,Diag_A(:,3)./max(Diag_A(:,3)),'k-o','LineWidth',1,'MarkerSize',5);
hold on
plot(1:12,Diag_A(:,1)./max(Diag_A(:,1)),'k-->','LineWidth',1,'MarkerSize',5);
plot(1:12,Diag_A(:,2)./max(Diag_A(:,2)),'-->','LineWidth',1,'MarkerSize',5,'color',[0.5 0.5 0.5]);
plot(1:12,PC./max(PC),'r-','LineWidth',1);
xlabel('Month','fontsize',14);
ylabel('An','fontsize',14);
set(gca,'fontsize',12);
axis([0 13 0.7 1.05]);


% display FAPAR seasonality
figure('color','white');
plot(1:12,Diag_I(:,3)./max(Diag_I(:,3)),'k-o','LineWidth',1,'MarkerSize',5);
hold on
plot(1:12,Diag_I(:,1)./max(Diag_I(:,1)),'k-->','LineWidth',1,'MarkerSize',5);
plot(1:12,Diag_I(:,2)./max(Diag_I(:,2)),'-->','LineWidth',1,'MarkerSize',5,'color',[0.5 0.5 0.5]);
plot(1:12,FAPAR./max(FAPAR),'r-','LineWidth',1);
xlabel('Month','fontsize',14);
ylabel('An','fontsize',14);
set(gca,'fontsize',12);
axis([0 13 0.7 1.05]);


%% Step 7: Disaply seasonality to allow for differential leaf phenology turnover in top vs. bottom canopy
ftop=0.8;
LAI_cut=2.5;
JW=To_Cohort_To_Layer_Model(LAI_t, Lage_t, PC_age, ftop, LAI_cut, Flag_Scale, Vcmax0_25);
for i=1:12
    Vc25=40;
    sf_sun=JW.V25_top(i,1)/Vc25; 
    sf_shade=JW.V25_bottom(i,1)/Vc25;
    LAI=LAI_t(i,1);
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vc25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);   
    An_DF(i,1)=DF.An_tot;
    An_ML(i,1)=ML.An_tot;    
end

figure('color','white');
plot(1:12,An_DF(:,1)./max(An_DF(:,1)),'m-+','LineWidth',1,'MarkerSize',6);
hold on
plot(1:12,An_ML(:,1)./max(An_ML(:,1)),'b-o','LineWidth',1,'MarkerSize',5);
plot(1:12,PC./max(PC),'k-s','LineWidth',1,'MarkerSize',5);

xlabel('Month','fontsize',14);
ylabel('GPP (umol/m2/s)','fontsize',14);
set(gca,'fontsize',12);
axis([0 13 0.7 1.05]);

% different cutoff
ftop=0.4;
LAI_cut=2;
JW=To_Cohort_To_Layer_Model(LAI_t, Lage_t, PC_age, ftop, LAI_cut, Flag_Scale, Vcmax0_25);
for i=1:12
    Vc25=40;
    sf_sun=JW.V25_top(i,1)/Vc25; 
    sf_shade=JW.V25_bottom(i,1)/Vc25;
    LAI=LAI_t(i,1);
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vc25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);   
    An_DF20(i,1)=DF.An_tot;
    An_ML20(i,1)=ML.An_tot;    
end

LAI_cut=2.5;
JW=To_Cohort_To_Layer_Model(LAI_t, Lage_t, PC_age, ftop, LAI_cut, Flag_Scale, Vcmax0_25);
for i=1:12
    Vc25=40;
    sf_sun=JW.V25_top(i,1)/Vc25; 
    sf_shade=JW.V25_bottom(i,1)/Vc25;
    LAI=LAI_t(i,1);
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vc25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);   
    An_DF25(i,1)=DF.An_tot;
    An_ML25(i,1)=ML.An_tot;    
end

LAI_cut=3.0;
JW=To_Cohort_To_Layer_Model(LAI_t, Lage_t, PC_age, ftop, LAI_cut, Flag_Scale, Vcmax0_25);
for i=1:12
    Vc25=40;
    sf_sun=JW.V25_top(i,1)/Vc25; 
    sf_shade=JW.V25_bottom(i,1)/Vc25;
    LAI=LAI_t(i,1);
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vc25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);   
    An_DF30(i,1)=DF.An_tot;
    An_ML30(i,1)=ML.An_tot;    
end

figure('color','white');
plot(1:12,An_ML25(:,1)./max(An_ML25(:,1)),'k-s','LineWidth',1,'MarkerSize',5);
hold on
plot(1:12,An_ML20(:,1)./max(An_ML20(:,1)),'r-*','LineWidth',1,'MarkerSize',6);
plot(1:12,An_ML30(:,1)./max(An_ML30(:,1)),'b-o','LineWidth',1,'MarkerSize',5);

xlabel('Month','fontsize',14);
ylabel('GPP (umol/m2/s)','fontsize',14);
set(gca,'fontsize',12);
axis([0 13 0.7 1.05]);



%% Step 8: the effect of ftop on canopy photosythesis
ftop=0.2;
LAI_cut=2.5;
JW=To_Cohort_To_Layer_Model(LAI_t, Lage_t, PC_age, ftop, LAI_cut, Flag_Scale, Vcmax0_25);
for i=1:12
    Vc25=30;
    sf_sun=JW.V25_top(i,1)/Vc25; 
    sf_shade=JW.V25_bottom(i,1)/Vc25;
    LAI=LAI_t(i,1);
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vc25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);   
    An_DF_02(i,1)=DF.An_tot;
    An_ML_02(i,1)=ML.An_tot;    
    
    clear DF ML
end

ftop=0.4;
JW=To_Cohort_To_Layer_Model(LAI_t, Lage_t, PC_age, ftop, LAI_cut, Flag_Scale, Vcmax0_25);
for i=1:12
    Vc25=30;
    sf_sun=JW.V25_top(i,1)/Vc25; 
    sf_shade=JW.V25_bottom(i,1)/Vc25;
    LAI=LAI_t(i,1);
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vc25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);   
    An_DF_04(i,1)=DF.An_tot;
    An_ML_04(i,1)=ML.An_tot;    
    
    clear DF ML
end


ftop=0.6;
JW=To_Cohort_To_Layer_Model(LAI_t, Lage_t, PC_age, ftop, LAI_cut, Flag_Scale, Vcmax0_25);
for i=1:12
    Vc25=30;
    sf_sun=JW.V25_top(i,1)/Vc25; 
    sf_shade=JW.V25_bottom(i,1)/Vc25;
    LAI=LAI_t(i,1);
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vc25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);   
    An_DF_06(i,1)=DF.An_tot;
    An_ML_06(i,1)=ML.An_tot;    
    
    clear DF ML
end


ftop=0.8;
JW=To_Cohort_To_Layer_Model(LAI_t, Lage_t, PC_age, ftop, LAI_cut, Flag_Scale, Vcmax0_25);
for i=1:12
    Vc25=30;
    sf_sun=JW.V25_top(i,1)/Vc25; 
    sf_shade=JW.V25_bottom(i,1)/Vc25;
    LAI=LAI_t(i,1);
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vc25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);   
    An_DF_08(i,1)=DF.An_tot;
    An_ML_08(i,1)=ML.An_tot;    
    
    clear DF ML
end


ftop=1.0;
JW=To_Cohort_To_Layer_Model(LAI_t, Lage_t, PC_age, ftop, LAI_cut, Flag_Scale, Vcmax0_25);
for i=1:12
    Vc25=30;
    sf_sun=JW.V25_top(i,1)/Vc25; 
    sf_shade=JW.V25_bottom(i,1)/Vc25;
    LAI=LAI_t(i,1);
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vc25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);   
    An_DF_10(i,1)=DF.An_tot;
    An_ML_10(i,1)=ML.An_tot;    
    
    clear DF ML
end

figure('color','white');
plot(1:12,An_ML_02(:,1),'r-o','LineWidth',1,'MarkerSize',5);
hold on
plot(1:12,An_ML_04(:,1),'g-o','LineWidth',1,'MarkerSize',5);
plot(1:12,An_ML_06(:,1),'b-o','LineWidth',1,'MarkerSize',5);
plot(1:12,An_ML_08(:,1),'c-o','LineWidth',1,'MarkerSize',5);
plot(1:12,An_ML_10(:,1),'m-o','LineWidth',1,'MarkerSize',5);
%plot(1:12,PC(:,1)./max(PC(:,1)),'k-s','LineWidth',1,'MarkerSize',5);

xlabel('Month','fontsize',14);
ylabel('GPP (umol/m2/s)','fontsize',14);
set(gca,'fontsize',12);
%axis([0 13 0.7 1.05]);


figure('color','white');
plot(1:12,An_ML_02(:,1)./max(An_ML_02(:,1)),'r-o','LineWidth',1,'MarkerSize',5);
hold on
plot(1:12,An_ML_04(:,1)./max(An_ML_04(:,1)),'g-o','LineWidth',1,'MarkerSize',5);
plot(1:12,An_ML_06(:,1)./max(An_ML_06(:,1)),'b-o','LineWidth',1,'MarkerSize',5);
plot(1:12,An_ML_08(:,1)./max(An_ML_08(:,1)),'c-o','LineWidth',1,'MarkerSize',5);
plot(1:12,An_ML_10(:,1)./max(An_ML_10(:,1)),'m-o','LineWidth',1,'MarkerSize',5);
plot(1:12,PC(:,1)./max(PC(:,1)),'k-s','LineWidth',1,'MarkerSize',5);

xlabel('Month','fontsize',14);
ylabel('GPP (umol/m2/s)','fontsize',14);
set(gca,'fontsize',12);


ftop=0.1;
JW=To_Cohort_To_Layer_Model(LAI_t, Lage_t, PC_age, ftop, LAI_cut, Flag_Scale, Vcmax0_25);
for i=1:12
    Vc25=30;
    sf_sun=JW.V25_top(i,1)/Vc25; 
    sf_shade=JW.V25_bottom(i,1)/Vc25;
    LAI=LAI_t(i,1);
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vc25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);   
    An_DF_01(i,1)=DF.An_tot;
    An_ML_01(i,1)=ML.An_tot;    
    
    clear DF ML
end

ftop=0.3;
JW=To_Cohort_To_Layer_Model(LAI_t, Lage_t, PC_age, ftop, LAI_cut, Flag_Scale, Vcmax0_25);
for i=1:12
    Vc25=30;
    sf_sun=JW.V25_top(i,1)/Vc25; 
    sf_shade=JW.V25_bottom(i,1)/Vc25;
    LAI=LAI_t(i,1);
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vc25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);   
    An_DF_03(i,1)=DF.An_tot;
    An_ML_03(i,1)=ML.An_tot;    
    
    clear DF ML
end

ftop=0.5;
JW=To_Cohort_To_Layer_Model(LAI_t, Lage_t, PC_age, ftop, LAI_cut, Flag_Scale, Vcmax0_25);
for i=1:12
    Vc25=30;
    sf_sun=JW.V25_top(i,1)/Vc25; 
    sf_shade=JW.V25_bottom(i,1)/Vc25;
    LAI=LAI_t(i,1);
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vc25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);   
    An_DF_05(i,1)=DF.An_tot;
    An_ML_05(i,1)=ML.An_tot;    
    
    clear DF ML
end

ftop=0.7;
JW=To_Cohort_To_Layer_Model(LAI_t, Lage_t, PC_age, ftop, LAI_cut, Flag_Scale, Vcmax0_25);
for i=1:12
    Vc25=30;
    sf_sun=JW.V25_top(i,1)/Vc25; 
    sf_shade=JW.V25_bottom(i,1)/Vc25;
    LAI=LAI_t(i,1);
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vc25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);   
    An_DF_07(i,1)=DF.An_tot;
    An_ML_07(i,1)=ML.An_tot;    
    
    clear DF ML
end


ftop=0.9;
JW=To_Cohort_To_Layer_Model(LAI_t, Lage_t, PC_age, ftop, LAI_cut, Flag_Scale, Vcmax0_25);
for i=1:12
    Vc25=30;
    sf_sun=JW.V25_top(i,1)/Vc25; 
    sf_shade=JW.V25_bottom(i,1)/Vc25;
    LAI=LAI_t(i,1);
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vc25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);   
    An_DF_09(i,1)=DF.An_tot;
    An_ML_09(i,1)=ML.An_tot;    
    
    clear DF ML
end


ftop=0.0;
JW=To_Cohort_To_Layer_Model(LAI_t, Lage_t, PC_age, ftop, LAI_cut, Flag_Scale, Vcmax0_25);
for i=1:12
    Vc25=30;
    sf_sun=JW.V25_top(i,1)/Vc25; 
    sf_shade=JW.V25_bottom(i,1)/Vc25;
    LAI=LAI_t(i,1);
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vc25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);   
    An_DF_00(i,1)=DF.An_tot;
    An_ML_00(i,1)=ML.An_tot;    
    
    clear DF ML
end


mean_ML=[mean(An_ML_00(:,1))
    mean(An_ML_01(:,1))
    mean(An_ML_02(:,1))
    mean(An_ML_03(:,1))
    mean(An_ML_04(:,1))
    mean(An_ML_05(:,1))
    mean(An_ML_06(:,1))
    mean(An_ML_07(:,1))
    mean(An_ML_08(:,1))
    mean(An_ML_09(:,1))
    mean(An_ML_10(:,1))
    ];

figure('color','white');
plot(0.0:0.1:1,mean_ML,'k-s','MarkerSize',5);
xlabel('ftop','fontsize',14);
ylabel('mean GPP','fontsize',14);
set(gca,'fontsize',12);


StatA(1,:)=Estimate_Rsqure1(PC, An_ML_00);
StatA(2,:)=Estimate_Rsqure1(PC, An_ML_01);
StatA(3,:)=Estimate_Rsqure1(PC, An_ML_02);
StatA(4,:)=Estimate_Rsqure1(PC, An_ML_03);
StatA(5,:)=Estimate_Rsqure1(PC, An_ML_04);
StatA(6,:)=Estimate_Rsqure1(PC, An_ML_05);
StatA(7,:)=Estimate_Rsqure1(PC, An_ML_06);
StatA(8,:)=Estimate_Rsqure1(PC, An_ML_07);
StatA(9,:)=Estimate_Rsqure1(PC, An_ML_08);
StatA(10,:)=Estimate_Rsqure1(PC, An_ML_09);
StatA(11,:)=Estimate_Rsqure1(PC, An_ML_10);

figure('color','white');
plot(0.0:0.1:1,StatA(:,1),'k-s','MarkerSize',5);
xlabel('ftop','fontsize',14);
ylabel('mean GPP','fontsize',14);
set(gca,'fontsize',12);


%% for Fig. 6
% allow for both quality and quantity
ftop=0.8;
JW=To_Cohort_To_Layer_Model(LAI_t, Lage_t, PC_age, ftop, LAI_cut, Flag_Scale, Vcmax0_25);
for i=1:12
    Vc25=30;
    sf_sun=JW.V25_top(i,1)/Vc25; 
    sf_shade=JW.V25_bottom(i,1)/Vc25;
    LAI=LAI_t(i,1);
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vc25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);   
    An_DF_08_A(i,1)=DF.An_tot;
    An_ML_08_A(i,1)=ML.An_tot;    
    
    clear DF ML
end

% allow for quantity
for i=1:12
    Vc25=30;
    sf_sun=mean(JW.V25_top(:,1))/Vc25; 
    sf_shade=mean(JW.V25_bottom(:,1))/Vc25;
    LAI=LAI_t(i,1);
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vc25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);   
    An_DF_08_B(i,1)=DF.An_tot;
    An_ML_08_B(i,1)=ML.An_tot;    
    
    clear DF ML
end

% allow for quality
for i=1:12
    Vc25=30;
    sf_sun=JW.V25_top(i,1)/Vc25; 
    sf_shade=JW.V25_bottom(i,1)/Vc25;
    LAI=mean(LAI_t(:,1));
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vc25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);   
    An_DF_08_C(i,1)=DF.An_tot;
    An_ML_08_C(i,1)=ML.An_tot;    
    
    clear DF ML
end

figure('color','white');
plot(1:12,An_DF_08_B/max(An_DF_08_B),'-o','color',[120 120 120]./255,'LineWidth',1,'MarkerSize',5);
hold on
plot(1:12,PC/max(PC),'k-s','LineWidth',1,'MarkerSize',5);
xlabel('Month','fontsize',14);
ylabel('GPP (umol/m2/s)','fontsize',14);
set(gca,'fontsize',12);
axis([0 13 0.7 1.05]);

figure('color','white');
plot(1:12,An_ML/max(An_ML),'-o','color',[120 120 120]./255,'LineWidth',1,'MarkerSize',5);
hold on
plot(1:12,PC/max(PC),'k-s','LineWidth',1,'MarkerSize',5);
xlabel('Month','fontsize',14);
ylabel('GPP (umol/m2/s)','fontsize',14);
set(gca,'fontsize',12);
axis([0 13 0.7 1.05]);

figure('color','white');
plot(1:12,An_ML_08_A/max(An_ML_08_A),'-o','color',[120 120 120]./255,'LineWidth',1,'MarkerSize',5);
hold on
plot(1:12,PC/max(PC),'k-s','LineWidth',1,'MarkerSize',5);
xlabel('Month','fontsize',14);
ylabel('GPP (umol/m2/s)','fontsize',14);
set(gca,'fontsize',12);
axis([0 13 0.7 1.05]);

%% for Fig. S1
% 1. allow for varying leaf quantity and quality
for i=1:12
    LAI=LAI_t(i,1);
    Vc25=V25_t(i,1)*0.5;
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vc25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);   
    An_DF_b05(i,1)=DF.An_tot;
    An_ML_b05(i,1)=ML.An_tot;
    clear DF ML
end

for i=1:12
    LAI=LAI_t(i,1);
    Vc25=V25_t(i,1)*1;
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vc25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);   
    An_DF_b10(i,1)=DF.An_tot;
    An_ML_b10(i,1)=ML.An_tot;
    clear DF ML
end

for i=1:12
    LAI=LAI_t(i,1);
    Vc25=V25_t(i,1)*1.5;
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vc25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);   
    An_DF_b15(i,1)=DF.An_tot;
    An_ML_b15(i,1)=ML.An_tot;
    clear DF ML
end

for i=1:12
    LAI=LAI_t(i,1);
    Vc25=V25_t(i,1)*2.0;
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vc25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut);   
    An_DF_b20(i,1)=DF.An_tot;
    An_ML_b20(i,1)=ML.An_tot;
    clear DF ML
end

figure('color','white');
hold on
plot(1:12,An_ML_b10(:,1)./max(An_ML_b10(:,1)),'b-o','LineWidth',1,'MarkerSize',5);
plot(1:12,An_ML_b15(:,1)./max(An_ML_b15(:,1)),'r-*','LineWidth',1,'MarkerSize',6);
plot(1:12,An_ML_b20(:,1)./max(An_ML_b20(:,1)),'->','LineWidth',1,'color',[0 0.5 0],'MarkerSize',5);
plot(1:12,PC./max(PC),'k-s','LineWidth',1,'MarkerSize',5);

xlabel('Month','fontsize',14);
ylabel('GPP (umol/m2/s)','fontsize',14);
set(gca,'fontsize',12);
axis([0 13 0.7 1.05]);