clc
clear

%**************************************************************************
% Fig. 3 fo Wu et al. GCB
% Jin Wu, June 2016, Brookhaven National Laboratory
%**************************************************************************

%% Step 1--default model parameters
FLAG=1; % model version control; 1--Lloyd et al. 2010 Vc-LAI relationship; 2--Mercado et al. 2006 Vc-LAI relationship;
SZA=30; % Solar Zenith Angle, in degree
Pres=10.^5; % Atmosphere Pressure, in pa
LAI=6; % Leaf Area Index

Tl=25; % leaf temperature for sunlit leaf
Tldiff=0; % leaf temperature difference between sunlit and shade leaves

ambCO2=380; % Ambient CO2 in ppm
Vcmax0_25=40; % Bonan et al., 2012 for the tropcis
CI=0.63; % Clumping index, from Chen etal., 2005 for tropical evergreen forests
Topt=35; % optimal leaf temperature for the tropics, from Lloyd and Farquhar, 2008

N=20; % number of layers for Multi-Layer Canopy Photosynthesis Modeling

% variables related to sun/shade leaf maximum intrinsic quantum yield
Phi_sun=0.75; 
PSII_sun=0.75; 
Phi_shade=0.75;
PSII_shade=0.75; 

% varaibles related to leaf age effect; scale factor 
sf_sun=1; 
sf_shade=1;
sf=1;

%% Step 2--Call Light Partitioning Function
PAR0=1320; % top canop irradiance, in umol/m2/s
LQ=Func_Light_Partitioning(SZA,Pres, PAR0);

%% Step3--Call Multi-Layer Photosynthesis Model
% 3.1 LAI + Leaf Demography
LAIx=[6.009206308	0.555485554	1.613756699	3.839964055
5.983849077	0.454831417	1.170120819	4.358896841
5.901438077	0.296062423	0.928998582	4.676377072
5.717598154	0.261296962	0.737115359	4.719185833
5.734714285	0.272848424	0.569264035	4.892601826
5.387954154	0.420247639	0.442926767	4.524779748
5.626312123	1.063671279	0.387281042	4.175359802
5.882420154	1.623401902	0.567111983	3.691906269
5.958491846	1.579914004	1.189655359	3.188922483
6.049777877	1.240243671	1.917665461	2.891868745
6.072599385	0.887836877	2.349321724	2.835440784
6.039001054	0.65579697	2.094237853	3.288966231];

% AY=14.1;
% AM=33.6;
% AO=21.7;

% AY=14*(1+0.5)/2;
AY=10;
AM=34;
AO=22;

% Phy=[9.5 32.7 20.8];

A_demo=LAIx(:,2)*AY+LAIx(:,3)*AM+LAIx(:,4)*AO;
A_demo=A_demo./LAIx(:,1);
A_demo=A_demo./mean(A_demo);

for i=1:12
    Vcmax0_25=40*A_demo(i,1);
    LAI=LAIx(i,1);
    [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vcmax0_25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf);
    
    DF_An_tot(i,1)=DF.An_tot;
    DF_Rd_tot(i,1)=DF.Rd_tot;
    DF_An_sun(i,1)=DF.An_sun;
    DF_An_shade(i,1)=DF.An_shade;
    DF_Vcmax_sun(i,1)=DF.Vcmax_sun;
    DF_Vcmax_shade(i,1)=DF.Vcmax_shade;
    DF_LAI_sun(i,1)=DF.LAI_sun;
    DF_LAI_shade(i,1)=DF.LAI_shade;   
    
    ML_An_tot(i,1)=ML.An_tot;
    ML_Rd_tot(i,1)=ML.Rd_tot;
    ML_An_sun(i,1)=ML.An_sun;
    ML_An_shade(i,1)=ML.An_shade;
    ML_V25_sun(i,1)=ML.V25_sun;
    ML_V25_shade(i,1)=ML.V25_shade;
    ML_LAI_sun(i,1)=ML.LAI_sun;
    ML_LAI_shade(i,1)=ML.LAI_shade;  
    ML_I_sun(i,1)=ML.I_sun;
    ML_I_shade(i,1)=ML.I_shade;  
    
    clear DF ML
end

DF_quan_qual=[DF_An_tot DF_Rd_tot DF_An_sun DF_An_shade DF_Vcmax_sun DF_Vcmax_shade DF_LAI_sun DF_LAI_shade];
clear DF_An_tot DF_Rd_tot DF_An_sun DF_An_shade DF_Vcmax_sun DF_Vcmax_shade DF_LAI_sun DF_LAI_shade

ML_quan_qual=[ML_An_tot ML_Rd_tot ML_An_sun ML_An_shade ML_V25_sun ML_V25_shade ML_LAI_sun ML_LAI_shade ML_I_sun ML_I_shade];
clear ML_An_tot ML_Rd_tot ML_An_sun ML_An_shade ML_V25_sun ML_V25_shade ML_LAI_sun ML_LAI_shade ML_I_sun ML_I_shade

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

FPAR=0.95-exp(-0.5./cos(30/180.0*pi)*LAIx(:,1));
A_demo1=A_demo.*FPAR;

DG2008=ML_quan_qual(:,7)*10.*A_demo+ML_quan_qual(:,8)*3.*A_demo;

MO=1:12;
figure('color','white');
plot(MO,DF_quan_qual(:,1)./max(DF_quan_qual(:,1)),'r-o','MarkerSize',5,'LineWidth',1);
hold on
plot(MO,ML_quan_qual(:,1)./max(ML_quan_qual(:,1)),'b-+','MarkerSize',6,'LineWidth',1);
plot(MO,A_demo1./max(A_demo1),'g->','MarkerSize',5,'LineWidth',1);
plot(MO,PC./max(PC),'k-s','MarkerSize',5,'LineWidth',1);
plot(MO,DG2008./max(DG2008),'c-d','MarkerSize',5,'LineWidth',1);


% %% Step4--Canopy Photosynthesis vs. ML
% NL=[1:1:30]'; % number of Layer
% for i=1:length(NL)
%     N=NL(i);
%     [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vcmax0_25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf);
%     An_ML(i,1)=ML.An_tot;
%     clear DF ML
% end
% 
% figure('color','white');
% plot(NL,An_ML,'k-o','LineWidth',2);
% xlabel('Number of Layer','fontsize',14);
% ylabel('Canopy Photosynthesis','fontsize',14);
% set(gca,'fontsize',12);


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
