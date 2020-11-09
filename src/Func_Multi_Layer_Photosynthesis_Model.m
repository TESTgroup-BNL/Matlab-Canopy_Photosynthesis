function [DF, ML]=Func_Multi_Layer_Photosynthesis_Model(FLAG, SZA, Pres, LQ, LAI, Tl, Tldiff, ambCO2, Vcmax0_25, CI, Topt, N, Phi_sun, PSII_sun, Phi_shade, PSII_shade, sf_sun, sf_shade, sf, LAI_cut)
%% Goals: perform Multiple-Layer Photosynthesis Model under given each given set of model input

%% Input
% FLAG: model version controller;  % 1--Lloyd et al. 2010 Vc-LAI relationship; 2--Mercado et al. 2006 Vc-LAI relationship;
% SZA: Solar zenith angle, in degree; 0 for nadir view
% Pres: atmosphere pressure, in pa
% LQ: Light partitioning following the function of "Light_Quality_Partitioning_up.m"
% LAI: leaf area index of a given month
% Tl: leaf temperature, in centi-degree
% Tldiff: leaf temperature difference between sunlit leaves and shade leaves; in centi-degree
% ambCO2: ambient CO2 concentration, in ppm 
% Vcmax0_25: the Vcmax for the canopy top leaves at reference temperature 25 centi-degree; value from Bonan et al., 2012
% CI: clumping index for tropical evergreen forests, from Chen et al. 2005  
% Topt: optimal leaf temperature for tropical evergreen forests, from Lloyd and Farquhar, 2008
% N--number of layers for Multi-Layer Photosynthesis Modeling
% Phi_sun: the curvature factor for light response curves for sunlit leaves 
% PSII_sun: maximum quantum yield for sunlit leaves
% Phi_shade: the curvature factor for light response curves for sunlit leaves 
% PSII_shade: maximum quantum yield for shade leaves
% sf_sun: scaling factor for sunlit leaves, due to leaf age effect
% sf_shade: scaling factor for shade leaves, due to leaf age effect
% sf: scaling factor due to leaf age effect, assuming no phenological partitioning across vertical canopy profile
% LAI_cut: cut off of the top and bottom canopy

%% Output
% DF.An_tot: total An for DF1997 model
% DF.Rd_tot: total Rd for DF1997 model
% DF.An_sun: Sunlit leaf An for DF1997 model
% DF.An_shade: Shade leaf An for DF1997 model
% DF.Vcmax_sun: Sunlit leaf Vcmax for DF1997 model
% DF.Vcmax_shade: Shade leaf Vcmax for DF1997 model
% DF.LAI_sun: Sunlit leaf LAI for DF1997 model
% DF.LAI_shade: Shade leaf LAI for DF1997 model
% 
% ML.Profile: All important variables along vertical profiles in MLCan
% ML.V25_sun: Sunlit leaf V25 for MLCan
% ML.V25_shade: Shade leaf V25 for MLCan
% ML.An_sun: Sunlit leaf An for MLCan
% ML.An_shade: Shade leaf An for MLCan
% ML.An_tot: Total An for MLCan
% ML.I_sun: Sunlit leaf total absorbed light for MLCan
% ML.I_shade: Shade leaf total absorbed light for MLCan
% ML.Rd_tot: Total dark respiration rate for MLCan
% ML.LAI_sun: Sunlit leaf LAI for MLCan
% ML.LAI_shade: Shade Leaf LAI for MLCan

%% Step 1--Calculate vertical distribution of LAI_sun/shade, PAR_sun/shade, Vcmax_sun/shade
Ib0=LQ.Model_DV;
Id0=LQ.Model_dV;

for i=1:N
    LAIi=i./N*LAI;
    LAIi_1=(i-1)./(N*LAI);
    LRT=Func_Canopy_Radiance_Transfer(FLAG, SZA, LAIi, Ib0, Id0, Vcmax0_25, CI);
    %test(i,:)=canopy_radiance_transfer_up(FLAG, SZA, LAIi, Ib0, Id0, Vcmax0_25, CI);
    test(i,:)=[LRT.PAR0 LRT.Ib0 LRT.Id0 LRT.Lsun LRT.Lshade LRT.Ic LRT.Isun LRT.Ishade LRT.Vc LRT.Vcsun LRT.Vcshade]; % test is a temporary matrix for data storage 
    %         1     2   3   4     5     6  7     8    9   10     11
    %test=[Ib0+Id0 Ib0 Id0 Lsun Lshade Ic Isun Ishade Vc Vcsun Vcshade];
    if FLAG==1
       kn=exp(0.00963*Vcmax0_25-2.43); % using the vertical distribution of Vcmax-LAI relationship from Lloyd et al., 2010
       Vc(i,1)=(Vcmax0_25*exp(-kn*LAIi)+Vcmax0_25*exp(-kn*LAIi_1))./2;  % calculate the average between the two layers
       
       %% using 2.5 as a scaling factor 
       if LAIi<=LAI_cut
          Vc(i,1)=Vc(i,1)*sf_sun;   
       else
          Vc(i,1)=Vc(i,1)*sf_shade; 
       end
    elseif FLAG==2
       kn=0.1823;
       Vc(i,1)=(Vcmax0_25*exp(-kn*LAIi)+Vcmax0_25*exp(-kn*LAIi_1))./2; 
       if LAIi<=LAI_cut
          Vc(i,1)=Vc(i,1)*sf_sun;   
       else
          Vc(i,1)=Vc(i,1)*sf_shade; 
       end
    elseif FLAG==3
        kn=0.1823*6;
        Vc(i,1)=Vcmax0_25*(exp(-kn*LAIi*CI/6)+exp(-kn*LAIi_1*CI/6))/2;
        if LAIi<=LAI_cut
          Vc(i,1)=Vc(i,1)*sf_sun;   
       else
          Vc(i,1)=Vc(i,1)*sf_shade; 
       end
    end
    
    clear LAIi LRT LAIi_1
end


%% Calculate the vertical distribution of light, and LAI
    %         1     2   3   4     5     6  7     8    9   10     11
    %test=[Ib0+Id0 Ib0 Id0 Lsun Lshade Ic Isun Ishade Vc Vcsun Vcshade];
    
for i=1:N
    if i==1
       Profile(i,1:4)=test(i,[4 5 7 8]); % 1-LAIsun; 2-LAIshade; 4-PARsun; 5-PARshade
       Profile(i,5)=Profile(i,3)./Profile(i,1); % sunlit PAR
       Profile(i,6)=Profile(i,4)./Profile(i,2); % shade PAR
       Profile(i,7)=Vc(i,1);
       Profile(i,8)=Profile(i,1)*Vc(i,1); % sunlit Vcmax
       Profile(i,9)=Profile(i,2)*Vc(i,1); % shade Vcmax     
    else
       Profile(i,1:4)=test(i,[4 5 7 8])-test(i-1,[4 5 7 8]); % 1-LAIsun; 2-LAIshade; 4-PARsun; 5-PARshade
       Profile(i,5)=Profile(i,3)./Profile(i,1); % sunlit PAR
       Profile(i,6)=Profile(i,4)./Profile(i,2); % shade PAR
       Profile(i,7)=Vc(i,1);
       Profile(i,8)=Profile(i,1)*Vc(i,1); % sunlit Vcmax25
       Profile(i,9)=Profile(i,2)*Vc(i,1); % shade Vcmax25
    end
end



%% Caculate the vertical distribution of photosynthesis rate
for i=1:N
 %% Version 1: assume sun/shade leaves have different physiology (quantum yield)  
%     %% Sunlit leaf at Layer i
%     V25=Profile(i,7); J25=1.67*V25; T=Tl; I=Profile(i,5); Ci=ambCO2*0.7; 
%     FvCB_sun=Func_Leaf_FvCB_Photosynthesis_Model(V25, J25, T, Topt, I, Ci, Pres, PSII_sun, Phi_sun);
%     Profile(i,10)=FvCB_sun.An; % sunlit An
%     Profile(i,11)=FvCB_sun.Rd; % sunlit Rd
%     Profile(i,12)=FvCB_sun.Wc; % sunlit Wc
%     Profile(i,13)=FvCB_sun.Wj; % sunlit Wj
%     Profile(i,14)=FvCB_sun.Vcmax; % sunlit Vmax   
%     
%     %% Shade leaf at Layer i
%     V25=Profile(i,7); J25=1.67*V25; T=Tl-Tldiff; I=Profile(i,6); Ci=ambCO2*0.7;
%     FvCB_shade=Func_Leaf_FvCB_Photosynthesis_Model(V25, J25, T, Topt, I, Ci, Pres, PSII_shade, Phi_shade);
%     Profile(i,15)=FvCB_shade.An; % shade An
%     Profile(i,16)=FvCB_shade.Rd; % shade Rd
%     Profile(i,17)=FvCB_shade.Wc; % shade Wc
%     Profile(i,18)=FvCB_shade.Wj; % shade Wj
%     Profile(i,19)=FvCB_shade.Vcmax; % shade Vmax
%     
%     Profile(i,20)=Profile(i,10)*Profile(i,1); % Sunlit leaf photosynthesis
%     Profile(i,21)=Profile(i,15)*Profile(i,2); % Shade leaf photosynthesis
%     Profile(i,22)=Profile(i,20)+Profile(i,21); % leaf level total photosynthesis
%     Profile(i,23)=Profile(i,14)*Profile(i,1); % Sunlit leaf Vmax
%     Profile(i,24)=Profile(i,19)*Profile(i,2); % Shade leaf Vmax  
%     Profile(i,25)=Profile(i,11)*Profile(i,1)+Profile(i,16)*Profile(i,2); % leaf level total respiration    
 
 %% Version 2: assume top canopy and understory have different physiology (quantum yield)
    if sum(Profile(1:i,1)+Profile(1:i,2))<=LAI_cut
        %% Sun leaf at Layer i
        V25=Profile(i,7); J25=1.67*V25; T=Tl; I=Profile(i,5); Ci=ambCO2*0.7; 
        FvCB_sun=Func_Leaf_FvCB_Photosynthesis_Model(V25, J25, T, Topt, I, Ci, Pres, PSII_sun, Phi_sun);
        Profile(i,10)=FvCB_sun.An; % sunlit An
        Profile(i,11)=FvCB_sun.Rd; % sunlit Rd
        Profile(i,12)=FvCB_sun.Wc; % sunlit Wc
        Profile(i,13)=FvCB_sun.Wj; % sunlit Wj
        Profile(i,14)=FvCB_sun.Vcmax; % sunlit Vmax   
    
        %% Shade leaf at Layer i
        V25=Profile(i,7); J25=1.67*V25; T=Tl-Tldiff; I=Profile(i,6); Ci=ambCO2*0.7;
        FvCB_shade=Func_Leaf_FvCB_Photosynthesis_Model(V25, J25, T, Topt, I, Ci, Pres, PSII_sun, Phi_sun);
        Profile(i,15)=FvCB_shade.An; % shade An
        Profile(i,16)=FvCB_shade.Rd; % shade Rd
        Profile(i,17)=FvCB_shade.Wc; % shade Wc
        Profile(i,18)=FvCB_shade.Wj; % shade Wj
        Profile(i,19)=FvCB_shade.Vcmax; % shade Vmax
    
        Profile(i,20)=Profile(i,10)*Profile(i,1); % Sunlit leaf photosynthesis
        Profile(i,21)=Profile(i,15)*Profile(i,2); % Shade leaf photosynthesis
        Profile(i,22)=Profile(i,20)+Profile(i,21); % leaf level total photosynthesis
        Profile(i,23)=Profile(i,14)*Profile(i,1); % Sunlit leaf Vmax
        Profile(i,24)=Profile(i,19)*Profile(i,2); % Shade leaf Vmax  
        Profile(i,25)=Profile(i,11)*Profile(i,1)+Profile(i,16)*Profile(i,2); % leaf level total respiration    
    else
        V25=Profile(i,7); J25=1.67*V25; T=Tl; I=Profile(i,5); Ci=ambCO2*0.7; 
        FvCB_sun=Func_Leaf_FvCB_Photosynthesis_Model(V25, J25, T, Topt, I, Ci, Pres, PSII_shade, Phi_shade);
        Profile(i,10)=FvCB_sun.An; % sunlit An
        Profile(i,11)=FvCB_sun.Rd; % sunlit Rd
        Profile(i,12)=FvCB_sun.Wc; % sunlit Wc
        Profile(i,13)=FvCB_sun.Wj; % sunlit Wj
        Profile(i,14)=FvCB_sun.Vcmax; % sunlit Vmax   
    
        %% Shade leaf at Layer i
        V25=Profile(i,7); J25=1.67*V25; T=Tl-Tldiff; I=Profile(i,6); Ci=ambCO2*0.7;
        FvCB_shade=Func_Leaf_FvCB_Photosynthesis_Model(V25, J25, T, Topt, I, Ci, Pres, PSII_shade, Phi_shade);
        Profile(i,15)=FvCB_shade.An; % shade An
        Profile(i,16)=FvCB_shade.Rd; % shade Rd
        Profile(i,17)=FvCB_shade.Wc; % shade Wc
        Profile(i,18)=FvCB_shade.Wj; % shade Wj
        Profile(i,19)=FvCB_shade.Vcmax; % shade Vmax
    
        Profile(i,20)=Profile(i,10)*Profile(i,1); % Sunlit leaf photosynthesis
        Profile(i,21)=Profile(i,15)*Profile(i,2); % Shade leaf photosynthesis
        Profile(i,22)=Profile(i,20)+Profile(i,21); % leaf level total photosynthesis
        Profile(i,23)=Profile(i,14)*Profile(i,1); % Sunlit leaf Vmax
        Profile(i,24)=Profile(i,19)*Profile(i,2); % Shade leaf Vmax  
        Profile(i,25)=Profile(i,11)*Profile(i,1)+Profile(i,16)*Profile(i,2); % leaf level total respiration           
    end
    
end

%% Canopy integrated photosynthesis rate
V25_sun=sum(Profile(:,8));
V25_shade=sum(Profile(:,9));
An_sun=sum(Profile(:,20));
An_shade=sum(Profile(:,21));
An_tot=sum(Profile(:,22));
I_sun=sum(Profile(:,3));
I_shade=sum(Profile(:,4));
Rd_tot=sum(Profile(:,25));
LAI_sun=sum(Profile(:,1));
LAI_shade=sum(Profile(:,2));

test_ML=[V25_sun V25_shade An_sun An_shade An_tot I_sun I_shade Rd_tot LAI_sun LAI_shade];


%% Depury and Farquhar Model 
V25=V25_sun; J25=1.67*V25; T=Tl; I=I_sun; Ci=ambCO2*0.7; 
FvCB_Csun=Func_Leaf_FvCB_Photosynthesis_Model(V25, J25, T, Topt, I, Ci, Pres, PSII_sun, Phi_sun);
     
V25=V25_shade; J25=1.67*V25; T=Tl-Tldiff; I=I_shade; Ci=ambCO2*0.7; 
FvCB_Cshade=Func_Leaf_FvCB_Photosynthesis_Model(V25, J25, T, Topt, I, Ci, Pres, PSII_shade, Phi_shade);
     
Summary1(1,1)=FvCB_Csun.An+FvCB_Cshade.An; % Total An
Summary1(1,2)=FvCB_Csun.Rd+FvCB_Cshade.Rd; % Total Rd
Summary1(1,3)=FvCB_Csun.An; % Asun
Summary1(1,4)=FvCB_Cshade.An; % Ashade
Summary1(1,5)=FvCB_Csun.Vcmax; % Vcmax_sun
Summary1(1,6)=FvCB_Cshade.Vcmax; % Vcmax_shade

DF.An_tot=Summary1(1,1);
DF.Rd_tot=Summary1(1,2);
DF.An_sun=Summary1(1,3);
DF.An_shade=Summary1(1,4);
DF.Vcmax_sun=Summary1(1,5);
DF.Vcmax_shade=Summary1(1,6);
DF.LAI_sun=LAI_sun;
DF.LAI_shade=LAI_shade;

ML.Profile=Profile;
ML.V25_sun=V25_sun;
ML.V25_shade=V25_shade;
ML.An_sun=An_sun;
ML.An_shade=An_shade;
ML.An_tot=An_tot;
ML.I_sun=I_sun;
ML.I_shade=I_shade;
ML.Rd_tot=Rd_tot;
ML.LAI_sun=LAI_sun;
ML.LAI_shade=LAI_shade;
% ML.Summary1=Summary1;
% ML.physiology1=physiology1;
% ML.test_ML=test_ML;
% ML.Profile=Profile;