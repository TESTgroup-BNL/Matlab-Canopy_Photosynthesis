function LRT=Func_Canopy_Radiance_Transfer(FLAG, SZA, LAI, Ib0, Id0, Vcmax0_25, CI)
% Goals: using revised DF1997 model to partition canopy LAI into sunlit/shade leaves LAI and partition canopy Vcmax into sunlit/shade leaves Vcmax
%        clumping index was added to original DF1997 model, following the reference from Ryu et al. 2011

%% Input
% FLAG: Model version controller;  1--Lloyd et al. 2010 Model for Vcmax-LAI relationship; 2--Mercado et al. 2006 Model for Vcmax-LAI relationship in the tropics
% SZA: solar zenith angle, in degree
% LAI: Canopy leaf area index
% Ib0: direct beam at canopy top
% Id0: diffuse irradiance at canopy top
% Vcmax0_25: Vcmax at reference 25 centi-degree for canopy top leaves
% CI: Clumping inedx; 0.63 for tropical evergreen forests (Chen et al, 2005)

%% Output
% output=[Ib0+Id0 Ib0 Id0 Lsun Lshade Ic Isun Ishade Vc Vcsun Vcshade];
% Ib0+Id0: total PAR above-canopy
% Ib0: direct beam above-canopy
% Id0: diffuse beam above-canopy
% Lsun: sunlit leaves LAI
% Lshade: shade leaves LAI
% Ic: canopy-scale total absorbed PAR
% Isun: canopy-scale sunlit leaves absorbed PAR
% Ishade: canopy-scale shade leaves absorbed PAR
% Vc: canopy scale integreated vcmax at 25 centi-degree
% Vcsun: canopy scale integrated sunlit leaves vcmax at 25 centi-degree
% Vcshade: canopy scale integreated shade leaves vcmax at 25 centi-degree

%  or

% LRT.Lsun: Sunlit LAI
% LRT.Lshade: Shade LAI
% LRT.Ic:  Canopy total absorbed irradiance
% LRT.Isun: Sunlit leaf absorbed irradiance
% LRT.Ishade: Shade leaf absorbed irradiance
% LRT.Vc: Canop total Vcmax
% LRT.Vcsun: Sunlit leaf Vcmax
% LRT.Vcshade: Shade leaf Vcmax


%% Reference: dePury and Farquhar, 1997; Ryu et al., 2011
pi=3.1415926;

%% 1 Sun/shade LAI paritioning
G=0.5; % the parameter used in Depury and Farquhar's model
kb=G./cos(SZA./180*pi); % extinction coefficient; G refers to G function, could be 0.5 to simplify it
Lsun=(1-exp(-kb.*LAI.*CI))./kb; % Sun LAI; consistent with Eq. 18a in Chen et al. 1999
Lshade=LAI-Lsun; % Shade LAI

%% 2 Total irradiance absorbed by the canopy (Eq. 13; dePury and Farquhar, 1997)
rocb=0.029; % canopy reflection coefficient for beam PAR; assuming it is the same as for diffuse PAR; 0.036 is apparently too low
rocd=0.036; % canopy reflection coefficient for diffuse PAR; assuming it is the same as for beam PAR; 0.036 is apparently too low
kb1=0.46./cos(SZA./180*pi); % beam and scattered beam PAR extinction coefficient
kd1=0.719; % diffuse and scattered diffuse PAR exintinction coefficient
Ic=(1-rocb).*Ib0.*(1-exp(-kb1.*LAI.*CI))+(1-rocd)*Id0*(1-exp(-kd1*LAI.*CI)); % canopy-scale irradiance absorbed; Eq. 2 in Ryu et al. 2011

%% 3 Calculating sunlit leave irradiance
sigma=0.15; % leaf scattering coefficient of pAR, ro1+tao1; where ro1=0.1 for leaf reflection; tao1=0.05 for leaf transmissivity
Isun1=Ib0.*(1-sigma).*(1-exp(-kb.*LAI.*CI));% direct beam, sunlit leaves; Eq. 3 from Ryu et al. 2011
Isun2=Id0.*(1-rocd).*(1-exp(-(kd1+kb).*LAI.*CI)).*kd1./(kd1+kb); % diffuse beam, sunlit leaves; Eq. 4 from Ryu et al. 2011
Isun3=Ib0.*((1-rocb).*(1-exp(-(kb1+kb).*LAI.*CI)).*kb1./(kb1+kb)-(1-sigma).*(1-exp(-2*kb.*LAI.*CI))/2); %scattering beam absored by sunlit leaves; Eq. 5 from Ryu et al. 2011
Isun=Isun1+Isun2+Isun3; 

% compared with Ryu et al. 2011, here we didn't consider the downward sunlight which is reflected by the soil/litter in the forest floor
% This is probably reasonable, as in the tropics, the radiance reached to
% forest floor is very small and negilible. 
Ishade=Ic-Isun;


%% 4 Photosynthetic capacity of sunlit and shaded leaf fraction
% Vertical distribution of leaf-level Vcmax
if FLAG==1
   % updated relationship based on the Lloyd et al. 2010, cieted by Bonan et al. 2014
   kn=exp(0.00963*Vcmax0_25-2.43);
   Vcsun=CI.*Vcmax0_25./(kn+kb.*CI).*(1-exp(-(kb.*CI+kn)*LAI)); % sunlit leaves Vcmax
   Vc=Vcmax0_25./kn.*(1-exp(-kn.*LAI)); % canopy scale vcmax; consistent with eqn. 32 in Ryu et al. 2011
   Vcshade=Vc-Vcsun;
elseif FLAG==2
   % updated relationship based on the Mercado et al. 2006
   kn=0.1823;
   Vcsun=CI.*Vcmax0_25./(kn+kb.*CI).*(1-exp(-(kb.*CI+kn)*LAI)); % sunlit leaves Vcmax
   Vc=Vcmax0_25./kn.*(1-exp(-kn.*LAI)); % canopy scale vcmax
   Vcshade=Vc-Vcsun; 
elseif FLAG==3
   kn=0.1823*6;
   Vcsun=LAI*Vcmax0_25*1/(kn+kb*LAI)*(1-exp(-CI*(kn+kb*LAI)));
   Vc=LAI*Vcmax0_25./kn.*(1-exp(-kn.*LAI)); % canopy scale vcmax
   Vcshade=Vc-Vcsun;       
end

% output=[Ib0+Id0 Ib0 Id0 Lsun Lshade Ic Isun Ishade Vc Vcsun Vcshade];
LRT.PAR0=Ib0+Id0;
LRT.Ib0=Ib0;
LRT.Id0=Id0;
LRT.Lsun=Lsun; % Sunlit LAI
LRT.Lshade=Lshade; % Shade LAI
LRT.Ic=Ic; % Canopy total absorbed irradiance
LRT.Isun=Isun; % Sunlit leaf absorbed irradiance
LRT.Ishade=Ishade; % Shade leaf absorbed irradiance
LRT.Vc=Vc; % Canop total Vcmax
LRT.Vcsun=Vcsun; % Sunlit leaf Vcmax
LRT.Vcshade=Vcshade; % Shade leaf Vcmax
