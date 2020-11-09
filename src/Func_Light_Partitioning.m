function LQ=Func_Light_Partitioning(SZA, P, PAR)
% Reference: Weiss and Norman, 1985
% SZA: Solar Zenith Angle
% P: air pressure
% PAR: measured total light 

RT=PAR*0.219./(0.46);

theta=SZA./180*3.1415926;
P0=101325; % unit: pa
m=1./cos(theta);

% follow Beer's law, the expected visble direct beam radiation under clear sky, RDV
RDV=600*exp(-0.185.*(P/P0).*m).*cos(theta); % unit: W/m2

% expected visible diffuse radiation under clear sky; RdV
RdV=0.4*(600-RDV./cos(theta)).*cos(theta);

% define parameter antilog10
antilog10=10.^(-1.1950+0.4459*log(m)./log(10)-0.0345.*(log(m)./log(10)).*(log(m)./log(10)));

% w--is water absorption in the near infreared for 10 mm of precipitable
% water, under clear sky (adapted from Wang, 1976)
w=1320*antilog10;

% expected direct-beam of near-infrared radiation under clear sky 
RDN=(720*exp(-0.06*(P./P0).*m)-w).*cos(theta);

% expected diffuse-beam of near-infrared radiation under clear sky
RdN=0.6*(720-RDN./cos(theta)-w).*cos(theta);


RDV=max(RDV,0); RdV=max(RdV,0);
RDN=max(RDN,0); RdN=max(RdN,0);
% total visible light expected under clear sky
RV=RDV+RdV; 

% total nir light expected under clear sky
RN=RDN+RdN;

% total Visible light 
SV=RT.*(RV./(RV+RN)); SV=max(0, SV);

% total NIR light
SN=RT.*(RN./(RV+RN)); SN=max(0, SN);

Ratio=RT./(RV+RN); % the ratio between total measured light and total modeled light


A=0.9;
B=0.7;
C=0.88;
D=0.68;

s1=(A-Ratio)./B; s1(s1<0)=0;
s1a=1-s1.^(2/3);

s2=(C-Ratio)./D; s2(s2<0)=0;
s2a=1-s2.^(2/3);

fV=RDV./RV.*s1a; % fraction of visble direct beam
fN=RDN./RN.*s2a; % fraction of NIR direct beam

LQ.SZA=SZA;
LQ.PAR=PAR;
LQ.SV=SV*4.57;
LQ.SN=SN*4.57;
LQ.Ratio=Ratio;
LQ.fV=fV;
LQ.fN=fN;

LQ.Model_DV=SV.*fV*4.57; % direct visible light
if LQ.Model_DV<0
    LQ.Model_DV=0;
end

LQ.Model_dV=SV*4.57-SV.*fV*4.57; % diffuse visible light
if LQ.Model_dV<0
    LQ.Model_dV=0;
end

LQ.Model_DN=SN.*fN*4.57; % direct NIR light
if LQ.Model_DN<0
    LQ.Model_DN=0;
end

LQ.Model_dN=SN*4.57-SN.*fN*4.57; % diffuse NIR light
if LQ.Model_dN<0
    LQ.Model_dN=0;
end

