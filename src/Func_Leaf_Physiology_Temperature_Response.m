function LPT=Func_Leaf_Physiology_Temperature_Response(V25, J25, T, Topt, Press)
% Goals: scale physiological parameters at reference 25 centi-degree to any given leaf temperature

%% Input
% V25: Vcmax standardized at 25 centi-degree
% J25: Jmax standardized at 25 centi-degree
% T: leaf temperature in centi-degree
% Topt: optimal temperature for Jmax in centi-degree

%% Output
% Jmax: Jmax after temperature adjustment from J25
% Vmax: Vmax after temperature adjustment from V25
% Tau_star: after temperature adjustment from 25 centi-degree
% Kc: after temperature adjustment from 25 centi-degree
% Ko: after temperature adjustment from 25 centi-degree
% PSII: after temperature adjustment from 25 centi-degree
% Phi: after temperature adjustment from 25 centi-degree 
% Rd: after temperature adjustment from 25 centi-degree
% Vomax: after temperature adjustment from 25 centi-degree

% or
% LPT.Jmax: Jmax after temperature adjustment from J25
% LPT.Vmax: Vmax after temperature adjustment from V25
% LPT.Tau_star: Tau_star after temperature adjustment from 25 centi-degree
% LPT.Kc: Kc after temperature adjustment from 25 centi-degree
% LPT.Ko: Ko after temperature adjustment from 25 centi-degree
% LPT.PSII: PSII after temperature adjustment from 25 centi-degree
% LPT.Phi: Phi after temperature adjustment from 25 centi-degree 
% LPT.Rd: Rd after temperature adjustment from 25 centi-degree
% LPT.Vomax: Vomax after temperature adjustment from 25 centi-degree

%% Key reference used for temperatuer scaling function includes:
% Bernacchi et al. 2013; June et al. 2004; Bernacchi et al. 2003; Bernacchi
% et al. 2001

% Specifically, Jmax temperature response function follows Jun et al. 2004
% Others follow Bernacchi et al. 2003


%% Temperature scaling function
Jmax_25=Func_Temperature_June(J25, Topt, 25);
Jmax=Func_Temperature_June(J25, Topt, T); % reference: June et al. 2004; Bernacchi et al. 2013
Jmax=Jmax./Jmax_25*J25;


% Temperature functions for Vcmax, Tau_star, Ko, Kc, Rd
% Vcmax
delta_H=65.33;
c=26.35; 
s1=Func_Temperature_Bernacchi(delta_H, c, T); % reference: Bernacchi et al. 2001; Bernacchi et al. 2013
sl_25=Func_Temperature_Bernacchi(delta_H, c, 25); % reference: Bernacchi et al. 2001; Bernacchi et al. 2013
Vmax=V25.*s1./sl_25; clear delta_H c s1 sl_25

% Gama_star; here I called it Tau_star 
delta_H=37.83;
c=19.02; 
s1=Func_Temperature_Bernacchi(delta_H, c, T); % reference: Bernacchi et al. 2001; Bernacchi et al. 2013
sl_25=Func_Temperature_Bernacchi(delta_H, c, 25); % reference: Bernacchi et al. 2001; Bernacchi et al. 2013
Tau_star=42.75.*s1./sl_25; clear delta_H c s1 sl_25 % umol/mol

% Kc
delta_H=79.43;
c=38.05; 
s1=Func_Temperature_Bernacchi(delta_H, c, T); % reference: Bernacchi et al. 2001; Bernacchi et al. 2013
sl_25=Func_Temperature_Bernacchi(delta_H, c, 25); % reference: Bernacchi et al. 2001; Bernacchi et al. 2013
Kc=404.9.*s1./sl_25; clear delta_H c s1 sl_25 % umol/mol

% Ko
delta_H=36.38;
c=20.30; 
s1=Func_Temperature_Bernacchi(delta_H, c, T); % reference: Bernacchi et al. 2001; Bernacchi et al. 2013
sl_25=Func_Temperature_Bernacchi(delta_H, c, 25); % reference: Bernacchi et al. 2001; Bernacchi et al. 2013
Ko=278.4.*s1./sl_25; clear delta_H c s1 sl_25 % mmol/mol

% Theta, curvature factor, which is noted as Phi here
% The Theta might be over-estimated by using Bernacchi et al. 2003 approach
Phi=0.76+0.018*T-3.7*(1e-4)*T.^2; % growth temperature=25 from Bernacchi et al. 2003 % Curvature factor between An-PAR
% Phi=0.54+0.023*T-3.3*(1e-4)*T.^2; % growth temperature=35 from Bernacchi et al. 2003 % Curvature factor between An-PAR

% Maximum quantumn yield, light adapated
PSII=0.352+0.022*T-3.4*(1e-4)*T.^2; % maximum quantumn yield; using light adapted version from Bernacchi et al. 2003

% Rd
delta_H=46.39;
c=18.72; 
s1=Func_Temperature_Bernacchi(delta_H, c, T); % reference: Bernacchi et al. 2001; Bernacchi et al. 2013
sl_25=Func_Temperature_Bernacchi(delta_H, c, 25); % reference: Bernacchi et al. 2001; Bernacchi et al. 2013
Rd=0.015*V25.*s1./sl_25; clear delta_H c s1 sl_25 % umol/mol

% Oxygen concentration
O=210.*Press./(101325); % mmol/mol

% not very sure for Vomax modeling
%Vomax=Vmax.*Ko.*Tau_star./(0.5*Kc.*O); % according to equation 7 in Bernacchi et al. 2001
Vomax25=V25.*278.4.*42.75./(0.5.*404.9.*O); % This is not confident

% Vomax_T
delta_H=60.11;
c=22.98; 
s1=Func_Temperature_Bernacchi(delta_H, c, T); % reference: Bernacchi et al. 2001; Bernacchi et al. 2013
sl_25=Func_Temperature_Bernacchi(delta_H, c, 25); % reference: Bernacchi et al. 2001; Bernacchi et al. 2013
Vomax=Vomax25.*s1./sl_25; clear delta_H c s1 sl_25 % umol/mol


LPT.Jmax=Jmax; %Jmax after temperature adjustment from J25
LPT.Vmax=Vmax; %Vmax after temperature adjustment from V25
LPT.Tau_star=Tau_star;  %Tau_star after temperature adjustment from 25 centi-degree
LPT.Kc=Kc;  %Kc after temperature adjustment from 25 centi-degree
LPT.Ko=Ko;  %Ko after temperature adjustment from 25 centi-degree
LPT.PSII=PSII;  %PSII after temperature adjustment from 25 centi-degree
LPT.Phi=Phi;  %Phi after temperature adjustment from 25 centi-degree 
LPT.Rd=Rd;    % Rd after temperature adjustment from 25 centi-degree
LPT.Vomax=Vomax; %Vomax after temperature adjustment from 25 centi-degree
