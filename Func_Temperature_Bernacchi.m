function scale=Func_Temperature_Bernacchi(delta_H, c, T)
% temperature response function and parameters follow the reference at
% Bernacchi et al., 2001

% delta_H: activation energey
% c: scaling constant
Tk=T+273.15;
R=8.3144598; % J/K/mol, universal gas constant
delta_H=delta_H*1000; % unit conversion from kJ/mol to J/mol
scale=exp(c-delta_H./(R.*Tk));