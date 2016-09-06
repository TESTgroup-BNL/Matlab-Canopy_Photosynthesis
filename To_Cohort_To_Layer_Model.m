function JW=To_Cohort_To_Layer_Model(LAI, Lage, PC_age, ftop, LAI_cut, Flag_Scale, Vcmax0_25)
% LAI--Phenology of LAI
% Lage--Phenology of leaf age
% PC_age--Age dependency of Vcmax
% ftop--the fraction of leaves allocated to top canopy
% Vcmax0_25=40; % Bonan et al., 2012 for the tropcis
% Flag_Scale=0; % 1--apply scale factor; 0--do not apply; to scale Vcmax into a given value

% LAI=0.85*LAI;

for i=1:12
    YMO_tot(i,:)=LAI(i,1)*Lage(i,:); % total turnover leaf
end

YMO_top=YMO_tot*ftop;
YMO_bottom=YMO_tot-YMO_top;

YMO_top(:,3)=LAI_cut-YMO_top(:,1)-YMO_top(:,2);
YMO_bottom(:,3)=LAI-LAI_cut;


for i=1:12
    if YMO_top(i,1)+YMO_top(i,2)>LAI_cut % leaf life span< 1 year
       V25_top(i,1)=(YMO_top(i,1)*PC_age(1,1)+YMO_top(i,2)*PC_age(1,2))./(YMO_top(i,1)+YMO_top(i,2)); 
       YMO_bottom(i,3)=LAI(i,1)-YMO_top(i,1)-YMO_top(i,2);       
    else
       V25_top(i,1)=(YMO_top(i,1)*PC_age(1,1)+YMO_top(i,2)*PC_age(1,2)+YMO_top(i,3)*PC_age(1,3))./(YMO_top(i,1)+YMO_top(i,2)+YMO_top(i,3)); 
    end
    
    V25_bottom(i,1)=(YMO_bottom(i,1)*PC_age(1,1)+YMO_bottom(i,2)*PC_age(1,2)+YMO_bottom(i,3)*PC_age(1,3))./(YMO_bottom(i,1)+YMO_bottom(i,2)+YMO_bottom(i,3)); 
end

if Flag_Scale==1
    V25_top=V25_top./34*Vcmax0_25;
    V25_bottom=V25_bottom./34*Vcmax0_25;
end

JW.V25_top=V25_top;
JW.V25_bottom=V25_bottom;