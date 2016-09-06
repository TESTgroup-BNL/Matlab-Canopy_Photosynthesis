function Stat=Estimate_Rsqure(yobs, ypr)
%% Case 1
% yobs=yobs./max(yobs);
% yobs=yobs-mean(yobs);
% 
% ypr(:,1)=ypr(:,1)./max(ypr(:,1));
% ypr(:,1)=ypr(:,1)-mean(ypr(:,1));
% 
% ypr(:,2)=ypr(:,2)./max(ypr(:,2));
% ypr(:,2)=ypr(:,2)-mean(ypr(:,2));
% 
% ypr(:,3)=ypr(:,3)./max(ypr(:,3));
% ypr(:,3)=ypr(:,3)-mean(ypr(:,3));

%% Case 2
yobs=yobs./mean(yobs);

ypr(:,1)=ypr(:,1)./mean(ypr(:,1));
ypr(:,2)=ypr(:,2)./mean(ypr(:,2));
ypr(:,3)=ypr(:,3)./mean(ypr(:,3));

SST=sum((yobs-mean(yobs)).^2);
SS1=sum((yobs-ypr(:,1)).^2);
SS2=sum((yobs-ypr(:,2)).^2);
SS3=sum((yobs-ypr(:,3)).^2);

R2_1=1-SS1./SST;
R2_2=1-SS2./SST;
R2_3=1-SS3./SST;

[R,P1] = corrcoef(yobs,ypr(:,1));
[R,P2] = corrcoef(yobs,ypr(:,2));
[R,P3] = corrcoef(yobs,ypr(:,3));

Stat=[R2_1 P1(1,2) R2_2 P2(1,2) R2_3 P3(1,2)];