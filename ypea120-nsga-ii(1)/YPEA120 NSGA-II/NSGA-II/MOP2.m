%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA120
% Project Title: Non-dominated Sorting Genetic Algorithm II (NSGA-II)
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

function z=MOP2(x)

%     n=numel(x);
%     
%     z1=1-exp(-sum((x-1/sqrt(n)).^2));
%     
%     z2=1-exp(-sum((x+1/sqrt(n)).^2));
%     
% %     z=[z1 z2]';
% v=x(1);
% f=x(2);
% ap=x(3);
% Kt=1;Ct=242;m=2.5;q=2.25;s=0.75;Nm=0.98;
% Cfz=270;Cfy=199;Cfx=294;
% Kfz=0.75;Kfy=1.35;Kfx=1.0;
% Xfz=1.0;Yfz=0.75;Nfz=-0.15;
% Xfx=0.9;Yfx=0.6;Nfx=-0.3;
% dw=97;lw=35;delta=0.2;Tct=2;Tot=5;T=120;Ra=1.25;r=0.3;
% a1=(pi*dw*lw*delta)/(1000*v*f*ap);
% a2=Tct/(Ct*Kt*power(v,-m)*power(f,-q)*power(ap,-s))
% z1=a1*(1+a2)+Tot
% z2=[];
%  z2(1)=Cfz*power(ap,Xfz)*power(f,Yfz)*power(v,Nfz)*Kfz-Nm*Pm;
%  z2(2)=power(v,m)*power(f,q)-(Kt*Ct)/T*power(ap,s);
%  z2(3)=f-power(8*Ra*r,0.5);
n=7;S=300;Sn=[48 50 48 45 49 46 45 ];hn=[1 1 1 1 1 1 1];
Dn=[1000 800 600 1200 700 800 900];
MOQ=2500;
T=x(1);
K=x(2:n+1);
z1=S/T+sum(Sn(1:n)./(K*T))+(0.5)*sum(K*T.*Dn.*hn);

% I=[0 0 1 0 1 1 0];
I=[1 0 1 0 1 1 1];

z2=-(sum(Dn*T.*I.*K)-MOQ);
z3=-(I-0.5).*(1.5-K);
z3=sum((z3>0).*z3);
z4=[z2,z3];
err=sum((z4>0).*z4);
z=[z1,err];
end