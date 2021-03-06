

function z=MOP2(x)

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