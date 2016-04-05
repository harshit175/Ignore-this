ec=exp(1);

T=50;

alpha=0.1;

k01=ec^20.2;
k02=ec^9.5;

Eatt1=64735;
Eatt2=31561;

mm=1;
nn=1.07;

p=6.96e-3;
q=4.57e-1;

R=8.314;

k1=k01*ec^(-Eatt1/(R*(T+273)));
k2=k02*ec^(-Eatt2/(R*(T+273)));

dk1bydT=k1*(Eatt1/R)*(1/(T+273))^2;
dk2bydT=k2*(Eatt2/R)*(1/(T+273))^2;

pc=(k1+k2*alpha^mm)*(1-alpha)^nn;

dpcbydT=(dk1bydT+dk2bydT*(alpha)^mm)*(1-alpha)^nn;

dpcbydalpha=(k2*mm*alpha^(mm-1))*(1-alpha)^nn - nn*(k1+k2*(alpha)^mm)*(1-alpha)^(nn-1);

dpcbydT/pc

dpcbydalpha/pc