Tstep=80;

ec=exp(1);

k01=600*ec^20.2;
k02=600*ec^9.5;

Eatt1=64735;
Eatt2=31561;

mm=1;
nn=1.07;

p=6.96e-3;
q=4.57e-1;

R=8.314;

k1=k01*ec^(-Eatt1/(R*(Tstep+273)));
k2=k02*ec^(-Eatt2/(R*(Tstep+273)));

dpcbydalpha(i)=(1/60)*((k2*mm*alphastep(a(i))^(mm-1))*(alphmax-alphastep(a(i)))^nn - nn*(k1+k2*(alphastep(a(i)))^mm)*(alphmax-alphastep(a(i)))^(nn-1));


dU=-inv(Im)*R(:,iteration);
dUrecord(:,iteration,step)=dU;
Uit(:,iteration+1)=Uit(:,iteration)+dU;
iteration