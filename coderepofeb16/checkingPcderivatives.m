clear all

%PCFUNCTION Summary of this function goes here
%   Detailed explanation goes here
%totalnodes =(m+1)*(n+1);

%global totalnodes n1 n3 n m bb 

%para

%bb
epsilon=1e-5;
Tstep(1)=200;
Tstepadd=epsilon*Tstep(1);
Tstep(2)=Tstep(1)+Tstepadd;

alphastep(1)=.995;
alphastepadd=0;%epsilon*alphastep(1);
alphastep(2)=alphastep(1)+alphastepadd;
ec=exp(1);


%k01=ec^20.2;
%k02=ec^9.5;

k01=600*ec^20.2;
k02=600*ec^9.5;

Eatt1=64735;
Eatt2=31561;

mm=1;
nn=1.07;

p=6.96e-3;
q=4.57e-1;

R=8.314;

lumpingcounter=1; mat=1;

alphmax=1;    
   
%     if mat==1
%     rleft=rinm(mat)+(c-1)*deltar(mat);
%     else
%         rleft=rinm(mat)+(c-n1-1)*deltar(mat);
%     end
%     rright=rleft+deltar(mat);
%     zbot=(r-1)*deltaz;
%     ztop=r*deltaz;   
%     
    
    %pc=channel*0e-3*[1;1;1;1]; %this line will be replaced by a detailed integration rule 
    %to calculate pc.
    

for i=1:2

k1=k01*ec^(-Eatt1/(R*(Tstep(i)+273)));
k2=k02*ec^(-Eatt2/(R*(Tstep(i)+273)));

dk1bydT=k1*(Eatt1/R)*(1/(Tstep(i)+273))^2;
dk2bydT=k2*(Eatt2/R)*(1/(Tstep(i)+273))^2;



%pc(i)=mlump(i)*(1/60)*(k1+k2*(alphastep(i)(a(i)))^mm)*(alphmax-alphastep(i)(a(i)))^nn;
pc(i)=(1/60)*(k1+k2*(alphastep(i))^mm)*(alphmax-alphastep(i))^nn;


dpcbydT(i)=(1/60)*(dk1bydT+dk2bydT*(alphastep(i))^mm)*(alphmax-alphastep(i))^nn;




dpcbydalpha(i)=(1/60)*((k2*mm*alphastep(i)^(mm-1))*(alphmax-alphastep(i))^nn - nn*(k1+k2*(alphastep(i))^mm)*(alphmax-alphastep(i))^(nn-1));


% if not(isreal(pc(i)))
%     Disp('Improve definition of Pcfunction')
% end

end
 


%lumpingcounter


