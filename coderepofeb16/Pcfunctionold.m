function [ Pc ] = Pcfunction(n1,n,m,totalnodes,Tstep,alphastep,mlump )
%PCFUNCTION Summary of this function goes here
%   Detailed explanation goes here
Pc=zeros(totalnodes,1);
ec=exp(1);
pc=[0 0 0 0];

k01=ec^20.2;
k02=ec^9.5;

Eatt1=64735;
Eatt2=31561;

mm=1;
nn=1.07;

p=6.96e-3;
q=4.57e-1;

R=8.314;

lumpingcounter=1; mat=1;
for c=1:n1
    if (c<=n1)  %will have to be modified
        mat=1; channel=1;
    else
        mat=2; channel=0;
    end 
for r=1:m    
%     if mat==1
%     rleft=rinm(mat)+(c-1)*deltar(mat);
%     else
%         rleft=rinm(mat)+(c-n1-1)*deltar(mat);
%     end
%     rright=rleft+deltar(mat);
%     zbot=(r-1)*deltaz;
%     ztop=r*deltaz;   
%     
    a=[(r-1)*(n+1)+c (r-1)*(n+1)+c+1 r*(n+1)+c r*(n+1)+c+1];
    if r==1
    mlumpelement=mlump(:,lumpingcounter);
    lumpingcounter=lumpingcounter+1;
%     r
%     c
%     lumpingcounter
    end
    %pc=channel*0e-3*[1;1;1;1]; %this line will be replaced by a detailed integration rule 
    %to calculate pc.
    


for i=1:4
k1=k01*ec^(-Eatt1/(R*(Tstep(a(i))+273)));
k2=k02*ec^(-Eatt2/(R*(Tstep(a(i))+273)));
if Tstep(a(i))<78
    alphmax=p*Tstep(a(i))+q;
    if alphmax-alphastep(a(i))<0
        alphmax=alphastep(a(i));
        disp('curing value had to be adjusted')
        alphastep(a(i))
        Tstep(a(i))
        disp('row number')
        r
        disp('column number')
        c
    end
else
    alphmax=1;
end

%pc(i)=mlump(i)*(1/60)*(k1+k2*(alphastep(a(i)))^mm)*(alphmax-alphastep(a(i)))^nn;
pc(i)=mlumpelement(i)*(1/60)*(k1+k2*(alphastep(a(i)))^mm)*(alphmax-alphastep(a(i)))^nn;
% if not(isreal(pc(i)))
%     Disp('Improve definition of Pcfunction')
% end
end
for i=1:4
       Pc(a(i))=Pc(a(i))+pc(i); 
end    
    
end
end
%lumpingcounter
end

