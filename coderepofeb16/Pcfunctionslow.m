
function [ Pc ] = Pcfunction(Tstep,alphastep,mlump, para )
%PCFUNCTION Summary of this function goes here
%   Detailed explanation goes here
%totalnodes =(m+1)*(n+1);
 
global totalnodes n1 n3 n m bb 
% for i=1:n1
%      mlump(1,i)=i;
%      mlump(2,i)=3*i;
%      mlump(3,i)=3*i;
%      mlump(4,i)=i;
% end
%  mlump
 
%para

%bb

Pc=zeros(totalnodes,1);
ec=exp(1);
pc=[0 0 0 0];

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
for c=n3+1:n1+n3
%     if (c>=n3+1 && c<=n1+n3)  %will have to be modified
%         mat=1; channel=1;
%     else
%         mat=2; channel=0;
%     end 
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
    a=[(r-1)*(n+1)+c (r-1)*(n+1)+c+1 r*(n+1)+c+1 r*(n+1)+c];  %These are 
    %the element node numbers as a function of row number, column number 
    %and total elements in one row
    if r==1
    mlumpelement=mlump(:,lumpingcounter);
    %mlumpelement
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

dk1bydT=k1*(Eatt1/R)*(1/(Tstep(a(i))+273))^2;
dk2bydT=k2*(Eatt2/R)*(1/(Tstep(a(i))+273))^2;

if Tstep(a(i))<78
    alphmax=p*Tstep(a(i))+q;
    if alphmax-alphastep(a(i))<0
        alphmax=alphastep(a(i));
        %disp('curing value had to be adjusted')
        %alphastep(a(i))
        %Tstep(a(i))
        %disp('row number')
        %r
        %disp('column number')
        %c
    end
else
    alphmax=1;
end

%pc(i)=mlump(i)*(1/60)*(k1+k2*(alphastep(a(i)))^mm)*(alphmax-alphastep(a(i)))^nn;
            if para==1
            %pc(i)=mlump(i)*(1/60)*(k1+k2*(alphastep(a(i)))^mm)*(alphmax-alphastep(a(i)))^nn;
            pc(i)=(1/60)*(k1+k2*(alphastep(a(i)))^mm)*(alphmax-alphastep(a(i)))^nn;
            
            elseif para==2
                pc(i)=(1/60)*(dk1bydT+dk2bydT*(alphastep(a(i)))^mm)*(alphmax-alphastep(a(i)))^nn;
                %pc(i)=dpcbydT(i);
            
            
            else
                pc(i)=(1/60)*((k2*mm*alphastep(a(i))^(mm-1))*(alphmax-alphastep(a(i)))^nn - nn*(k1+k2*(alphastep(a(i)))^mm)*(alphmax-alphastep(a(i)))^(nn-1));
                %pc(i)=dpcbydalpha(i);
            end
% if not(isreal(pc(i)))
%     Disp('Improve definition of Pcfunction')
% end
end

for i=1:4
       Pc(a(i))=Pc(a(i))+mlumpelement(i)*pc(i); 
end   

    Pc(bb)=0;  %don't forget to change when running actual code
end
end
%lumpingcounter 
disp('Pcfunction executed')
end

