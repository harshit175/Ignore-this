
function [ Pc ] = Pcfunction(Tstep,alphastep,mlump, para )
%PCFUNCTION Summary of this function goes here
%   Detailed explanation goes here
%totalnodes =(m+1)*(n+1);
 %para means parameter
 
global totalnodes n1 n3 n m bb multfactor

% for i=1:n1
%      mlump(1,i)=1%i;
%      mlump(2,i)=1%3*i;
%      mlump(3,i)=1%3*i;
%      mlump(4,i)=1%i;
%  end
% mlump
%para

%bb

Pc=zeros(totalnodes,1);
ec=exp(1);
pc=0;

%k01=ec^20.2;
%k02=ec^9.5;

k01=0;
k02=ec^44.22;

Eatt1=0;
Eatt2=127038;

mm=0.1127;
nn=1.6563;

p=0.02184;
q=-0.3107;

R=8.314;

lumpingcounter=0; 
for c=n3+1:n1+n3+1
    
%     if (c>=n3+1 && c<=n1+n3)  %will have to be modified
%         mat=1; channel=1;
%     else
%         mat=2; channel=0;
%     end
    for r=1:m+1
        %     if mat==1
        %     rleft=rinm(mat)+(c-1)*deltar(mat);
        %     else
        %         rleft=rinm(mat)+(c-n1-1)*deltar(mat);
        %     end
        %     rright=rleft+deltar(mat);
        %     zbot=(r-1)*deltaz;
        %     ztop=r*deltaz;
        %
        a=(r-1)*(n+1)+c;  %this is the bottom left node number related to element number
        
        if r==1
            
            lumpingcounter=lumpingcounter+1;
            %mlumpelement=mlump(:,lumpingcounter);
            %     r
            %     c
            
        end
        
%         if ((r==1 || r==m+1))
%             if (c==n3+1 )
%                 multfactor=mlump(1,lumpingcounter);
%             elseif(c==n3+n1+1)
%                 multfactor=mlump(2,lumpingcounter-1);
%             else
%                 multfactor=mlump(1,lumpingcounter)+mlump(2,lumpingcounter-1);
%             end
%             
%         else
%             if (c==n3+1 )
%                 multfactor=2*mlump(1,lumpingcounter);
%             elseif(c==n3+n1+1)
%                 multfactor=2*mlump(2,lumpingcounter-1);
%             else
%                 
%                 multfactor=2*(mlump(1,lumpingcounter)+mlump(2,lumpingcounter-1));
%             end
%             
%         end
        
        
        
        for i=1:1
            k1=k01*ec^(-Eatt1/(R*(Tstep(a)+273)));
            k2=k02*ec^(-Eatt2/(R*(Tstep(a)+273)));
            
            dk1bydT=k1*(Eatt1/R)*(1/(Tstep(a)+273))^2;
            dk2bydT=k2*(Eatt2/R)*(1/(Tstep(a)+273))^2;
            
            if Tstep(a)<60
                alphmax=p*Tstep(a)+q;
                if alphmax-alphastep(a)<0
                    alphmax=alphastep(a);
                    %disp('curing value had to be adjusted')
                    %alphastep(a)
                    %Tstep(a)
                    %disp('row number')
                    %r
                    %disp('column number')
                    %c
                end
            else
                alphmax=1;
            end
            if para==1
            %pc(i)=mlump(i)*(1/60)*(k1+k2*(alphastep(a))^mm)*(alphmax-alphastep(a))^nn;
            pc(i)=(k1+k2*(alphastep(a))^mm)*(alphmax-alphastep(a))^nn;
            
            elseif para==2
                pc(i)=(dk1bydT+dk2bydT*(alphastep(a))^mm)*(alphmax-alphastep(a))^nn;
                %pc(i)=dpcbydT(i);
                
                if alphmax<1 %this 'if' statement accounts for temperature dependence of alphamax, 
                    %by adding the extra derivativeterm
                    pc(i)=pc(i) + (k1+k2*(alphastep(a))^mm)*nn*p*(alphmax-alphastep(a))^(nn-1);
                end
            
            
            else
                pc(i)=((k2*mm*alphastep(a)^(mm-1))*(alphmax-alphastep(a))^nn - nn*(k1+k2*(alphastep(a))^mm)*(alphmax-alphastep(a))^(nn-1));
                %pc(i)=dpcbydalpha(i);
            end
            % if not(isreal(pc(i)))
            %     Disp('Improve definition of Pcfunction')
            % end
        end
        
        if (r==1 || r==m+1)
            Pc(a)=multfactor(lumpingcounter)*pc; 
        else
            Pc(a)=2*multfactor(lumpingcounter)*pc;
        end
        
        %Pc(bb)=0; 
        
    end
end

end
%lumpingcounter


