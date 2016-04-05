clear all

Tstep=200;

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
alphmax=p*Tstep+q;
if alphmax>.99
    alphmax=1;
end
beta=0.5;
dt=.2e-3;
totaltime=.4;
te=int64(totaltime/dt)
pcnr(te)=0;
alphastep(1)=0.001;
check=0;
alphmax=1;  %this line is to be removed
for i=2:double(te)
    
    pc(i-1)=(1/60)*(k1+k2*(alphastep(i-1))^mm)*(alphmax-alphastep(i-1))^nn;
    
    if beta~=0
        alphai(1)=alphastep(i-1);
        iteration=1;
        clear Residue
        Residue(1)=1;
        for iteration=1:100 
            dpcbydalpha(iteration)=(1/60)*((k2*mm*alphai(iteration)^(mm-1))*(alphmax-alphai(iteration))^nn - nn*(k1+k2*(alphai(iteration))^mm)*(alphmax-alphai(iteration))^(nn-1));
            pci(iteration)=(1/60)*(k1+k2*(alphai(iteration))^mm)*(alphmax-alphai(iteration))^nn;
            Residue(iteration)=(alphai(iteration)-alphastep(i-1))/dt-beta*pci(iteration)-(1-beta)*pc(i-1);
            dalpha(iteration)=-Residue(iteration)/(1/dt-beta*dpcbydalpha(iteration));
            alphai(iteration+1)=alphai(iteration)+dalpha(iteration);
            if not(isreal(dpcbydalpha(iteration)))
                disp('cure value complex')
                check=1;                
                break
            end
            Residuematrix(iteration,i)=Residue(iteration);
            if abs(Residue(iteration))<.0001*abs(Residue(1))                
                break
            end
            
        end
        
        iteration
        dalphanr(i)=alphai(iteration+1)-alphastep(i-1);
        %pcnr(i)=(1/60)*(k1+k2*(alphai(iteration+1))^mm)*(alphmax-alphai(iteration+1))^nn; %nr is newton raphson
        alphastep(i)=alphai(iteration+1);
    else
        alphastep(i)=alphastep(i-1)+pc(i-1)*dt;
        dalphastep(i)=alphastep(i)-alphastep(i-1);
    end
    i
    if check==1
        break
    end

time(i)=(i-1)*dt;
% 
% 
% dU=-inv(Im)*R(:,iteration);
% dUrecord(:,iteration,step)=dU;
% Uit(:,iteration+1)=Uit(:,iteration)+dU;
% iteration
end

plot(time,alphastep,'k')

xlabel('Time (sec.)','fontsize',14)
ylabel('\alpha','fontsize',14)