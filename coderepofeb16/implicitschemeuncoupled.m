clear Uit

Uit(:,1)=sparse(U);


%P_pre=sparse(vertcat(F', Pc_pre));
%gamma=B*U+(1-beta)*P_pre;
clear Residue
for iteration=1:1000
    iteration
    %Pc(:,iteration)=Pcfunction(Uit(i,iteration),Uit(j,iteration),mlump,1);
    %Pit=sparse(vertcat(F', Pc(:,iteration))); %it means iteration
    
    Residue(:,iteration)=Im11*Uit(i,iteration)-B11*U(i)-rho1*Hr*Pc_pre-F';
    
    if iteration==1
        firstitnorm=norm(Residue(1:totalnodes,1));
    end
    if norm(Residue(1:totalnodes,iteration))<.0001*firstitnorm
        break
    end    
    
    %Im=[Im1; Im21 Im22];
    disp('dU calculation starts')
    dU(i)=-Im11\Residue(:,iteration);
    disp('dU calculation ends')
    dUrecord(i,iteration,step)=dU(i);
    Uit(i,iteration+1)=Uit(i,iteration)+dU(i)';    
end

for iterationc=1:1000
    iterationc
    Pc(:,iterationc)=Pcfunction(Uit(i,iteration),Uit(j,iterationc),mlump,1);
    Residuec(:,iterationc)=A22*(Uit(j,iteration)-U(j))-beta*Pc(:,iterationc)-(1-beta)*Pc_pre;
    %have the temperature from above iteration, curing from previous step.
    %iterate for curing until residue goes to a small value.
    
    
    if iterationc==1
        firstitnormc=norm(Residuec(:,1));
    end
    if norm(Residuec(1:totalnodes,iterationc))<.0001*firstitnormc
        break
    end
    
    %Im21=-beta*sparse(1:totalnodes,1:totalnodes,(Pcfunction(Uit(i,iteration),Uit(j,iteration),mlump,2)));
    Im22=Cc-beta*sparse(1:totalnodes,1:totalnodes,(Pcfunction(Uit(i,iteration),Uit(j,iterationc),mlump,3)));
    
    disp('dUc calculation starts')
    dU(j)=-Im22\Residuec(:,iterationc);
    disp('dUc calculation ends')
    dUrecord(j,iteration,step)=dU(j);
    Uit(j,iterationc+1)=Uit(j,iterationc)+dU(j)'; 
    
end

%Pc_implicit=Pcfunction(Uit(i,iteration+1),Uit(j,iteration+1),mlump,1);
%above line is not needed

