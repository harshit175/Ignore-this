
clear Uit

Uit(:,1)=(U);


P_pre=sparse(vertcat(F', Pc_pre));
gamma=B*U+(1-beta)*P_pre;
clear Residue

for iteration=1:1000
    iteration
    
    Pc(:,iteration)=Pcfunction(Uit(i,iteration),Uit(j,iteration),mlump,1);
    Pit=(vertcat(F', Pc(:,iteration))); %it means iteration
    
    
    Residue(:,iteration)=A*Uit(:,iteration)-beta*Pit-gamma; 
    
    if check==0 
    Residue(bb,iteration)=0; %Sets temperature constant 
    Residue(totalnodes+bb,iteration)=0; %Sets cure constant
    end
    
    
    if iteration==1
        firstitnorm=norm(Residue(i,1));
    end
    if norm(Residue(i,iteration))<.0001*firstitnorm
        break
    end
    %Im21=-sparse(beta*diag(Pcfunction(Uit(i,iteration),Uit(j,iteration),mlump,2)));
    Im21=-beta*sparse(1:totalnodes,1:totalnodes,(Pcfunction(Uit(i,iteration),Uit(j,iteration),mlump,2)));
    if check==0 && exist('bb')
        Im21=modifymatrices(Im21,bb);
    end
    
    %Im22=sparse(Cc-beta*diag(Pcfunction(Uit(i,iteration),Uit(j,iteration),mlump,3)));
    Im22=Cc-beta*sparse(1:totalnodes,1:totalnodes,(Pcfunction(Uit(i,iteration),Uit(j,iteration),mlump,3)));
    
    %tt2=cputime;
    Im=[Im1; Im21 Im22];
    %timeelapsed2=cputime-tt2
    
    %code to impose constant Temperature starts
    if check==0    
    Im=modifymatrices(Im,bb);
    end
    %code to impose constant Temperature ends
    
    
    %tic
    %disp('dU calculation starts')
    dU=-(Im\Residue(:,iteration));  %this is the most time consuming line in the whole process !!
    %disp('dU calculation ends')
    %toc
    
    
    %dUrecord(:,iteration,step)=dU;
    Uit(:,iteration+1)=Uit(:,iteration)+dU;    
end



%Pc_implicit=Pcfunction(Uit(i,iteration+1),Uit(j,iteration+1),mlump,1);
%above line is not needed

