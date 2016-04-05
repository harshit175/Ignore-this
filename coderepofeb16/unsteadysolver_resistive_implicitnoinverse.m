%Data from input file will be used. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Processing starts

if gpratio~=1
n2int=int8(log(1-(rout-rinter)/(rinter-rin)*(1-gpratio)*n1)/log(gpratio));
n2int=n2int+1;
n2=double(n2int);
for i=1:n2+1
    r2(i)=rinter+ (rinter-rin)/n1*(1-gpratio^(i-1))/(1-gpratio); %i=1, r2=rinter; i=2, r2=rinter+a;  i=3, r2=rinter+a+a*gpratio
    deltar=[(rinter-rin)/n1 (rout-rinter)/n2 (rin-rc)/n3];
end
else
n2=0;  %n2 as to be specified if gpratio =1; the code will be erroneous otherwise.
if n2==0
    deltar=[(rinter-rin)/n1 0 (rin-rc)/n3];
else
    deltar=[(rinter-rin)/n1 (rout-rinter)/n2 (rin-rc)/n3];
end
end

n=n1+n2+n3;


global totalnodes;
totalnodes =(m+1)*(n+1);
curingnodes=(m+1)*(n1+1);



rinm=[rin rinter rc]; %in stands for inner radius, m stands for material
%So rinm means inner radius of either polymer or PDMS. 
%For instance, innr radius for PDMS (mat=2) is the interface radius i.e. .
%rinter
deltaz=hz/m;

%Assembling KT and CT matrix
K=sparse(totalnodes,totalnodes);
CT=sparse(totalnodes,totalnodes);
Cc=sparse(totalnodes,totalnodes);

%mlump=diag(ctfunction(0,(rinter-rin)/n1,0,deltaz,1,1,1)); 

F(totalnodes)=0; %load vector initialised at zero, this will be modified later and used in the
%ignition stage
F0=F;  %this will be used after hot boundary is insulated


disp('matrix assembly starting')
i=1:4; lumpingcounter=1;
mlump=[0; 0; 0; 0];
for c=1:n
    if (c<=n3)
        mat=3; channel=0;
    elseif (c>=n3+1 && c<=n3+n1)
        mat=1; channel=1;
    else
        mat=2; channel=0;
    end
for r=1:m
    if mat==1
        rleft=rinm(mat)+(c-n3-1)*deltar(mat);
        rright=rleft+deltar(mat);
    elseif mat==2 %i.e. PDMS, has two cases, one for a gp in the r direction, and other is no gp
        if gpratio~=1
            rleft=r2(c-n1-n3);
            rright=r2(c-n1-n3+1);
        else
            rleft=rinm(mat)+(c-n1-n3-1)*deltar(mat);
            rright=rleft+deltar(mat);
        end
    elseif mat==3
        rleft=rinm(mat)+(c-1)*deltar(mat);
        rright=rleft+deltar(mat);
    end
    
    zbot=(r-1)*deltaz;
    ztop=r*deltaz;
    
    if r==1 %to be used only when mesh in z direction is uniform
    k=stiffnessfunction(rleft,rright,zbot,ztop,mat,kc);
    ct=rhocp(mat)*ctfunction(rleft,rright,zbot,ztop,mat,dt);    
    end
    
    a=[(r-1)*(n+1)+c (r-1)*(n+1)+c+1 r*(n+1)+c+1 r*(n+1)+c];
   
    K(a(i),a(i))=K(a(i),a(i))+k(i,i);
    CT(a(i),a(i))=CT(a(i),a(i))+ct(i,i);
    
   
    cc=channel*ctfunction(rleft,rright,zbot,ztop,mat,dt);
    Cc(a(i),a(i))=Cc(a(i),a(i))+cc(i,i);
    
    if mat==1 && r==1  %capitalizes on the structured nature of microchannel mesh, so is only computed for the first row of elements.
        mlump(:,lumpingcounter)=diag(ctfunction(rleft,rright,zbot,ztop,mat,1));
        lumpingcounter=lumpingcounter+1;        
    end
    
    if mat==3
        qhe=heatsourcevectorfunction( rleft,rright,zbot,ztop,Q);
        F(a(i))=F(a(i))+qhe(i);
    end
    
end
c/n
end

Fresistive=F;

lumpingcounter=0;  %calculation of multiplication factors for Pcfunction
for cc=n3+1:n1+n3+1
    lumpingcounter=lumpingcounter+1;
    if (cc==n3+1 )
        multfactor(cc-n3)=mlump(1,lumpingcounter);
    elseif(cc==n3+n1+1)
        multfactor(cc-n3)=mlump(2,lumpingcounter-1);
    else
        multfactor(cc-n3)=mlump(1,lumpingcounter)+mlump(2,lumpingcounter-1);
    end
end

%The zero diagonal elements have to be changed to 1, else the Cc matrix is
%singular, this does not change the solution as the alpha corresonding to
%those elements is in PDMS (mat=2)or copper (mat=3) and is a constant (equal to zero)
for i=1:totalnodes
    if (Cc(i,i)==0)
        Cc(i,i)=1;
    end
end
%Assembly of KT, CT and Cc done
disp('assembly of matrices complete')
%initialising alpha to a small value<<<1 at t=0
% i=1:4;
% for c=1:n
%     if (c>=n3+1 && c<=(n1+n3))
%         mat=1; curingnode=1;
%     elseif (c>=n1+n3+1)
%         mat=2; curingnode=0;
%     elseif (c<=n3)
%         mat=3; curingnode=0;
%     end
%     
% for r=1:m
%     a=[(r-1)*(n+1)+c (r-1)*(n+1)+c+1 r*(n+1)+c r*(n+1)+c+1];
%     alpha(a(i),1)=curingnode*0.000;
% end
% end
%alpha=alpha';
%initialization complete






%Assembling load vector (to ignite the monomer)
i=1:4;
for c=1:(n3+n1)
    if (c<=n3)
        mat=3;
        rleft=rinm(mat)+(c-1)*deltar(mat);
    elseif (c>=n3+1 && c<=n3+n1)
        mat=1; channel=1;
        rleft=rinm(mat)+(c-n3-1)*deltar(mat);
    end
    
    rright=rleft+deltar(mat);
    Qve=heatvectorfunction(rleft,rright,q);
    
    a=[c c+1 c+2 c+3];  %could give an error
    F(a(i))=F(a(i))+Qve(i);
    
    if (c<=n3)
        %F0(a(i))=F0(a(i))+Qve(i);  %This step is for continued heating
    end
end
Fresistiveplusignition=F;
Fignition=Fresistiveplusignition-Fresistive;
%F(c+1)=0;
%Load vector assembled


%Applying constant temperature boundary conditions on right boundary

j=1:(m+1);
global bb
bb=1:n1+n3+1;   %bb=ignition boundary nodes

if Q~=0
   clear bb 
end

%bb

br=(n+1)*j; %br=right boundary nodes
bl=br-n;    %bl=left boundary nodes


%Tl=000; %messing with this alone will give unphysical results, other changes to 
%left boundary have to be made as well.
%Tr=000;  %Temperatures at left and right boundaries defined

%DT=CT-K; %equation to be solved is CT*T(tn)=DT*T(t(n-1))+F;

%CT0=CT; DT0=DT; K0=K; Cc0=Cc; %these are needed after hot boundary is switched off
%disp('post ignition matrices created')

%Adjusting F values for right boundary
% for p=1:totalnodes
%     F(p)=F(p)+Tr*( sum(DT(p,br))-sum(CT(p,br)) ) + Tl*( sum(DT(p,bl))-sum(CT(p,bl)) );
%     F0(p)=F0(p)+Tr*( sum(DT(p,br))-sum(CT(p,br)) ) + Tl*( sum(DT(p,bl))-sum(CT(p,bl)) );
% end


nodesnearresist=horzcat(bl+n3,bl+n3+1,bl+n3+2);
nodesnearignition=n+2+n3:n+2+n1+n3;

%if (exist('bb'))
%Adjusting F values for bottom boundary
%for p=1:totalnodes
    %F(p)=F(p)+Tb*( sum(DT(p,bb))-sum(CT(p,bb)) );
    %F0(p)=F0(p)+Tb*( sum(DT(p,bb))-sum(CT(p,bb)) );
%end

%Setting nodes in F matrix to fixed temperature values corresponding to
%constant temperature boundary conditions
%F(br)=Tr;  
%F0(br)=Tr;
%F(bl)=Tl;
%F(bb)=Tb;

%Fss=F; %this is the force vector for steady state

%F(br)=0; %this is needed for solving in unsteady state for ignition on
%F0(br)=0; %this is needed for solving in unsteady state for ignition off

%F(bb)=0; %this is needed for solving in unsteady state for ignition on
%F0(bb)=0; %this is needed for solving in unsteady state for ignition off,
          %this step is not necessary, as F0 is already zero.

%CT=modifymatrices(CT,br); %see the function for explanation
%DT=modifymatrices(DT,br);
%K=modifymatrices(K,br);

%CT=modifymatrices(CT,bl); %see the function for explanation
%DT=modifymatrices(DT,bl);
%K=modifymatrices(K,bl);

%CT=modifymatrices(CT,bb); %see the function for explanation
%DT=modifymatrices(DT,bb);
%K=modifymatrices(K,bb);
%Cc=modifymatrices(Cc,bb);  
%disp('matrices modified for bottom boundary')

%end



%initialising T vector for t=0;

%T(br,1)=Tr;
%T(bl,1)=Tl;


%solving for steady state
%Tss=K\Fss';

%alpha=zeros(totalnodes,1)+.0001;
Pc_pre=zeros(totalnodes,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%solving for unsteady state

%implicit solver

A11=sparse(CT+beta*K);  %has to be modified for pre ignition

A12=sparse(-rho1*Hr*Cc);

A21=sparse(totalnodes,totalnodes);

A22=sparse(Cc);

A=sparse([A11 A12; A21 A22]);  %has to be modified for pre ignition
disp('pre ignition matrices inversion started')
%invA=inv(A);




B11=sparse(CT-(1-beta)*K);  %has to be modified for pre ignition

B12=sparse(-rho1*Hr*Cc);

B21=sparse(totalnodes,totalnodes);

B22=sparse(Cc);

B=[B11 B12; B21 B22];  %has to be modified for pre ignition

B011=B11;%since B11 will be modified for pre ignition
B0=B;%since B will be modified for pre ignition

%clear CT  K


%IA=speye(size(A));
%invA = A\IA;
%invA0=invA;  %since it will be modified for pre ignition
A0=A;        %since it will be modified for pre ignition
A011=A11;    %since it will be modified for pre ignition


%invAB0=invA0*B0;%since it will be modified for pre ignition

%invAB=invA*B;%since it will be modified for pre ignition

T(totalnodes,te/recordstep)=0; alpha(totalnodes,te/recordstep)=0;


T(1:totalnodes,1)=Tambient;

Tinitialrows=double(int64(m/10));

for count=1:Tinitialrows
   %T(count*(n+1)+1:count*(n+1)+n3+1)=Tambient+(Tinitialrows-count)*(Tb-Tambient)/Tinitialrows; 
end

if (1==0)%exist('bb')
    for p=1:totalnodes
        %F(p)=F(p)+Tb*( sum(B11(p,bb))-sum(A11(p,bb)) ); %needs to be commented if ignition 
        %is to be done using heating and not temperature bc
        
    end
    
    %F(bb)=0; %needs to be commented if ignition is to be done using heating and not temperature bc
    
    %T(bb,1)=Tb;  %needs to be commented if ignition is to be done using heating and not temperature bc
    
    
    %nodesnearignition=n+2+n3;
        
    A11=modifymatrices(A11,bb);
    %A12=modifymatrices(A12,bb);
    %A22=modifymatrices(A22,bb);
    
    A=sparse([A11 A12; A21 A22]);
    disp('post ignition matrices inversion started') %isn't this pre ignition?
    %invA=inv(A);
    %invA = A\IA;
    B11=modifymatrices(B11,bb);
    %B12=modifymatrices(B12,bb);
    %B22=modifymatrices(B22,bb);
    
    B=[B11 B12; B21 B22];
    
    %invAB=invA*B;
end



% for step=2:te
%     T(:,step)=inv(CT)*DT*T(:,step-1)+inv(CT)*F';
% end

%A=inv([sparse(CT) sparse(-rho1*Hr*Cc); sparse(totalnodes,totalnodes) sparse(Cc)]);
%B= [sparse(DT) sparse(-rho1*Hr*Cc); sparse(totalnodes,totalnodes) sparse(Cc)];
%AB=A*B;


%A0=inv(sparse([CT0 -rho1*Hr*Cc0; sparse(totalnodes,totalnodes) sparse(Cc0)]));
%B0=sparse([DT0 -rho1*Hr*Cc0; sparse(totalnodes,totalnodes) sparse(Cc0)]);
%AB0=A0*B0;


i=1:totalnodes;
j=totalnodes+1:2*totalnodes;
%ignitionnodes=1:n1;



Im1=sparse([A11 A12]);
Im01=sparse([A011 A12]);

U=vertcat(T(:,1), alpha(:,1));

check=0;
%Pc=zeros(totalnodes,1);
disp('time marching starting')

for step=2:te
tic    
Pc_pre=Pcfunction(U(i),U(j),mlump,1);
    
    if beta~=0
    disp('entering implicit algorithm')
    implicitscheme%uncoupled
    disp('exiting implicit algorithm')
    U(i)=Uit(i,iteration); U(j)=Uit(j,iteration);
    else        
    U = A\(B*U)+A\sparse(vertcat(F', Pc_pre));
    end
    if not(isreal(U))
        disp('curing values imaginary at time step');
        disp(step)
        complexcurestep=step;
        break
    end
    %T(:,step)=U(i);
    %alpha(:,step)=U(j);
        
    %if check==0 && (exist('bb'))
    %Pc(bb)=0;
    %end
    %the implicit scheme will be coded here, it will have U as input
    %and Pc_implicit as output
    
    
    if mod(step,recordstep)==0
        T(:,step/recordstep)=U(i);
        alpha(:,step/recordstep)=U(j);
        time(step/recordstep)=step*dt;
        if (mean(alpha(nodesnearignition,step/recordstep))>curelimit && check==0)
            A=A0; A11=A011; Im1=Im01; B=B0; B11=B011;   %these conversions not needed for ignition by heating
            %F=F0; 
            check=1;
            disp('ignition stopped at time step');
            disp(step)
            ignoffstp=step;
        end
        
        if(mean(alpha(nodesnearignition,step/recordstep))>curelimit && mean(alpha(nodesnearresist,step/recordstep))>curelimit)
        F=F0;
        elseif(mean(alpha(nodesnearignition,step/recordstep))>curelimit && mean(alpha(nodesnearresist,step/recordstep))<curelimit)
            F=Fresistive;
        elseif(mean(alpha(nodesnearignition,step/recordstep))<curelimit && mean(alpha(nodesnearresist,step/recordstep))>curelimit)
            F=Fignition;
        end
    end
    
    if mod(step,recordstep/10)==0        
        step
    end
   toc 
end



%for only the thermal problem,
%equation to be solved is CT*T(n)=DT*T(n-1)+F, where DT=CT-K
%the above loop generates a matrix with each column containing the nodal
%values at certain time step.

%Deleting variables not needed any more

%clear check i j
%clear K CT Cc DT 
%clear K0 CT0 Cc0 DT0 
%clear A B AB A0 B0 AB0

%clear A A0 A011 A11 A12 A21 A22 B B0 B11 B12 B21 B22 Im1 Im01 Im1 Im21 Im22

Axisdistanceset=1e3*(0:deltaz:hz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%visualising Temperature distribution
%The following command rearranges the temperature vector to a matrix for 
%clearer visualisation
%m in Tm and alpham stands for matrix
% Tmss=zeros(m+1,n+1);
% for i=1:te/recordstep
% Tm(:,:,i)=zeros(m+1,n+1);
% alpham(:,:,i)=zeros(m+1,n+1);
% end

%changing temperature and curing vector to matrix
% for i=1:te/recordstep
% for ii=1:m+1
%     jj=(ii-1)*(n+1)+1:(ii-1)*(n+1)+1+n;
%     Tm(ii,:,i)=T(jj,i);
%     %Tmss(ii,:)=Tss(jj);
%     alpham(ii,:,i)=alpha(jj,i);
% end
% end
%temperature and curing vector changed to a matrix

%clear T alpha
%clear U j jj
%clear bl br Fss

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%postprocessing in 'postprocessor' file







