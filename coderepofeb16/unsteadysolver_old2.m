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
n2=20;  %n2 as to be specified if gpratio =1; the code will be erroneous otherwise.
deltar=[(rinter-rin)/n1 (rout-rinter)/n2 (rin-rc)/n3];
end

n=n1+n2+n3;


%global totalnodes;
totalnodes =(m+1)*(n+1);
curingnodes=(m+1)*(n1+1);

T(totalnodes,te)=0; alpha(totalnodes,te)=0;


rinm=[rin rinter rc]; %in stands for inner radius, m stands for material
%So rinm means inner radius of either polymer or PDMS. 
%For instance, innr radius for PDMS (mat=2) is the interface radius i.e. .
%rinter
deltaz=hz/m;

%Assembling KT and CT matrix
K=zeros(totalnodes,totalnodes);
CT=zeros(totalnodes,totalnodes);
Cc=zeros(totalnodes,totalnodes);

%mlump=diag(ctfunction(0,(rinter-rin)/n1,0,deltaz,1,1,1)); 


i=1:4; lumpingcounter=1;
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
    
    k=stiffnessfunction(rleft,rright,zbot,ztop,mat,kc);
    ct=rhocp(mat)*ctfunction(rleft,rright,zbot,ztop,mat,dt);    
    
    a=[(r-1)*(n+1)+c (r-1)*(n+1)+c+1 r*(n+1)+c r*(n+1)+c+1];
   
    K(a(i),a(i))=K(a(i),a(i))+k(i,i);
    CT(a(i),a(i))=CT(a(i),a(i))+ct(i,i);
    
   
    cc=channel*ctfunction(rleft,rright,zbot,ztop,mat,dt);
    Cc(a(i),a(i))=Cc(a(i),a(i))+cc(i,i);
    
    if mat==1 && r==1  %capitalizes on the structured nature of microchannel mesh, so is only computed for the first row of elements.
        mlump(:,lumpingcounter)=diag(ctfunction(rleft,rright,zbot,ztop,mat,1));
        lumpingcounter=lumpingcounter+1;        
    end
    
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
i=1:4;
for c=1:n
    if (c>=n3+1 && c<=(n1+n3))
        mat=1; curingnode=1;
    elseif (c>=n1+n3+1)
        mat=2; curingnode=0;
    elseif (c<=n3)
        mat=3; curingnode=0;
    end
    
for r=1:m
    a=[(r-1)*(n+1)+c (r-1)*(n+1)+c+1 r*(n+1)+c r*(n+1)+c+1];
    alpha(a(i),1)=curingnode*0.000;
end
end
%alpha=alpha';
%initialization complete

F(totalnodes)=0; %load vector initialised at zero, this will be used in the
%ignition stage

F0=F;  %this will be used after hot boundary is insulated


%Assembling load vector (to ignite the monomer)
% i=1:4;
% for c=1:(n-1)
%     if (c<=n1)
%         mat=1;
%     else
%         mat=2;
%     end
%     rleft=rinm(mat)+(c-1)*deltar(mat);
%     rright=rleft+deltar(mat);
%     Qve=heatvectorfunction(rleft,rright,Q);
%     
%     a=[c c+1 c+2 c+3];  %could give an error
%     F(a(i))=F(a(i))+Qve(i);
% end
%Load vector assembled


%Applying constant temperature boundary conditions on right boundary

j=1:(m+1);

bb=1:n1+n3+1;   %bb=ignition boundary nodes, modification maybe needed for conducting material
br=(n+1)*j; %br=right boundary nodes
bl=br-n;    %bl=left boundary nodes


%Tl=000; %messing with this alone will give unphysical results, other changes to 
%left boundary have to be made as well.
%Tr=000;  %Temperatures at left and right boundaries defined

DT=CT-K; %equation to be solved is CT*T(tn)=DT*T(t(n-1))+F;

%Adjusting F values for right boundary
% for p=1:totalnodes
%     F(p)=F(p)+Tr*( sum(DT(p,br))-sum(CT(p,br)) ) + Tl*( sum(DT(p,bl))-sum(CT(p,bl)) );
%     F0(p)=F0(p)+Tr*( sum(DT(p,br))-sum(CT(p,br)) ) + Tl*( sum(DT(p,bl))-sum(CT(p,bl)) );
% end

%Adjusting F values for bottom boundary
for p=1:totalnodes
    F(p)=F(p)+Tb*( sum(DT(p,bb))-sum(CT(p,bb)) );
    %F0(p)=F0(p)+Tb*( sum(DT(p,bb))-sum(CT(p,bb)) );
end

%Setting nodes in F matrix to fixed temperature values corresponding to
%constant temperature boundary conditions
%F(br)=Tr;  
%F0(br)=Tr;
%F(bl)=Tl;
F(bb)=Tb;

Fss=F; %this is the force vector for steady state

%F(br)=0; %this is needed for solving in unsteady state for ignition on
%F0(br)=0; %this is needed for solving in unsteady state for ignition off

F(bb)=0; %this is needed for solving in unsteady state for ignition on
F0(bb)=0; %this is needed for solving in unsteady state for ignition off

%CT=modifymatrices(CT,br); %see the function for explanation
%DT=modifymatrices(DT,br);
%K=modifymatrices(K,br);

%CT=modifymatrices(CT,bl); %see the function for explanation
%DT=modifymatrices(DT,bl);
%K=modifymatrices(K,bl);
CT0=CT; DT0=DT; K0=K; Cc0=Cc; %these are needed after hot boundary is switched off
disp('post ignition matrices created')
CT=modifymatrices(CT,bb); %see the function for explanation
DT=modifymatrices(DT,bb);
K=modifymatrices(K,bb);
%Cc=modifymatrices(Cc,bb);  %not modified as the ignition curing nodes
%remain free, they are not constrained
disp('matrices modified for bottom boundary')
%initialising T vector for t=0;
T(1:totalnodes,1)=Tambient;
%T(br,1)=Tr;
%T(bl,1)=Tl;
T(bb,1)=Tb;

%solving for steady state
%Tss=K\Fss';

%alpha=zeros(totalnodes,1)+.0001;
Pc=zeros(totalnodes,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%solving for unsteady state



% for step=2:te
%     T(:,step)=inv(CT)*DT*T(:,step-1)+inv(CT)*F';
% end
disp('pre ignition matrices inversion started')
A=inv([CT -rho1*Hr*Cc; zeros(totalnodes,totalnodes) Cc]);
B= [DT -rho1*Hr*Cc; zeros(totalnodes,totalnodes) Cc];
AB=A*B;

disp('post ignition matrices inversion started')
A0=inv([CT0 -rho1*Hr*Cc; zeros(totalnodes,totalnodes) Cc]);
B0= [DT0 -rho1*Hr*Cc; zeros(totalnodes,totalnodes) Cc];
AB0=A0*B0;


i=1:totalnodes;
j=totalnodes+1:2*totalnodes;
%ignitionnodes=1:n1;

mu=vertcat(T(:,step-1), alpha(:,step-1));

check=0;
%Pc=zeros(totalnodes,1);
disp('time marching starting')
for step=2:te
    mu = AB*vertcat(T(:,step-1), alpha(:,step-1))+A*vertcat(F', Pc);
    if not(isreal(mu))
        disp('curing values imaginary at time step');
        disp(step)
        complexcurestep=step;
        break
    end
    T(:,step)=mu(i);
    alpha(:,step)=mu(j);
    Pc=Pcfunction(n1,n3,n,m,totalnodes,T(:,step),alpha(:,step),mlump);
    %Pc(bb)=0;
    if (mean(alpha(n+2+n3:n+2+n1+n3,step))>0.95 && check==0)
        AB=AB0; A=A0; F=F0; 
        check=1;
        disp('ignition stopped at time step');
        disp(step)
        ignoffstp=step;
    end
    step
end
clear check i j
%for only the thermal problem,
%equation to be solved is CT*T(n)=DT*T(n-1)+F, where DT=CT-K
%the above loop generates a matrix with each column containing the nodal
%values at certain time step.

%Deleting variables not needed any more
clear K CT Cc DT 
clear K0 CT0 Cc0 DT0 
clear A B AB A0 B0 AB0



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%visualising Temperature distribution
%The following command rearranges the temperature vector to a matrix for 
%clearer visualisation
%m in Tm and alpham stands for matrix
Tmss=zeros(m+1,n+1);
for i=1:te
Tm(:,:,i)=zeros(m+1,n+1);
alpham(:,:,i)=zeros(m+1,n+1);
end

for i=1:te
for ii=1:m+1
    jj=(ii-1)*(n+1)+1:(ii-1)*(n+1)+1+n;
    Tm(ii,:,i)=T(jj,i);
    %Tmss(ii,:)=Tss(jj);
    alpham(ii,:,i)=alpha(jj,i);
end
end
%temperature and curing vector changed to a matrix
clear T alpha
clear mu j jj
clear bl br Fss

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%postprocessing in 'postprocessor' file







