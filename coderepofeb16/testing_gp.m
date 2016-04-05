clear all

rin=0e-3; %inner radius
rinter=1e-3; %interface radius
rout=15e-3;   %outer radius

Tambient=50;
Tb=200;
%conductivity, (c suffix added to k, to distinguish from k)
kc=[.19 .27]; %arbitrary values added for now

rho1=1.19e3;
cp1=3118;

rho2=1.03e3;
cp2=1100;

Hr=550e3;
%density and enthalpy, lumped together as rhocp
rhocp=[rho1*cp1 rho2*cp2]; %arbitrary values added for now

Q=0; %heat supplied to start the process

hz=05e-3; %height in z direction;

%elements along radial direction
n1=10;

gpratio=1.05;

if gpratio~=1
nint=int8(log(1-(rout-rinter)/(rinter-rin)*(1-gpratio)*n1)/log(gpratio));
nint=nint+1;
n2=double(nint);
for i=1:n2+1
    r2(i)=rinter+ (rinter-rin)/n1*(1-gpratio^(i-1))/(1-gpratio); %i=1, r2=rinter; i=2, r2=rinter+a;  i=3, r2=rinter+a+a*gpratio
    deltar=[(rinter-rin)/n1 (rout-rinter)/n2];
end
else
n2=70;
deltar=[(rinter-rin)/n1 (rout-rinter)/n2];
end

n=n1+n2;

%elements along z direction
m=40;
%global totalnodes;
totalnodes =(m+1)*(n+1);
curingnodes=(m+1)*(n1+1);
%time step in seconds;
%dt=2;
dt=.01;
te=3;
T(totalnodes,te)=0; alpha(totalnodes,te)=0;


rinm=[rin rinter]; %in stands for inner radius, m stands for material
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
    if (c<=n1)
        mat=1; channel=1;
    else
        mat=2; channel=0;
    end
for r=1:m
    if mat==1
        rleft=rinm(mat)+(c-1)*deltar(mat);
        rright=rleft+deltar(mat);
    else
        if gpratio~=1
            rleft=r2(c-n1);
            rright=r2(c-n1+1);
        else
            rleft=rinm(mat)+(c-n1-1)*deltar(mat);
            rright=rleft+deltar(mat);
        end
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
    
    if mat==1 && r==1
        mlump(:,lumpingcounter)=diag(ctfunction(rleft,rright,zbot,ztop,mat,1));
        lumpingcounter=lumpingcounter+1;
        
    end
    
end
end