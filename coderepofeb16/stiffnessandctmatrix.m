%This code will calculate the value of kt and ct which will be used
%in the stiffnessfunction and ctfunction.

%rleft, rright, zbot, ztop is input, stiffness matrix is output.

%rleft=radius value at left node
%rright=radius value at right node

%zbot=z value at bottom node
%ztop=z value at top node

%xi is along x axis
%eta is along y axis

syms rleft rright zbot ztop xi eta real ;

Deltar=rright-rleft;
Deltaz=ztop-zbot;

J=[2/Deltar 0;
    0 2/Deltaz]; %this is actually the inverse

j=det(J); %this is actually the inverse of determinant of jacobian

B=J*[-(1-eta)/4  (1-eta)/4 (1+eta)/4 -(1+eta)/4;
     -(1-xi)/4  -(1+xi)/4 (1+xi)/4 (1-xi)/4];  %this is correct
 
r=0.5*[rleft*(1-xi)+rright*(1+xi)];

integrandK=r*j*B'*B;

f=1/sqrt(3);

K=subs(integrandK, [xi, eta],[f, f])+subs(integrandK, [xi, eta],[-f, -f])+subs(integrandK, [xi, eta],[f, -f])+subs(integrandK, [xi, eta],[-f, f]);

%N=0.25*[(1-xi)*(1-eta) (1+xi)*(1-eta) (1+xi)*(1+eta) (1-xi)*(1+eta)];

%integrandCT=r*j*N'*N;

%ct=subs(integrandCT, [xi, eta],[f, f])+subs(integrandCT, [xi, eta],[-f, -f])+subs(integrandCT, [xi, eta],[f, -f])+subs(integrandCT, [xi, eta],[-f, f]);




