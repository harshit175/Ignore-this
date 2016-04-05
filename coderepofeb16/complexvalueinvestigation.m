Tc=U(1:totalnodes);
alphac=U(totalnodes+1:2*totalnodes);

for ii=1:m+1
    jj=(ii-1)*(n+1)+1:(ii-1)*(n+1)+1+n;
    Tmc(ii,:)=Tc(jj);
    %Tmss(ii,:)=Tss(jj);
    alphamc(ii,:)=alphac(jj);
end