clear K CT Cc DT 
clear K0 CT0 Cc0 DT0 
clear A B AB A0 B0 AB0



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%visualising Temperature distribution
%The following command rearranges the temperature vector to a matrix for 
%clearer visualisation
%m in Tm and alpham stands for matrix
%Tmss=zeros(m+1,n+1);
for i=1:1222000/recordstep
Tm(:,:,i)=zeros(m+1,n+1);
alpham(:,:,i)=zeros(m+1,n+1);
end

for i=1:1222000/recordstep
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