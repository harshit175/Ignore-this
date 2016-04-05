if gpratio~=1
nint=int8(log(1-(rout-rinter)/(rinter-rin)*(1-gpratio)*n1)/log(gpratio));
%nint=nint+1;
n2=double(nint);
for i=1:n2+1
    r2(i)=rinter+ (rinter-rin)/n1*(1-gpratio^(i-1))/(1-gpratio); %i=1, r2=rinter; i=2, r2=rinter+a;  i=3, r2=rinter+a+a*gpratio
    deltar=[(rinter-rin)/n1 (rout-rinter)/n2];
end
else
n2=20;  %n2 as to be specified if gpratio =1; the code will be erroneous otherwise.
deltar=[(rinter-rin)/n1 (rout-rinter)/n2];
end

n=n1+n2;


%global totalnodes;
totalnodes =(m+1)*(n+1);
curingnodes=(m+1)*(n1+1);

rinm=[rin rinter]; %in stands for inner radius, m stands for material
deltaz=hz/m;

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
    
    
    
    a=[(r-1)*(n+1)+c (r-1)*(n+1)+c+1 r*(n+1)+c r*(n+1)+c+1]; %These are node numbers corresponding to r and c
    rcoor=[rleft rright rright rleft rleft];
    zcoor=[zbot zbot ztop ztop zbot];
    
    
    
    if mat==1 
        plot(rcoor,zcoor,'r')
    else
        plot(rcoor,zcoor,'k')
    end
    %
    hold on
    
end
end
axis equal
axis off