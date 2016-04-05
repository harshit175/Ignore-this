lumpingcounter=0;
global multfactor
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