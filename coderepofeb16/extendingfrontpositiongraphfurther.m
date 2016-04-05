axisinterval=Axisdistanceset(length(Axisdistanceset))-Axisdistanceset(length(Axisdistanceset)-1);
Lengthinmm=100;
while(Axisdistanceset(length(Axisdistanceset))<=Lengthinmm)
    Axisdistanceset(length(Axisdistanceset)+1)=Axisdistanceset(length(Axisdistanceset))+axisinterval;
end

slope=axisinterval/(Fronttimeextend(length(Fronttimeextend))-Fronttimeextend(length(Fronttimeextend)-1));

start=length(Fronttimeextend);

for i=start:length(Axisdistanceset)
    Fronttimeextend(i)=Fronttimeextend(start)+((Axisdistanceset(i)-Axisdistanceset(start)))/slope;

end

plot(Fronttimeextend,Axisdistanceset,'*')