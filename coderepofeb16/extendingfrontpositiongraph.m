for count=2:length(Fronttime)
    if Fronttime(count)==Fronttime(count-1) && Fronttime(count)==Fronttime(count-2)
        start=count-3;
        break;
    end
end

slope=(Axisdistanceset(start)-Axisdistanceset(start-5))/(Fronttime(start)-Fronttime(start-5));
Fronttimeextend=Fronttime;
for i=start:length(Fronttime)
    Fronttimeextend(i)=Fronttime(start)+((Axisdistanceset(i)-Axisdistanceset(start)))/slope;

end

plot(Fronttimeextend,Axisdistanceset,'*')