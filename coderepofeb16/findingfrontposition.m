X=horzcat(0:deltar(1):rinter,r2(1:length(r2)-1));
Y=0:deltaz:hz;
Y=10e2*Y;
X=10e2*X;

tp=[1500 3000 4500];
sizealpha=size(alpham);
for plots=1:3
k=1;
for i=1:sizealpha(2)-1
    for j=1:sizealpha(1)
        %js=[j-1 j j+1];
        %is=[i-1 i i+1];
        if alpham(j,i,tp(plots))<.999
            fp(k,1)=j;
            fp(k,2)=i;
            k=k+1;
            break;
        end
    end
end
plot(X(fp(:,2)),Y(fp(:,1)));
hold on
end

axis('equal',[0 10e2*rinter 0 10e2*hz])
ylabel('length (mm)','FontSize', 12)
xlabel('radius (mm)','FontSize', 12)