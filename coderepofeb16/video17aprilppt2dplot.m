clear FF
te=double(te);
iii=1:5:te/recordstep;
%iii=1:10:788;
loops = length(iii);
FF(loops) = struct('cdata',[],'colormap',[]);

%load('rad3point6te6000')
%X=horzcat(0:deltar(3):rin-rc-deltar(3),rin-rc:deltar(1):rinter-deltar(1),r2);
X=horzcat(0:deltar(3):rin-rc-deltar(3),rin-rc:deltar(1):rinter);
%X=horzcat(0:deltar(3):rin-rc);
Y=0:deltaz:hz;


j=1;
for i=iii
%for i=1:10:250

figure(1)
%contourf(10e2*X(1:46),10e2*Y,alpham(:,1:46,i))
plot(10e2*Y,Tm(:,1,i))
%title('Degree of cure/Temperature for r_{ch}= 0.17 mm','FontSize', 12)
%title('Degree of cure/Temperature for resistive heating','FontSize', 12)
%title('Degree of cure/Temperature','FontSize', 12)

%title('Temperature (C)','FontSize', 12)
%title('Degree of Cure','FontSize', 12)
%surf(Tm(:,:,i))
xlabel('z (mm)','FontSize', 12)
ylabel('T (C)','FontSize', 12)
%ylabel('Cure','FontSize', 12)
%axis([-inf inf 0 7])
%axis([-inf inf -inf inf 25 250])

%axis([-inf inf 20 220])
axis([-inf inf 0 1.1])
%axis equal
%colormap gray
%colorbar
%caxis([25 280])

%zlim manual;
%zlabel('Temperature (C)')
%zlabel('Degree of cure')
%drawnow
FF(j)=getframe(figure(1));
j=j+1;
%figure(2)
%surf(alpham(:,:,i))
%pause(0.01)
end