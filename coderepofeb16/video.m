clear FF
iii=1:70:te;

loops = length(i);
FF(loops) = struct('cdata',[],'colormap',[]);

j=1;
for i=iii
%for i=1:10:250

figure(1)
surf(alpham(:,:,i))
%surf(Tm(:,:,i))
xlabel('nodes along r-direction')
ylabel('nodes along z-direction')
%zlim manual;
%zlabel('Temperature (C)')
zlabel('Degree of cure')
%drawnow
FF(j)=getframe(figure(1));
j=j+1;
%figure(2)
%surf(alpham(:,:,i))
%pause(0.01)
end