%vectortomatrixconvert
plot(Axisdistanceset,Tm(:,1,1/recordstep))

xlabel('Distance (mm)','fontsize',12)
ylabel('Temperature','FontSize',12)
hold on
load('analytic_axialtemp_dist')

%plot(x*1e3,Tanalytical,'-.','color','r')

%legend({strcat('FE soln., {\Delta}t =',num2str(dt),' s'),'analytical soln.'},'fontsize',12)


title('Comparison of Temperature at t=0.2 s.','fontsize',12)
title('Changed initial condition','fontsize',12)