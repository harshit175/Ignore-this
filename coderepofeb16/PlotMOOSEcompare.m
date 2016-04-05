plot(1e3*(0:deltaz:hz),Tm(:,1,1000),'r')
hold on
plot(1e3*MOOSEdist2,MOOSETemp2,'-.')

xlabel('x (mm)')

ylabel('T (C)')

legend('Result from MATLAB','Result from MOOSE')

title('Comparison at t=0.04 s')