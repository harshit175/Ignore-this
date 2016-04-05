clear all
load('36x8.mat')
postprocesser

plot(time,dTnode(9,:))
hold on

xlabel('Time (sec.)')
ylabel ('Temperature (C)')
title('Temperature variation with time at a point')
%title('Temperature derivative variation at a point')
%legend('80 nodes')

clear all

load('72x16.mat')
postprocesser

plot(time,dTnode(17,:),'red')

clear all

% load('160x8x05.mat')
% postprocesser
% 
% plot(time,dTnode(81,:),'green')
% clear all

% load('160x8x05.mat')
% postprocesser
% 
% plot(time,dTnode(81,:),'k')
% clear all

%legend('160 nodes')