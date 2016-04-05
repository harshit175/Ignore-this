function [ qhe ] = heatsourcevectorfunction( rleft,rright,zbot,ztop,Q)
%HEATSOURCEVECTORFUNCTION Summary of this function goes here
%   Detailed explanation goes here
qhe=Q*[ ((rleft - rright)*(2*rleft + rright)*(zbot - ztop))/12, 
    ((rleft - rright)*(rleft + 2*rright)*(zbot - ztop))/12, 
    ((rleft - rright)*(rleft + 2*rright)*(zbot - ztop))/12, 
    ((rleft - rright)*(2*rleft + rright)*(zbot - ztop))/12];


end

