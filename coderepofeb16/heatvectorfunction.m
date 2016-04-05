function [ Qve ] = heatvectorfunction( rleft,rright,Q )
%HEATVECTORFUNCTION Summary of this function goes here
%   Detailed explanation goes here


Qve=Q*[ -((rleft - rright)*(2*rleft + rright))/6, -((rleft - rright)*(rleft + 2*rright))/6, 0, 0];
end

