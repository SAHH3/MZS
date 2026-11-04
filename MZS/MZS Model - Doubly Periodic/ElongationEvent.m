function [value, isterminal, direction] = ElongationEvent(t,x)
N = length(x)/4;
value = max(x(2*N+1:3*N).^2 + x(3*N+1:4*N).^2) - 4/3;
isterminal = 1;
direction = 1;