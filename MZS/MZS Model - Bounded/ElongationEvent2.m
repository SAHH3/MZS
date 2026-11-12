function [value, isterminal, direction] = ElongationEvent2(t, x)
%ELONGATIONEVENT2  ODE event: stop when any vortex reaches |X|^2 = 25.
%
%   As ElongationEvent.m, but with a much larger threshold: |X|^2 =
%   lambda + 1/lambda - 2 = 25 is aspect ratio lambda of about 27. Used by
%   the gap runs, where a single vortex is strongly stretched passing
%   through the gap. This lets the simulation run for longer, to follow the
%   vortex's evolution further: it stops before the integration breaks down,
%   rather than when the vortex becomes unstable. Change 25 as you prefer.
%
%   t  time (unused)
%   x  state vector [x; y; X; Y], each block of length N

N = length(x)/4;
value = max(x(2*N+1:3*N).^2 + x(3*N+1:4*N).^2) - 25;    % crosses zero when max |X|^2 = 25
isterminal = 1;                                         % stop the integration
direction = 1;                                          % only when |X|^2 is increasing