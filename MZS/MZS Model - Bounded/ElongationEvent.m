function [value, isterminal, direction] = ElongationEvent(t, x)
%ELONGATIONEVENT  ODE event: stop when any vortex reaches aspect ratio 3.
%
%   For use with odeset('Events', @ElongationEvent). With X + iY =
%   (sqrt(lambda) - 1/sqrt(lambda)) exp(2 i phi), |X|^2 = lambda + 1/lambda - 2,
%   so |X|^2 = 4/3 is aspect ratio lambda = 3 (the elongation beyond which
%   an isolated elliptical vortex is unstable - Love's criterion for instability).
%
%   t  time (unused)
%   x  state vector [x; y; X; Y], each block of length N

N = length(x)/4;
value = max(x(2*N+1:3*N).^2 + x(3*N+1:4*N).^2) - 4/3;   % crosses zero when max |X|^2 = 4/3
isterminal = 1;                                         % stop the integration
direction = 1;                                          % only when |X|^2 is increasing