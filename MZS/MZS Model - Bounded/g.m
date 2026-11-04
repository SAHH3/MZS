function g = g(xy)
%G  Robin function of the unit disc (regular part of the Green's function at Z' = Z).
%
%   g = g(xy)
%
%     g(Z) = -log((1 - |Z|^2)^2) / (4 pi) = -log(1 - |Z|^2) / (2 pi)
%
%   This is the interaction of a point vortex at Z with its own image in the
%   unit circle; it gives the h3 term of energy_bounded_MZS.m.
%
%   xy  [x; y] of length 2N (image-plane positions); returns an N x 1 vector

g = -log(a(xy, xy)) ./ (4*pi);