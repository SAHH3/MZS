function GPV = GPV(xy, xyp)
%GPV  Point-vortex Green's function of the unit disc.
%
%   GPV = GPV(xy, xyp)
%
%     G(Z, Z') = ( log|Z - Z'|^2 - log|1 - Z conj(Z')|^2 ) / (4 pi)
%              = log| (Z - Z') / (1 - Z conj(Z')) | / (2 pi)
%
%   i.e. a vortex at Z' plus its image at 1/conj(Z'), so that the unit circle
%   is a streamline. Because the map is conformal, G(Z(z), Z(z')) is also
%   the Green's function of the physical domain. Gives the h2 term of
%   energy_bounded_MZS.m.
%
%   xy   [x; y] of length 2N as a COLUMN (image-plane positions Z_i)
%   xyp  [x'; y'] of length 2N as a ROW (positions Z_j)
%   Returns the N x N matrix G(Z_i, Z_j); the diagonal is -Inf and is not used.

N = size(xy,1)/2;
x = xy(1:N); y = xy(N+1:2*N);
xp = xyp(1:N); yp = xyp(N+1:2*N);

GPV = (log((x-xp).^2 + (y-yp).^2) - log(a(xy, xyp))) ./ (4*pi);