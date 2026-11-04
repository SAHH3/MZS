function a = a(xy, xyp)
%A  The image-vortex factor |1 - Z conj(Z')|^2 for the unit circle, real form.
%
%   a = a(xy, xyp)
%
%     a = |1 - Z conj(Z')|^2 = 1 - 2 (x x' + y y') + |Z|^2 |Z'|^2,
%     Z = x + i y,  Z' = x' + i y'
%
%   Appears in the unit-circle Green's function (GPV.m) and its Robin
%   function (g.m). Same quantity as a_complex.m, but from real coordinates.
%
%   xy   [x; y] of length 2N (image-plane positions)
%   xyp  [x'; y'] of length 2N. Passing xy as a column and xyp as a row
%        gives the N x N matrix a(i,j) for every pair (implicit expansion).

N = length(xy)/2;
x = xy(1:N); y = xy(N+1:2*N);
xp = xyp(1:N); yp = xyp(N+1:2*N);

a = 1 - 2*(x.*xp + y.*yp) + (x.^2 + y.^2) .* (xp.^2 + yp.^2);