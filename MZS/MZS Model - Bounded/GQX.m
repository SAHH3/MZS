function GQX = GQX(xy, xyp, XY, parvec)
%GQX  Coupling between the ellipticity of vortex j and a point vortex i.
%
%   GQX = GQX(xy, xyp, XY, parvec)
%
%     GQX(i,j) = Re[ (X_j + i Y_j) ( BB1(Z_i, Z_j) + BB2(Z_i, Z_j) ) ] / (2 pi)
%
%   BB1 + BB2 is the second derivative of the Green's function with respect
%   to the physical position of vortex j, i.e. the strain at j due to a
%   point vortex at i. Multiplied by eps_j*sqrt(|X_j|^2+4) it gives the
%   h4 term of energy_bounded_MZS.m.
%
%   xy      [x; y] of length 2N as a COLUMN (image-plane positions Z_i)
%   xyp     the same positions as a ROW (Z_j)
%   XY      [X; Y] of length 2N as a ROW (shape of vortex j)
%   parvec  [kappas; domain_var]; only domain_var is used
%   Returns the N x N matrix GQX(i,j); the diagonal is not used.

N = length(xy)/2;
domain_var = parvec(N+1:length(parvec))';
X = XY(1:N); Y = XY(N+1:2*N);
z = xy(1:N) + 1i*xy(N+1:2*N); zp = xyp(1:N) + 1i*xyp(N+1:2*N);

GQX = real((X + 1i.*Y) .* (BB1(z, zp, domain_var) + BB2(z, zp, domain_var))) ./ (2.*pi);