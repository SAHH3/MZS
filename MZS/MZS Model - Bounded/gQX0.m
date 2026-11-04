function gQX0 = gQX0(xy, XY, parvec)
%GQX0  Coupling between the ellipticity of each vortex and its own image.
%
%   gQX0 = gQX0(xy, XY, parvec)
%
%     gQX0(i) = Re[ (X_i + i Y_i) BB2(Z_i, Z_i) ] / (2 pi)
%
%   The self-interaction counterpart of GQX.m: only the image part BB2
%   contributes. Multiplied by eps_i*sqrt(|X_i|^2+4) it gives the h5 term
%   of energy_bounded_MZS.m.
%
%   xy      [x; y] of length 2N (image-plane positions)
%   XY      [X; Y] of length 2N (vortex shapes)
%   parvec  [kappas; domain_var]; only domain_var is used
%   Returns an N x 1 vector.

N = length(xy)/2;
domain_var = parvec(N+1:length(parvec))';
X = XY(1:N); Y = XY(N+1:2*N);
z = xy(1:N) + 1i*xy(N+1:2*N);

gQX0 = real((X + 1i.*Y) .* BB2(z, z, domain_var)) ./ (2.*pi);