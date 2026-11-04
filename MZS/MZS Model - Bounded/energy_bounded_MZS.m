function H = energy_bounded_MZS(~, xyvor, XY, eps, parvec, varargin)
%ENERGY_BOUNDED_MZS  Hamiltonian (energy) of N elliptical vortices (MZS model)
%                    in a bounded domain given by a conformal map of the unit circle.
%
%   H     = energy_bounded_MZS(t, xyvor, XY, eps, parvec)
%   Hlist = energy_bounded_MZS(t, xyvor, XY, eps, parvec, 1)
%
%   Positions are image-plane (unit-circle) positions Z_i = x_i + i y_i; the
%   physical positions are F(Z_i), with F the conformal map (heart_map.m or
%   gap_map.m, whichever Conformal folder is on the path). Each vortex is an
%   ellipse of aspect ratio lambda_i and orientation phi_i, stored as
%     X_i + i Y_i = (sqrt(lambda_i) - 1/sqrt(lambda_i)) exp(2 i phi_i),
%   with semi-axes sqrt(8 eps_i lambda_i), sqrt(8 eps_i / lambda_i) (area 8 pi eps_i).
%   Below, s_i = sqrt(|X_i|^2 + 4).
%
%     H = h1 + h2 + ... + h7, where
%     h1 = -sum_i k_i^2 log((|X_i|^2 + 4)/4) / (8 pi)       self-energy of each ellipse
%     h2 = -1/2 sum_{i~=j} k_i k_j G(Z_i, Z_j)              point-vortex interactions (GPV.m)
%     h3 = -1/2 sum_i k_i^2 g(Z_i)                          interaction with own image (g.m)
%     h4 = -sum_{i~=j} k_i k_j eps_j s_j GQX(i,j)           ellipticity of j with point vortex i (GQX.m)
%     h5 = -sum_i k_i^2 eps_i s_i gQX0(i)                   ellipticity of i with own image (gQX0.m)
%     h6 =  sum_i k_i^2 log|F'(Z_i)| / (4 pi)               conformal-map correction to h3 (jacob.m)
%     h7 = -sum_i k_i^2 eps_i s_i Re((X_i + i Y_i) C1(Z_i)) / (2 pi)
%                                                           conformal-map correction to h5 (C1.m)
%
%   The main scripts group these as self-interaction h1, point vortex
%   h2+h3+h6 and elliptic h4+h5+h7.
%
%   t        time (unused; present so the call matches the ODE right-hand side)
%   xyvor    [x; y], 2N x 1, image-plane positions
%   XY       [X; Y], 2N x 1, vortex shapes
%   eps      N x 1 vortex sizes
%   parvec   [kappas; domain_var]: N circulations, then the map parameters
%   (6th argument) 1 to return [h1, ..., h7, H] instead of just H
%
%   Based on energy.m by JGE and TLA.

N = length(xyvor)/2;
kappas = parvec(1:N); kappaij = kappas*kappas';
domain_var = parvec(N+1:length(parvec))';
z = xyvor(1:N) + 1i*xyvor(N+1:2*N);             % complex image-plane positions Z_i

% linear indices of all off-diagonal (i ~= j) entries of an N x N matrix
ultri = find(repmat(1:N,N,1) ~= repmat((1:N)',1,N));
kappaij = kappaij(ultri);                        % k_i k_j for each pair

X = XY(1:N); Y = XY(N+1:2*N);
absX2 = X.^2 + Y.^2;
% for each pair (i,j), the shape and size of vortex j (matrix column index)
Xmat = repmat(X',N,1); Ymat = repmat(Y',N,1);
X1 = Xmat(ultri); Y1 = Ymat(ultri);
absX2_ultri = X1.^2 + Y1.^2;

epsmat = repmat(eps',N,1); eps_ultri = epsmat(ultri);

% N x N matrices over all pairs: rows i (column inputs), columns j (row inputs)
GPV_mat = GPV(xyvor, xyvor'); GQX_mat = GQX(xyvor, xyvor', XY', parvec);

h1 = -sum((kappas.^2) .* log((absX2 + 4)/4)) ./ (8.*pi);
h2 = -sum(kappaij .* GPV_mat(ultri)) ./ 2;
h3 = -sum((kappas.^2) .* g(xyvor)) ./ 2;
h4 = -sum(kappaij .* eps_ultri .* sqrt(absX2_ultri + 4) .* GQX_mat(ultri));
h5 = -sum((kappas.^2) .* eps .* sqrt(absX2 + 4) .* gQX0(xyvor, XY, parvec));
h6 = sum((kappas.^2) .* log(jacob(z, domain_var))) ./ (4.*pi);
h7 = -sum((kappas.^2) .* eps .* sqrt(absX2 + 4) .* real((X + 1i.*Y) .* C1(z, domain_var))) ./ (2.*pi);
H = h1 + h2 + h3 + h4 + h5 + h6 + h7;

% optionally return the individual terms as well
cond = 0;
if nargin == 6
    cond = varargin{1};
end
if cond == 1
    H = [h1, h2, h3, h4, h5, h6, h7, H];
end