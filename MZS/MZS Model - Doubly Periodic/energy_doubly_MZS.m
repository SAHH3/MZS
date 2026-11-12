function [ h ] = energy_doubly_MZS( ~, xyvor, XY, eps, parvec, varargin)
%ENERGY_DOUBLY_MZS  Hamiltonian (energy) of N elliptical vortices (MZS model)
%                   in the doubly periodic domain [0, 2pi)^2.
%
%   H     = energy_doubly_MZS(t, xyvor, XY, eps, parvec)
%   Hlist = energy_doubly_MZS(t, xyvor, XY, eps, parvec, 1)
%
%   Each vortex is an ellipse of aspect ratio lambda_i and orientation phi_i
%   centred at (x_i, y_i), stored as
%     X_i + i Y_i = (sqrt(lambda_i) - 1/sqrt(lambda_i)) exp(2 i phi_i),
%   with semi-axes sqrt(8 eps_i lambda_i), sqrt(8 eps_i / lambda_i) (area 8 pi eps_i).
%
%     H = h1 + h2 + h3 + h4 + h5, where
%     h1 = -sum_i k_i^2 log((|X_i|^2 + 4)/4) / (8 pi)                       self-energy of each ellipse
%     h2 = -1/2 sum_{i~=j} k_i k_j G(x_ij, y_ij)                            point-vortex interactions (GPV.m)
%     h3 = -sum_{i~=j} k_i k_j eps_i s_i GQX(i,j)                           ellipticity of i with point vortex j (GQX.m)
%     h4 = (sum_i k_i) sum_j k_j eps_j (|X_j|^2 + 2) / (8 pi^2)             second moment of each ellipse with the uniform background vorticity; zero for a neutral gas
%     h5 = -1/2 sum_i k_i^2 g                                               interaction of each vortex with its own periodic images
%   with s_i = sqrt(|X_i|^2 + 4), and g = lim_{r->0} [ G(x,y) - log(r)/(2 pi) ]
%   (r^2 = x^2 + y^2) the regular part of the Green's function of GPV.m at
%   zero separation; the counterpart of h3 in energy_bounded_MZS.m. h5 is a
%   constant, so it does not affect the motion.
%
%   t        time (unused; present so the call matches the ODE right-hand side)
%   xyvor    [x; y], 2N x 1, positions
%   XY       [X; Y], 2N x 1, vortex shapes
%   eps      N x 1 vortex sizes
%   parvec   [NI; kappas]: number of image rows (see GPV.m), then N circulations
%   (6th argument) 1 to return [h1, h2 + h5, h3, h4, H] instead of just H
%                  (the point-vortex terms h2 and h5 are returned together)

N = size(xyvor,1)/2;                     % number of vortices
NI = parvec(1); kappas = parvec(2:N+1);  % extract params

% linear indices of all off-diagonal (i ~= j) entries of an N x N matrix
ultri = find(repmat(1:N,N,1) ~= repmat((1:N)',1,N));
kappaij = kappas*kappas'; kappaij = kappaij(ultri);  % note scaling

X = XY(1:N); Y = XY(N+1:2*N);
absX2 = X.^2 + Y.^2;
% for each pair (i,j), the shape and size of vortex i (matrix row index)
Xmat = repmat(X,1,N); Ymat = repmat(Y,1,N);
X = Xmat(ultri); Y = Ymat(ultri);
absX2_ultri = X.^2 + Y.^2;
g = -0.083929;                           % regular part of G at zero separation (see header)

epsmat = repmat(eps,1,N); eps_ultri = epsmat(ultri);

h1 = -sum((kappas.^2) .* log((absX2+4)/4)) ./ (8*pi);
h2 = -sum(kappaij .* GPV(xyvor, NI)) ./ 2;
h3 = -sum(kappaij .* eps_ultri .* sqrt(absX2_ultri + 4) .* GQX(xyvor, XY, NI, 1));
h4 = -sum(-(eps' .* (absX2' + 2)) * (kappas * kappas') ./ (8.*(pi^2)));
h5 = -sum(kappas.^2) .* g ./ 2;
h = h1 + h2 + h3 + h4 + h5;

% optionally return the individual terms as well
cond = 0;
if nargin == 6
    cond = varargin{1};
end
if cond == 1
    h = [h1, h2 + h5, h3, h4, h];
end