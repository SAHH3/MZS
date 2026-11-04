function H = energy_bounded_MZS(~, xyvor, XY, eps, parvec, varargin)

% energy_bounded_MZS:
% calculates energy in arbitrary bounded domain using
% conformal map Jacobian in jacob.m
% etavors are the complex coordinates in the circular (image) domain
% based on energy.m by JGE and TLA

N = length(xyvor)/2;
kappas = parvec(1:N); kappaij = kappas*kappas';
domain_var = parvec(N+1:length(parvec))';
z = xyvor(1:N) + 1i*xyvor(N+1:2*N);

ultri = find(repmat(1:N,N,1) ~= repmat((1:N)',1,N));
kappaij = kappaij(ultri);

X = XY(1:N); Y = XY(N+1:2*N);
absX2 = X.^2 + Y.^2;
Xmat = repmat(X',N,1); Ymat = repmat(Y',N,1);
X1 = Xmat(ultri); Y1 = Ymat(ultri);
absX2_ultri = X1.^2 + Y1.^2;

epsmat = repmat(eps',N,1); eps_ultri = epsmat(ultri);

GPV_mat = GPV(xyvor, xyvor'); GQX_mat = GQX(xyvor, xyvor', XY', parvec);

h1 = -sum((kappas.^2) .* log((absX2 + 4)/4)) ./ (8.*pi);
h2 = -sum(kappaij .* GPV_mat(ultri)) ./ 2;
h3 = -sum((kappas.^2) .* g(xyvor)) ./ 2;
h4 = -sum(kappaij .* eps_ultri .* sqrt(absX2_ultri + 4) .* GQX_mat(ultri));
h5 = -sum((kappas.^2) .* eps .* sqrt(absX2 + 4) .* gQX0(xyvor, XY, parvec));
h6 = sum((kappas.^2) .* log(jacob(z, domain_var))) ./ (4.*pi);
h7 = -sum((kappas.^2) .* eps .* sqrt(absX2 + 4) .* real((X + 1i.*Y) .* C1(z, domain_var))) ./ (2.*pi);
H = h1 + h2 + h3 + h4 + h5 + h6 + h7;

cond = 0;
if nargin == 6
    cond = varargin{1};
end
if cond == 1
    H = [h1, h2, h3, h4, h5, h6, h7, H];
end