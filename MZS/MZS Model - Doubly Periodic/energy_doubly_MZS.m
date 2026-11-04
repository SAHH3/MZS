function [ h ] = energy_doubly_MZS( ~, xyvor, XY, eps, parvec, varargin)

% function to return value of hamiltonian from vortex positions (x,y),
% (X,Y) and circulations in the doubly periodic domain.

N = size(xyvor,1)/2;                     % number of vortices
NI = parvec(1); kappas = parvec(2:N+1);  % extract params

ultri = find(repmat(1:N,N,1) ~= repmat((1:N)',1,N));
kappaij = kappas*kappas'; kappaij = kappaij(ultri);  % note scaling

X = XY(1:N); Y = XY(N+1:2*N);
absX2 = X.^2 + Y.^2;
Xmat = repmat(X,1,N); Ymat = repmat(Y,1,N);
X = Xmat(ultri); Y = Ymat(ultri);
absX2_ultri = X.^2 + Y.^2;

epsmat = repmat(eps,1,N); eps_ultri = epsmat(ultri);

h1 = -sum((kappas.^2) .* log((absX2+4)/4)) ./ (8*pi);
h2 = -sum(kappaij .* GPV(xyvor, NI)) ./ 2;
h3 = -sum(kappaij .* eps_ultri .* sqrt(absX2_ultri + 4) .* GQX(xyvor, XY, NI, 1));
h4 = -sum(-(eps' .* (absX2' + 2)) * (kappas * kappas') ./ (8.*(pi^2)));
h = h1 + h2 + h3 + h4;

cond = 0;
if nargin == 6
    cond = varargin{1};
end
if cond == 1
    h = [h1, h2, h3, h4, h];
end