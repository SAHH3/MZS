function GPV = GPV(xy, NI)
%GPV  Point-vortex Green's function of the doubly periodic domain [0, 2pi)^2.
%
%   GPV = GPV(xy, NI)
%
%     G(x, y) = [ sum_{n=-NI}^{NI} log( (cosh(x - 2 pi n) - cos y) / cosh(2 pi n) )
%                 - x^2 / (2 pi) - g00 ] / (4 pi)
%
%   Each log term is a row of images along y (period 2 pi); the sum over n
%   adds the rows at x = 2 pi n (Weiss & McWilliams summation), dividing by
%   cosh(2 pi n) so the sum converges. The -x^2/(2 pi) term makes G periodic
%   in x and accounts for the uniform background vorticity that a periodic
%   domain needs. g00 is a constant fixing the additive constant of G.
%   The terms decay like exp(-2 pi n), so NI = 6 is ample.
%
%   xy  [x; y] of length 2N (positions in [0, 2pi))
%   NI  number of image rows on each side
%   Returns G(x_i - x_j, y_i - y_j) for all i ~= j, as a column vector in
%   the order given by the off-diagonal index list ultri.

g00 = 0.354043395939778;
N = size(xy,1)/2;
% linear indices of all off-diagonal (i ~= j) entries of an N x N matrix
ultri = find(repmat(1:N,N,1) ~= repmat((1:N)',1,N));

% separations x_i - x_j, y_i - y_j for every pair
xdists = repmat(xy(1:N),1,N) - repmat(xy(1:N),1,N)';
ydists = repmat(xy(N+1:2*N),1,N) - repmat(xy(N+1:2*N),1,N)';
xdists = xdists(ultri); ydists = ydists(ultri);

GPV = 0;
for II=NI:-1:1          % image rows at x = -/+ 2 pi II

    GPV = GPV + log((cosh(xdists - 2*pi*II) - cos(ydists)) ./ cosh(2*pi*II));
    GPV = GPV + log((cosh(xdists + 2*pi*II) - cos(ydists)) ./ cosh(2*pi*II));

end

% central row, background term and constant
GPV = GPV + log(cosh(xdists) - cos(ydists)) - xdists.^2./(2 * pi) - g00;
GPV = GPV ./ (4*pi);