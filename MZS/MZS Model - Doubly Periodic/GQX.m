function GQX = GQX(xy, XY, NI, cond)
%GQX  Coupling between the ellipticity of vortex i and a point vortex j,
%     in the doubly periodic domain.
%
%   GQX = GQX(xy, XY, NI, cond)
%
%     GQX(i,j) = sum_{n=-NI}^{NI} [ X_i F1(x_ij - 2 pi n, y_ij) + Y_i F2(x_ij - 2 pi n, y_ij) ] / (4 pi)
%                - X_i / (8 pi^2),          x_ij = x_i - x_j,  y_ij = y_i - y_j
%
%   i.e. Re[(X_i + i Y_i) * (second complex derivative of G)] summed over the
%   image rows of GPV.m; -X_i/(8 pi^2) comes from the -x^2/(2 pi) background
%   term of G. Multiplied by eps_i*sqrt(|X_i|^2+4) it gives the h3 term of
%   energy_doubly_MZS.m.
%
%   xy    [x; y] of length 2N (positions)
%   XY    [X; Y] of length 2N (vortex shapes)
%   NI    number of image rows on each side
%   cond  1: return a column over all pairs i ~= j (order of the index list
%            ultri), using unwrapped separations. This is what energy_doubly_MZS.m uses.
%         2: return the N x N matrix, with separations wrapped to [-pi, pi)
%            and the diagonal set to zero (not currently used).

N = size(xy,1)/2;
X = XY(1:N); Y = XY(N+1:2*N);

if cond == 1

    % linear indices of all off-diagonal (i ~= j) entries of an N x N matrix
    ultri = find(repmat(1:N,N,1) ~= repmat((1:N)',1,N));

    xdists = repmat(xy(1:N),1,N) - repmat(xy(1:N),1,N)';
    ydists = repmat(xy(N+1:2*N),1,N) - repmat(xy(N+1:2*N),1,N)';
    xdists = xdists(ultri); ydists = ydists(ultri);
    xdistsp = xdists; ydistsp = ydists;

    % shape of vortex i (matrix row index) for each pair
    Xmat = repmat(X,1,N); Ymat = repmat(Y,1,N);
    X = Xmat(ultri); Y = Ymat(ultri);

elseif cond == 2

    xdists = mod(repmat(xy(1:N),1,N) - repmat(xy(1:N),1,N)' + pi*ones(N,N), 2*pi) - pi;
    ydists = mod(repmat(xy(N+1:2*N),1,N) - repmat(xy(N+1:2*N),1,N)' + pi*ones(N,N), 2*pi) - pi;
    % shift the diagonal (i = j) to avoid 0/0; it is removed at the end
    xdistsp = xdists + diag(ones(N,1)); ydistsp = ydists + diag(ones(N,1));

end

GQX = 0;
for II=NI:-1:1          % image rows at x = -/+ 2 pi II

    GQX = GQX + X .* F1(xdists + 2*pi*II, ydists) + Y .* F2(xdists + 2*pi*II, ydists);
    GQX = GQX + X .* F1(xdists - 2*pi*II, ydists) + Y .* F2(xdists - 2*pi*II, ydists);

end

% central row, then the background term
GQX = GQX + X .* F1(xdists, ydists, xdistsp) + Y .* F2(xdists, ydists, xdistsp);
GQX = GQX ./ (4*pi);
GQX = GQX - X ./ (8 .* (pi.^2));

if cond == 2
    GQX = GQX - diag(diag(GQX));
end