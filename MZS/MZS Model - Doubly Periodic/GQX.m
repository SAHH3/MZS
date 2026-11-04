function GQX = GQX(xy, XY, NI, cond)

N = size(xy,1)/2;
X = XY(1:N); Y = XY(N+1:2*N);

if cond == 1

    ultri = find(repmat(1:N,N,1) ~= repmat((1:N)',1,N));

    xdists = repmat(xy(1:N),1,N) - repmat(xy(1:N),1,N)';
    ydists = repmat(xy(N+1:2*N),1,N) - repmat(xy(N+1:2*N),1,N)';
    xdists = xdists(ultri); ydists = ydists(ultri);
    xdistsp = xdists; ydistsp = ydists;

    Xmat = repmat(X,1,N); Ymat = repmat(Y,1,N);
    X = Xmat(ultri); Y = Ymat(ultri);

elseif cond == 2

    xdists = mod(repmat(xy(1:N),1,N) - repmat(xy(1:N),1,N)' + pi*ones(N,N), 2*pi) - pi;
    ydists = mod(repmat(xy(N+1:2*N),1,N) - repmat(xy(N+1:2*N),1,N)' + pi*ones(N,N), 2*pi) - pi;
    xdistsp = xdists + diag(ones(N,1)); ydistsp = ydists + diag(ones(N,1));

end

GQX = 0;
for II=NI:-1:1

    GQX = GQX + X .* F1(xdists + 2*pi*II, ydists) + Y .* F2(xdists + 2*pi*II, ydists);
    GQX = GQX + X .* F1(xdists - 2*pi*II, ydists) + Y .* F2(xdists - 2*pi*II, ydists);

end

GQX = GQX + X .* F1(xdists, ydists, xdistsp) + Y .* F2(xdists, ydists, xdistsp);
GQX = GQX ./ (4*pi);
GQX = GQX - X ./ (8 .* (pi.^2));

if cond == 2
    GQX = GQX - diag(diag(GQX));
end