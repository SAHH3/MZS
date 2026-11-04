function GPV = GPV(xy, NI)

g00 = 0.354043395939778;
N = size(xy,1)/2;
ultri = find(repmat(1:N,N,1) ~= repmat((1:N)',1,N));

xdists = repmat(xy(1:N),1,N) - repmat(xy(1:N),1,N)';
ydists = repmat(xy(N+1:2*N),1,N) - repmat(xy(N+1:2*N),1,N)';
xdists = xdists(ultri); ydists = ydists(ultri);

GPV = 0;
for II=NI:-1:1

    GPV = GPV + log((cosh(xdists - 2*pi*II) - cos(ydists)) ./ cosh(2*pi*II));
    GPV = GPV + log((cosh(xdists + 2*pi*II) - cos(ydists)) ./ cosh(2*pi*II));

end

GPV = GPV + log(cosh(xdists) - cos(ydists)) - xdists.^2./(2 * pi) - g00;
GPV = GPV ./ (4*pi);