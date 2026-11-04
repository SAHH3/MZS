function GPV = GPV(xy, xyp)

N = size(xy,1)/2;
x = xy(1:N); y = xy(N+1:2*N);
xp = xyp(1:N); yp = xyp(N+1:2*N);

GPV = (log((x-xp).^2 + (y-yp).^2) - log(a(xy, xyp))) ./ (4*pi);