function a = a(xy, xyp)

N = length(xy)/2;
x = xy(1:N); y = xy(N+1:2*N);
xp = xyp(1:N); yp = xyp(N+1:2*N);

a = 1 - 2*(x.*xp + y.*yp) + (x.^2 + y.^2) .* (xp.^2 + yp.^2);