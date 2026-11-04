function GQX = GQX(xy, xyp, XY, parvec)

N = length(xy)/2;
domain_var = parvec(N+1:length(parvec))';
X = XY(1:N); Y = XY(N+1:2*N);
z = xy(1:N) + 1i*xy(N+1:2*N); zp = xyp(1:N) + 1i*xyp(N+1:2*N);

GQX = real((X + 1i.*Y) .* (BB1(z, zp, domain_var) + BB2(z, zp, domain_var))) ./ (2.*pi);