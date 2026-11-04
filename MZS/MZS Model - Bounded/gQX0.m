function gQX0 = gQX0(xy, XY, parvec)

N = length(xy)/2;
domain_var = parvec(N+1:length(parvec))';
X = XY(1:N); Y = XY(N+1:2*N);
z = xy(1:N) + 1i*xy(N+1:2*N);

gQX0 = real((X + 1i.*Y) .* BB2(z, z, domain_var)) ./ (2.*pi);