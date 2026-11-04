function [ uvUV ] = uvUV_bounded_MZS( ~, xyXYvor, eps, parvec)
%
% uvUV_bounded_MZS:
% returns values of (u,v) from vortex positions
% xyvor in an arbitrary bounded domain
% needs conformal map routines DF1, DF2
% nb xyvor must be row vector size=(2*N,1)

N = length(xyXYvor)/4; xy = xyXYvor(1:2*N); XY = xyXYvor(2*N+1:4*N);
kappas = parvec(1:N); domain_vars = parvec(N+1:length(parvec))';

x = xy(1:N); y = xy(N+1:2*N);
X = XY(1:N); Y = XY(N+1:2*N); modX2 = X.^2 + Y.^2;

z = x + 1i*y; zp = transpose(z); zbar = conj(z); zpbar = conj(zp);
modz2 = abs(z).^2; df1 = DF1(z, domain_vars); df2 = DF2(z, domain_vars);
df2df1 = df2./df1; J2 = abs(df1).^2; aij = a_complex(z, zp);
C = C1(z, domain_vars); dC = dC1(z, domain_vars); C2 = BB2(z, z, domain_vars);
B1 = BB1(zp, z, domain_vars); B2 = BB2(zp, z, domain_vars);
c1 = b44(z, z, domain_vars); c2 = b33(z, z, domain_vars);
b1 = b11(z, zp, domain_vars); b2 = b22(z, zp, domain_vars);
b3 = b33(z, zp, domain_vars); b4 = b44(z, zp, domain_vars);

u = 0*xyXYvor(1:N); v = 0*xyXYvor(1:N); UV1 = 0*xyXYvor(1:N); UV2 = 0*xyXYvor(1:N);

% u & v eqns
uiv1 = zbar ./ (1 - modz2) - df2df1 ./ 2;
XYi = eps .* (X + 1i.*Y) .* sqrt(modX2 + 4); XYj = transpose(XYi);
u = u + kappas .* imag(uiv1 + XYi .* (dC + c1 - c2));
v = v + kappas .* real(uiv1 + XYi .* (dC + c1 + c2));

uv_mat = 1./(z - zp) + zpbar .* (1 - zbar .* zp) ./ aij + XYi .* (b2 + b4);
u_mat = uv_mat + XYj .* (b1 - b3); u_mat(1:N+1:end) = zeros(1,N);
v_mat = uv_mat + XYj .* (b1 + b3); v_mat(1:N+1:end) = zeros(1,N);
u = u + imag(u_mat) * kappas; v = v + real(v_mat) * kappas;

u = u ./ (2.*pi .* J2); v = v ./ (2.*pi .* J2);

% U & V eqns
XYU = X.^2 + 2.*(Y.^2) + 4 - 1i.*X.*Y;
XYV = 2.*(X.^2) + Y.^2 + 4 + 1i.*X.*Y;

UV1 = UV1 + kappas .* (X - 1i.*Y) ./ (2.*eps .* sqrt(modX2 + 4));
UV_mat = B1 + B2; UV_mat(1:N+1:end) = zeros(1,N);
UV2 = UV2 + kappas .* (C + C2) + UV_mat * kappas;

U = imag(UV1 + XYU .* UV2); V = real(UV1 + XYV .* UV2);

U = U ./ (2.*pi .* sqrt(modX2 + 4)); V = V ./ (2.*pi .* sqrt(modX2 + 4));

uvUV = vertcat(u, v, U, V);