function [ uvUV ] = uvUV_bounded_MZS( ~, xyXYvor, eps, parvec)
%UVUV_BOUNDED_MZS  Right-hand side of the MZS equations in a bounded domain.
%
%   uvUV = uvUV_bounded_MZS(t, xyXYvor, eps, parvec)
%
%   Returns d/dt of the state [x; y; X; Y] (see energy_bounded_MZS.m for the
%   variables). These are Hamilton's equations for the energy H in
%   energy_bounded_MZS.m:
%
%     k_i J_i^2 dx_i/dt =  dH/dy_i,    k_i J_i^2 dy_i/dt = -dH/dx_i,
%     k_i eps_i dX_i/dt =  dH/dY_i,    k_i eps_i dY_i/dt = -dH/dX_i,
%
%   where J_i = |F'(Z_i)| and (x_i, y_i) are image-plane (unit-disc) positions.
%
%   t        time (unused; the equations are autonomous)
%   xyXYvor  state [x; y; X; Y], 4N x 1 column vector
%   eps      N x 1 vortex sizes
%   parvec   [kappas; domain_var]: N circulations, then the map parameters
%   Needs DF1-DF4 from one of the Conformal folders on the path.

N = length(xyXYvor)/4; xy = xyXYvor(1:2*N); XY = xyXYvor(2*N+1:4*N);
kappas = parvec(1:N); domain_vars = parvec(N+1:length(parvec))';

x = xy(1:N); y = xy(N+1:2*N);
X = XY(1:N); Y = XY(N+1:2*N); modX2 = X.^2 + Y.^2;

% z is a column (vortex i), zp the same positions as a row (vortex j), so
% two-argument functions below return N x N matrices over all pairs (i,j)
z = x + 1i*y; zp = transpose(z); zbar = conj(z); zpbar = conj(zp);
modz2 = abs(z).^2; df1 = DF1(z, domain_vars); df2 = DF2(z, domain_vars);
df2df1 = df2./df1; J2 = abs(df1).^2; aij = a_complex(z, zp);
% strain terms for the (X,Y) equations: self (C, C2) and from vortex j (B1, B2)
C = C1(z, domain_vars); dC = dC1(z, domain_vars); C2 = BB2(z, z, domain_vars);
B1 = BB1(zp, z, domain_vars); B2 = BB2(zp, z, domain_vars);
% derivatives of the strain terms, for the (x,y) equations: self (c1, c2),
% from the ellipticity of j (b1, b3) and from the ellipticity of i (b2, b4)
c1 = b44(z, z, domain_vars); c2 = b33(z, z, domain_vars);
b1 = b11(z, zp, domain_vars); b2 = b22(z, zp, domain_vars);
b3 = b33(z, zp, domain_vars); b4 = b44(z, zp, domain_vars);

u = 0*xyXYvor(1:N); v = 0*xyXYvor(1:N); UV1 = 0*xyXYvor(1:N); UV2 = 0*xyXYvor(1:N);

% u & v eqns
% self-induced motion: own image (from h3) and conformal-map term (from h6)
uiv1 = zbar ./ (1 - modz2) - df2df1 ./ 2;
% XYi = eps_i s_i (X_i + i Y_i), s_i = sqrt(|X_i|^2 + 4) (column, vortex i);
% XYj the same as a row (vortex j)
XYi = eps .* (X + 1i.*Y) .* sqrt(modX2 + 4); XYj = transpose(XYi);
u = u + kappas .* imag(uiv1 + XYi .* (dC + c1 - c2));
v = v + kappas .* real(uiv1 + XYi .* (dC + c1 + c2));

% pair terms: point vortex j and its image, plus the ellipticity terms
uv_mat = 1./(z - zp) + zpbar .* (1 - zbar .* zp) ./ aij + XYi .* (b2 + b4);
u_mat = uv_mat + XYj .* (b1 - b3); u_mat(1:N+1:end) = zeros(1,N);   % no i = j terms
v_mat = uv_mat + XYj .* (b1 + b3); v_mat(1:N+1:end) = zeros(1,N);
u = u + imag(u_mat) * kappas; v = v + real(v_mat) * kappas;          % sum over j, weighted by k_j

% convert from dH/dZ to image-plane velocity
u = u ./ (2.*pi .* J2); v = v ./ (2.*pi .* J2);

% U & V eqns
% shape-dependent weights from differentiating s_i (X_i + i Y_i) with respect to Y_i and X_i
XYU = X.^2 + 2.*(Y.^2) + 4 - 1i.*X.*Y;
XYV = 2.*(X.^2) + Y.^2 + 4 + 1i.*X.*Y;

% UV1: rotation of each ellipse by its own vorticity (from h1)
UV1 = UV1 + kappas .* (X - 1i.*Y) ./ (2.*eps .* sqrt(modX2 + 4));
% UV2: total strain on vortex i: self (own image and map) plus from every point vortex j
UV_mat = B1 + B2; UV_mat(1:N+1:end) = zeros(1,N);
UV2 = UV2 + kappas .* (C + C2) + UV_mat * kappas;

U = imag(UV1 + XYU .* UV2); V = real(UV1 + XYV .* UV2);

U = U ./ (2.*pi .* sqrt(modX2 + 4)); V = V ./ (2.*pi .* sqrt(modX2 + 4));

uvUV = vertcat(u, v, U, V);