function [ uvUV ] = uvUV_doubly_MZS( ~ , xyXYvor, eps, parvec)
%UVUV_DOUBLY_MZS  Right-hand side of the MZS equations in the doubly periodic domain.
%
%   uvUV = uvUV_doubly_MZS(t, xyXYvor, eps, parvec)
%
%   Returns d/dt of the state [x; y; X; Y] (see energy_doubly_MZS.m for the
%   variables). These are Hamilton's equations for the energy H in
%   energy_doubly_MZS.m:
%
%     k_i dx_i/dt =  dH/dy_i,          k_i dy_i/dt = -dH/dx_i,
%     k_i eps_i dX_i/dt =  dH/dY_i,    k_i eps_i dY_i/dt = -dH/dX_i.
%
%   t        time (unused; the equations are autonomous)
%   xyXYvor  state [x; y; X; Y], 4N x 1 column vector
%   eps      N x 1 vortex sizes
%   parvec   [NI; kappas]: number of image rows (see GPV.m), then N circulations

N = size(xyXYvor,1)/4;
NI = parvec(1); kappas = parvec(2:N+1);
xyvor = xyXYvor(1:2*N); XY = xyXYvor(2*N+1:4*N);

% N x N separations x_i - x_j, y_i - y_j, wrapped to the nearest image in [-pi, pi)
xdists = mod(repmat(xyvor(1:N),1,N) - repmat(xyvor(1:N),1,N)' + pi*ones(N,N), 2*pi) - pi;
ydists = mod(repmat(xyvor(N+1:2*N),1,N) - repmat(xyvor(N+1:2*N),1,N)' + pi*ones(N,N), 2*pi) - pi;
% the same with the diagonal (i = j) shifted, to avoid 0/0 there; the
% diagonal entries are removed below
xdistsp = xdists + diag(ones(N,1)); ydistsp = ydists + diag(ones(N,1));

X = XY(1:N); Y = XY(N+1:2*N);

u_sum = 0*xyXYvor(1:N); v_sum = 0*xyXYvor(1:N); U_sum = 0*xyXYvor(1:N); V_sum = 0*xyXYvor(1:N);

% Sum over image rows. Each *_sum becomes an N x N matrix: row i is the
% vortex being moved, column j the vortex doing the moving.
%   C:  point-vortex velocity. For v, the rows of images are taken along y
%       instead (C(ydists, xdists)): by the x <-> y symmetry of the square
%       domain this is the same G, and it avoids the -x^2/(2pi) background term.
%   D1, D2: velocity from the ellipticity of i (eps, X, Y) and of j (eps', X', Y').
%   B1, B2: strain from point vortex j acting on the shape of vortex i.
for II=NI:-1:1
    u_sum = u_sum + C(xdists + 2.*II.*pi, ydists) + eps .* D1(xdists + 2.*II.*pi, ydists, X, Y) + eps' .* D1(xdists + 2.*II.*pi, ydists, X', Y');
    u_sum = u_sum + C(xdists - 2.*II.*pi, ydists) + eps .* D1(xdists - 2.*II.*pi, ydists, X, Y) + eps' .* D1(xdists - 2.*II.*pi, ydists, X', Y');

    v_sum = v_sum + C(ydists + 2.*II.*pi, xdists) + eps .* D2(xdists + 2.*II.*pi, ydists, X, Y) + eps' .* D2(xdists + 2.*II.*pi, ydists, X', Y');
    v_sum = v_sum + C(ydists - 2.*II.*pi, xdists) + eps .* D2(xdists - 2.*II.*pi, ydists, X, Y) + eps' .* D2(xdists - 2.*II.*pi, ydists, X', Y');

    U_sum = U_sum - B1(xdists + 2.*II.*pi, ydists, X, Y) - B1(xdists - 2.*II.*pi, ydists, X, Y);

    V_sum = V_sum - B2(xdists + 2.*II.*pi, ydists, X, Y) - B2(xdists - 2.*II.*pi, ydists, X, Y);
end

% central row of images, then remove the i = j entries
u_sum = u_sum + C(xdists, ydists, xdistsp) + eps .* D1(xdists, ydists, X, Y, xdistsp) + eps' .* D1(xdists, ydists, X', Y', xdistsp);
u_sum = u_sum - diag(diag(u_sum));

v_sum = v_sum + C(ydists, xdists, ydistsp) + eps .* D2(xdists, ydists, X, Y, xdistsp) + eps' .* D2(xdists, ydists, X', Y', xdistsp);
v_sum = v_sum - diag(diag(v_sum));

% for U and V, also the strain from the background term (-X/(8 pi^2) in GQX.m)
U_sum = U_sum - B1(xdists, ydists, X, Y, xdistsp);
U_sum = U_sum + X.*Y ./ (2.*pi .* sqrt(X.^2 + Y.^2 + 4));
U_sum = U_sum - diag(diag(U_sum));

V_sum = V_sum - B2(xdists, ydists, X, Y, xdistsp);
V_sum = V_sum + (2.*(X.^2) + Y.^2 + 4) ./ (2.*pi .* sqrt(X.^2 + Y.^2 + 4));
V_sum = V_sum - diag(diag(V_sum));

u = 0*xyXYvor(1:N); v = 0*xyXYvor(1:N); U = 0*xyXYvor(1:N); V = 0*xyXYvor(1:N);

% sum over j, weighted by k_j
u = u - u_sum * kappas ./ (4.*pi);
v = v + v_sum * kappas ./ (4.*pi);

% first term: rotation of each ellipse by its own vorticity (from h1);
% Y./pi and X./pi: the h4 (background) term, summed over all j including j = i
U = U - (kappas .* Y ./ (eps .* (X.^2 + Y.^2 + 4)) - (U_sum + Y./pi) * kappas) ./ (4.*pi);
V = V + (kappas .* X ./ (eps .* (X.^2 + Y.^2 + 4)) - (V_sum + X./pi) * kappas) ./ (4.*pi);

uvUV = vertcat(u, v, U, V);