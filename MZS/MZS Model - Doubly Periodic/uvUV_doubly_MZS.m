function [ uvUV ] = uvUV_doubly_MZS( ~ , xyXYvor, eps, parvec)

% function to return values of (u,v) and (U,V), derivatives of (x,y) and (X,Y)
% from vortex positions (x,y), (X,Y) and circulations in the doubly periodic domain.

N = size(xyXYvor,1)/4;
NI = parvec(1); kappas = parvec(2:N+1);
xyvor = xyXYvor(1:2*N); XY = xyXYvor(2*N+1:4*N);

xdists = mod(repmat(xyvor(1:N),1,N) - repmat(xyvor(1:N),1,N)' + pi*ones(N,N), 2*pi) - pi;
ydists = mod(repmat(xyvor(N+1:2*N),1,N) - repmat(xyvor(N+1:2*N),1,N)' + pi*ones(N,N), 2*pi) - pi;
xdistsp = xdists + diag(ones(N,1)); ydistsp = ydists + diag(ones(N,1));

X = XY(1:N); Y = XY(N+1:2*N);

u_sum = 0*xyXYvor(1:N); v_sum = 0*xyXYvor(1:N); U_sum = 0*xyXYvor(1:N); V_sum = 0*xyXYvor(1:N);

for II=NI:-1:1
    u_sum = u_sum + C(xdists + 2.*II.*pi, ydists) + eps .* D1(xdists + 2.*II.*pi, ydists, X, Y) + eps' .* D1(xdists + 2.*II.*pi, ydists, X', Y');
    u_sum = u_sum + C(xdists - 2.*II.*pi, ydists) + eps .* D1(xdists - 2.*II.*pi, ydists, X, Y) + eps' .* D1(xdists - 2.*II.*pi, ydists, X', Y');

    v_sum = v_sum + C(ydists + 2.*II.*pi, xdists) + eps .* D2(xdists + 2.*II.*pi, ydists, X, Y) + eps' .* D2(xdists + 2.*II.*pi, ydists, X', Y');
    v_sum = v_sum + C(ydists - 2.*II.*pi, xdists) + eps .* D2(xdists - 2.*II.*pi, ydists, X, Y) + eps' .* D2(xdists - 2.*II.*pi, ydists, X', Y');

    U_sum = U_sum - B1(xdists + 2.*II.*pi, ydists, X, Y) - B1(xdists - 2.*II.*pi, ydists, X, Y);

    V_sum = V_sum - B2(xdists + 2.*II.*pi, ydists, X, Y) - B2(xdists - 2.*II.*pi, ydists, X, Y);
end

u_sum = u_sum + C(xdists, ydists, xdistsp) + eps .* D1(xdists, ydists, X, Y, xdistsp) + eps' .* D1(xdists, ydists, X', Y', xdistsp);
u_sum = u_sum - diag(diag(u_sum));

v_sum = v_sum + C(ydists, xdists, ydistsp) + eps .* D2(xdists, ydists, X, Y, xdistsp) + eps' .* D2(xdists, ydists, X', Y', xdistsp);
v_sum = v_sum - diag(diag(v_sum));

U_sum = U_sum - B1(xdists, ydists, X, Y, xdistsp);
U_sum = U_sum + X.*Y ./ (2.*pi .* sqrt(X.^2 + Y.^2 + 4));
U_sum = U_sum - diag(diag(U_sum));

V_sum = V_sum - B2(xdists, ydists, X, Y, xdistsp);
V_sum = V_sum + (2.*(X.^2) + Y.^2 + 4) ./ (2.*pi .* sqrt(X.^2 + Y.^2 + 4));
V_sum = V_sum - diag(diag(V_sum));

u = 0*xyXYvor(1:N); v = 0*xyXYvor(1:N); U = 0*xyXYvor(1:N); V = 0*xyXYvor(1:N);

u = u - u_sum * kappas ./ (4.*pi);
v = v + v_sum * kappas ./ (4.*pi);

U = U - (kappas .* Y ./ (eps .* (X.^2 + Y.^2 + 4)) - (U_sum + Y./pi) * kappas) ./ (4.*pi);
V = V + (kappas .* X ./ (eps .* (X.^2 + Y.^2 + 4)) - (V_sum + X./pi) * kappas) ./ (4.*pi);

uvUV = vertcat(u, v, U, V);