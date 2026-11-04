function D1 = D1(x, y, X, Y, varargin)
%D1  y-derivative of the ellipticity coupling: s d/dy (X F1 + Y F2).
%
%   D1 = D1(x, y, X, Y)       D1 = sqrt(X^2 + Y^2 + 4) (X a(x,y) - Y b(x,y))
%   D1 = D1(x, y, X, Y, xp)   with xp passed on to a.m and b.m
%
%   Gives the contribution of an ellipse's shape (X, Y) to the x-velocity
%   (u_sum in uvUV_doubly_MZS.m), for one row of periodic images at
%   separation (x, y). The optional xp avoids 0/0 on the diagonal (see C.m).

xp = x;
if nargin == 5
    xp = varargin{1};
end

D1 = sqrt(X.^2 + Y.^2 + 4) .* (X .* a(x, y, xp) - Y .* b(x, y, xp));