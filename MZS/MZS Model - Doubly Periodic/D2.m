function D2 = D2(x, y, X, Y, varargin)
%D2  x-derivative of the ellipticity coupling: s d/dx (X F1 + Y F2).
%
%   D2 = D2(x, y, X, Y)       D2 = sqrt(X^2 + Y^2 + 4) (X b(x,y) + Y a(x,y))
%   D2 = D2(x, y, X, Y, xp)   with xp passed on to a.m and b.m
%
%   Gives the contribution of an ellipse's shape (X, Y) to the y-velocity
%   (v_sum in uvUV_doubly_MZS.m), for one row of periodic images at
%   separation (x, y). The optional xp avoids 0/0 on the diagonal (see C.m).

xp = x;
if nargin == 5
    xp = varargin{1};
end

D2 = sqrt(X.^2 + Y.^2 + 4) .* (X .* b(x, y, xp) + Y .* a(x, y, xp));