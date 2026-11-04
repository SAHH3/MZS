function b = b(x, y, varargin)
%B  Derivative of the strain functions: b = dF1/dx = -dF2/dy.
%
%   b = b(x, y)       b = sinh(x) (cos(y) (cosh(x) + cos(y)) - 2) / (cosh(x) - cos(y))^3
%   b = b(x, y, xp)   the same, but with cosh(xp) in the denominator
%
%   With a.m, gives the velocity induced by an ellipse's shape (D1.m, D2.m).
%   The optional xp avoids 0/0 on the diagonal (see C.m).

xp = x;
if nargin == 3
    xp = varargin{1};
end

b = sinh(x) .* (cos(y) .* (cosh(x) + cos(y)) - 2) ./ ((cosh(xp) - cos(y)).^3);