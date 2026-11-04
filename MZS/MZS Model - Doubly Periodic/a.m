function a = a(x, y, varargin)
%A  Derivative of the strain functions: a = dF1/dy = dF2/dx.
%
%   a = a(x, y)       a = sin(y) (cosh(x) (cosh(x) + cos(y)) - 2) / (cosh(x) - cos(y))^3
%   a = a(x, y, xp)   the same, but with cosh(xp) in the denominator
%
%   With b.m, gives the velocity induced by an ellipse's shape (D1.m, D2.m).
%   The optional xp avoids 0/0 on the diagonal (see C.m).

xp = x;
if nargin == 3
    xp = varargin{1};
end

a = sin(y) .* (cosh(x) .* (cosh(x) + cos(y)) - 2) ./ ((cosh(xp) - cos(y)).^3);