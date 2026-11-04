function F2 = F2(x, y, varargin)
%F2  Imaginary part (with a minus sign) of the strain from one row of periodic images.
%
%   F2 = F2(x, y)       F2 = -sinh(x) sin(y) / (cosh(x) - cos(y))^2
%   F2 = F2(x, y, xp)   the same, but with cosh(xp) in the denominator
%
%   F1 - i F2 = 1 / (1 - cosh(x + i y)); see F1.m.
%
%   The optional xp avoids 0/0 on the diagonal (see C.m).

xp = x;
if nargin == 3
    xp = varargin{1};
end

F2 = -sinh(x) .* sin(y) ./ ((cosh(xp) - cos(y)).^2);