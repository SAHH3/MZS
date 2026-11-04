function b = b(x, y, varargin)

xp = x;
if nargin == 3
    xp = varargin{1};
end

b = sinh(x) .* (cos(y) .* (cosh(x) + cos(y)) - 2) ./ ((cosh(xp) - cos(y)).^3);