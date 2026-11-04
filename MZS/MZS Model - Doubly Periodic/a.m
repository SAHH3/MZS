function a = a(x, y, varargin)

xp = x;
if nargin == 3
    xp = varargin{1};
end

a = sin(y) .* (cosh(x) .* (cosh(x) + cos(y)) - 2) ./ ((cosh(xp) - cos(y)).^3);