function F2 = F2(x, y, varargin)

xp = x;
if nargin == 3
    xp = varargin{1};
end

F2 = -sinh(x) .* sin(y) ./ ((cosh(xp) - cos(y)).^2);