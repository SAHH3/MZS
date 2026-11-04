function F1 = F1(x, y, varargin)

xp = x;
if nargin == 3
    xp = varargin{1};
end

F1 = (1 - cosh(x) .* cos(y)) ./ ((cosh(xp) - cos(y)).^2);