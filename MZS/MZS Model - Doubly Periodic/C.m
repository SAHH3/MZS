function C = C(x, y, varargin)

xp = x;
if nargin == 3
    xp = varargin{1};
end

C = sin(y) ./ (cosh(xp) - cos(y));