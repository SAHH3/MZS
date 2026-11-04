function D1 = D1(x, y, X, Y, varargin)

xp = x;
if nargin == 5
    xp = varargin{1};
end

D1 = sqrt(X.^2 + Y.^2 + 4) .* (X .* a(x, y, xp) - Y .* b(x, y, xp));