function D2 = D2(x, y, X, Y, varargin)

xp = x;
if nargin == 5
    xp = varargin{1};
end

D2 = sqrt(X.^2 + Y.^2 + 4) .* (X .* b(x, y, xp) + Y .* a(x, y, xp));