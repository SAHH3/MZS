function B1 = B1(x, y, X, Y, varargin)

xp = x;
if nargin == 5
    xp = varargin{1};
end

B1 = (X.*Y .* F1(x, y, xp) + (X.^2 + 2.*(Y.^2) + 4) .* F2(x, y, xp)) ./ sqrt(X.^2 + Y.^2 + 4);