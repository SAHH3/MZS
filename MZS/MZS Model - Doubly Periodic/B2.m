function B2 = B2(x, y, X, Y, varargin)

xp = x;
if nargin == 5
    xp = varargin{1};
end

B2 = ((2.*(X.^2) + Y.^2 + 4) .* F1(x, y, xp) + X.*Y .* F2(x, y, xp)) ./ sqrt(X.^2 + Y.^2 + 4);