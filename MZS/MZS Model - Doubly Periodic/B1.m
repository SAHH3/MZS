function B1 = B1(x, y, X, Y, varargin)
%B1  Y-derivative of the ellipticity coupling: d/dY [ s (X F1 + Y F2) ].
%
%   B1 = B1(x, y, X, Y)       B1 = (X Y F1 + (X^2 + 2Y^2 + 4) F2) / s
%   B1 = B1(x, y, X, Y, xp)   with xp passed on to F1.m and F2.m
%
%   where s = sqrt(X^2 + Y^2 + 4). Gives the rate of change of X due to the
%   strain from a point vortex at separation (x, y) (U_sum in
%   uvUV_doubly_MZS.m). The optional xp avoids 0/0 on the diagonal (see C.m).

xp = x;
if nargin == 5
    xp = varargin{1};
end

B1 = (X.*Y .* F1(x, y, xp) + (X.^2 + 2.*(Y.^2) + 4) .* F2(x, y, xp)) ./ sqrt(X.^2 + Y.^2 + 4);