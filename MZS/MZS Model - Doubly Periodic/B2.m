function B2 = B2(x, y, X, Y, varargin)
%B2  X-derivative of the ellipticity coupling: d/dX [ s (X F1 + Y F2) ].
%
%   B2 = B2(x, y, X, Y)       B2 = ((2X^2 + Y^2 + 4) F1 + X Y F2) / s
%   B2 = B2(x, y, X, Y, xp)   with xp passed on to F1.m and F2.m
%
%   where s = sqrt(X^2 + Y^2 + 4). Gives the rate of change of Y due to the
%   strain from a point vortex at separation (x, y) (V_sum in
%   uvUV_doubly_MZS.m). The optional xp avoids 0/0 on the diagonal (see C.m).

xp = x;
if nargin == 5
    xp = varargin{1};
end

B2 = ((2.*(X.^2) + Y.^2 + 4) .* F1(x, y, xp) + X.*Y .* F2(x, y, xp)) ./ sqrt(X.^2 + Y.^2 + 4);