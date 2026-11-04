function F1 = F1(x, y, varargin)
%F1  Real part of the strain from one row of periodic images.
%
%   F1 = F1(x, y)       F1 = (1 - cosh(x) cos(y)) / (cosh(x) - cos(y))^2
%   F1 = F1(x, y, xp)   the same, but with cosh(xp) in the denominator
%
%   F1 - i F2 = 1 / (1 - cosh(x + i y)), the second complex derivative of
%   2 log sinh((x + i y)/2), where Re[2 log sinh((x+iy)/2)] = log(cosh x - cos y)
%   + const is one row of the periodic Green's function (GPV.m). So
%       X F1 + Y F2 = Re[ (X + i Y) / (1 - cosh(x + i y)) ]
%   is the coupling of an ellipse with shape (X, Y) to a point vortex at
%   separation (x, y). Used in GQX.m, B1.m and B2.m.
%
%   The optional xp avoids 0/0 on the diagonal (see C.m).

xp = x;
if nargin == 3
    xp = varargin{1};
end

F1 = (1 - cosh(x) .* cos(y)) ./ ((cosh(xp) - cos(y)).^2);