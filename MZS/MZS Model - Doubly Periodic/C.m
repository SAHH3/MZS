function C = C(x, y, varargin)
%C  Derivative of one row of the periodic Green's function: d/dy log(cosh x - cos y).
%
%   C = C(x, y)       C = sin(y) / (cosh(x) - cos(y))
%   C = C(x, y, xp)   the same, but with cosh(xp) in the denominator
%
%   log(cosh x - cos y) is the stream function of a vortex and its images
%   along the y direction (period 2 pi); see GPV.m. C gives the velocity
%   it induces (u_sum and v_sum in uvUV_doubly_MZS.m).
%
%   The optional xp is used on the diagonal (i = j, where x = y = 0) to avoid
%   0/0: uvUV_doubly_MZS.m passes xp = x + 1 there and discards the result.

xp = x;
if nargin == 3
    xp = varargin{1};
end

C = sin(y) ./ (cosh(xp) - cos(y));