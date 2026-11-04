function[alpha] = alpha_heart(q, g, A)
%ALPHA_HEART  Scale factor a(q,g) of the 'heart' conformal map.
%
%   alpha = alpha_heart(q, g, A)
%
%   The heart map is z = F(Z) = a Z / ((1 - q^2 Z^2)(1 - i g Z)), taking the
%   unit circle |Z| < 1 to the physical domain. alpha is the value of a that
%   makes the physical domain have area A. Throughout this code A = pi, so
%   the heart has the same area as the unit circle.
%
%   q, g  heart shape parameters (q = g = 0 gives the unit circle, alpha = 1)
%   A     required area of the physical domain

alpha = sqrt(A/pi)*(1-g^2)*(1+q^2*g^2)*(1-q^4)                   ...
 /sqrt(1+q^4+2*g^6*q^8-g^4*q^2*(q^6+4*q^4+q^2-4)+g^2*(4*q^6-2*q^4-4*q^2));