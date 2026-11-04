function DF3 = DF3(Z, domain_var)
%DF3  Third derivative F'''(Z) of the 'heart' conformal map.
%
%   The map is z = F(Z) = a Z / ((1 - q^2 Z^2)(1 - i g Z)), a = alpha_heart(q,g,pi).
%   Derivative of the partial-fraction form in DF2.m:
%
%     F'''(Z) = 3a [ -2 g^4 / ((q^2 + g^2)(i + g Z)^4)
%                    + q^3 ( 1/((q - i g)(q Z - 1)^4) + 1/((q + i g)(q Z + 1)^4) ) ]
%
%   Z           image-plane position(s), any array size (evaluated elementwise)
%   domain_var  [q, g], the heart shape parameters

q = domain_var(1); g = domain_var(2); a = alpha_heart(q,g,pi);

DF3 = 3.*a .* (-(2 .* (g.^4) ./ ((q.^2 + g.^2) .* (1i + g.*Z).^4)) ...
               + (q.^3) .* ((1./((q - 1i.*g) .* (q.*Z - 1).^4)) + (1./((q + 1i.*g) .* (q.*Z + 1).^4))));

% unit disc (F(Z) = Z): the formula above is 0/0, but F''' = 0
if (q == 0) && (g == 0)
    DF3 = Z.*0;
end