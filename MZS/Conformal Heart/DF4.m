function DF4 = DF4(Z, domain_var)
%DF4  Fourth derivative F''''(Z) of the 'heart' conformal map.
%
%   The map is z = F(Z) = a Z / ((1 - q^2 Z^2)(1 - i g Z)), a = alpha_heart(q,g,pi).
%   Derivative of the partial-fraction form in DF3.m:
%
%     F''''(Z) = 12a [ 2 g^5 / ((q^2 + g^2)(i + g Z)^5)
%                      - q^4 ( 1/((q - i g)(q Z - 1)^5) + 1/((q + i g)(q Z + 1)^5) ) ]
%
%   Z           image-plane position(s), any array size (evaluated elementwise)
%   domain_var  [q, g], the heart shape parameters

q = domain_var(1); g = domain_var(2); a = alpha_heart(q,g,pi);

DF4 = 12.*a .* ((2 .* (g.^5) ./ ((q.^2 + g.^2) .* (1i + g.*Z).^5)) ...
                - (q.^4) .* ((1./((q - 1i.*g) .* (q.*Z - 1).^5)) + (1./((q + 1i.*g) .* (q.*Z + 1).^5))));

% unit disc (F(Z) = Z): the formula above is 0/0, but F'''' = 0
if (q == 0) && (g == 0)
    DF4 = Z.*0;
end