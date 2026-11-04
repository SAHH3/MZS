function DF2 = DF2(Z, domain_var)
%DF2  Second derivative F''(Z) of the 'heart' conformal map.
%
%   The map is z = F(Z) = a Z / ((1 - q^2 Z^2)(1 - i g Z)), a = alpha_heart(q,g,pi).
%   Writing F in partial fractions over its poles at Z = 1/q, -1/q and -i/g
%   and differentiating twice gives
%
%     F''(Z) = a [ 2 g^3 / ((q^2 + g^2)(i + g Z)^3)
%                  - q^2 ( 1/((q - i g)(q Z - 1)^3) + 1/((q + i g)(q Z + 1)^3) ) ]
%
%   Z           image-plane position(s), any array size (evaluated elementwise)
%   domain_var  [q, g], the heart shape parameters

q = domain_var(1); g = domain_var(2); a = alpha_heart(q,g,pi);

DF2 = a .* ((2 .* (g.^3) ./ ((q.^2 + g.^2) .* (1i + g.*Z).^3)) ...
            - (q.^2) .* ((1./((q - 1i.*g) .* (q.*Z - 1).^3)) + (1./((q + 1i.*g) .* (q.*Z + 1).^3))));

% unit disc (F(Z) = Z): the formula above is 0/0, but F'' = 0
if (q == 0) && (g == 0)
    DF2 = Z.*0;
end