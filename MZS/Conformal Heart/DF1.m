function DF1 = DF1(Z, domain_var)
%DF1  First derivative F'(Z) of the 'heart' conformal map.
%
%   The map is z = F(Z) = a Z / ((1 - q^2 Z^2)(1 - i g Z)), a = alpha_heart(q,g,pi),
%   which takes the unit disc |Z| < 1 (image plane) to the heart (physical plane).
%
%     F'(Z) = a (1 + q^2 Z^2 - 2 i g q^2 Z^3) / ((1 - q^2 Z^2)(1 - i g Z))^2
%
%   Z           image-plane position(s), any array size (evaluated elementwise)
%   domain_var  [q, g], the heart shape parameters

q = domain_var(1); g = domain_var(2); a = alpha_heart(q,g,pi);

DF1 = a .* (1 + (q.^2) .* (Z.^2) - 2i .*g .* (q.^2) .* (Z.^3)) ./ ...
           (((1 - (q.^2) .* (Z.^2)) .* (1 - 1i.*g .* Z)).^2);