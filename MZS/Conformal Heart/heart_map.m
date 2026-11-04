function z = heart_map(Z, domain_var)

% Conformal 'heart' mapping z = F(Z) = a(q,g) Z / (1 - q^2 Z^2)(1 - igZ)

q = domain_var(1); g = domain_var(2); a = alpha_heart(q, g, pi);

z = a.*Z ./ ((1 - (q.*Z).^2) .* (1 - 1i .* g.*Z));