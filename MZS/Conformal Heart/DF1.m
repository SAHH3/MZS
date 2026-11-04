function DF1 = DF1(Z, domain_var)

% DF1:
% returns the value of F'(Z)
% of the 'heart' mapping z = F(Z) = a(q,g)Z / ((1 - (qZ)^2) (1 - igZ))
% Note that Z will typically be vector valued

q = domain_var(1); g = domain_var(2); a = alpha_heart(q,g,pi);

DF1 = a .* (1 + (q.^2) .* (Z.^2) - 2i .*g .* (q.^2) .* (Z.^3)) ./ ...
           (((1 - (q.^2) .* (Z.^2)) .* (1 - 1i.*g .* Z)).^2);