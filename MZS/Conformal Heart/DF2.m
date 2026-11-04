function DF2 = DF2(Z, domain_var)

% DF2:
% returns the value of F''(Z)
% of the 'heart' mapping z = F(Z) = a(q,g)Z / ((1 - (qZ)^2) (1 - igZ))
% Note that Z will typically be vector valued

q = domain_var(1); g = domain_var(2); a = alpha_heart(q,g,pi);

DF2 = a .* ((2 .* (g.^3) ./ ((q.^2 + g.^2) .* (1i + g.*Z).^3)) ...
            - (q.^2) .* ((1./((q - 1i.*g) .* (q.*Z - 1).^3)) + (1./((q + 1i.*g) .* (q.*Z + 1).^3))));

if (q == 0) && (g == 0)
    DF2 = Z.*0;
end