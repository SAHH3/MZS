function DF4 = DF4(Z, domain_var)

% DF4:
% returns the value of F''''(Z)
% of the 'heart' mapping z = F(Z) = a(q,g)Z / ((1 - (qZ)^2) (1 - igZ))
% Note that Z will typically be vector valued

q = domain_var(1); g = domain_var(2); a = alpha_heart(q,g,pi);

DF4 = 12.*a .* ((2 .* (g.^5) ./ ((q.^2 + g.^2) .* (1i + g.*Z).^5)) ...
                - (q.^4) .* ((1./((q - 1i.*g) .* (q.*Z - 1).^5)) + (1./((q + 1i.*g) .* (q.*Z + 1).^5))));

if (q == 0) && (g == 0)
    DF4 = Z.*0;
end