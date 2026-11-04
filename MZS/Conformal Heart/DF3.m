function DF3 = DF3(Z, domain_var)

% DF3:
% returns the value of F'''(Z)
% of the 'heart' mapping z = F(Z) = a(q,g)Z / ((1 - (qZ)^2) (1 - igZ))
% Note that Z will typically be vector valued

q = domain_var(1); g = domain_var(2); a = alpha_heart(q,g,pi);

DF3 = 3.*a .* (-(2 .* (g.^4) ./ ((q.^2 + g.^2) .* (1i + g.*Z).^4)) ...
               + (q.^3) .* ((1./((q - 1i.*g) .* (q.*Z - 1).^4)) + (1./((q + 1i.*g) .* (q.*Z + 1).^4))));

if (q == 0) && (g == 0)
    DF3 = Z.*0;
end