function DF4 = DF4(Z, a)

% DF4:
% returns the value of F''''(Z)
% of the 'gap' mapping z = F(Z) = 2iZ / (Z^2 - 1)
% Note that Z will typically be vector valued

DF4 = 48i.*Z .* (Z.^4 + 10.*(Z.^2) + 5) ./ ((Z.^2 - 1).^5);