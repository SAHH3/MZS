function DF2 = DF2(Z, a)

% DF2:
% returns the value of F''(Z)
% of the 'gap' mapping z = F(Z) = 2iZ / (Z^2 - 1)
% Note that Z will typically be vector valued

DF2 = 4i.*Z .* (Z.^2 + 3) ./ ((Z.^2 - 1).^3);