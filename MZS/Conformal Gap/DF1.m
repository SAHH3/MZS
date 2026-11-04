function DF1 = DF1(Z, a)

% DF1:
% returns the value of F'(Z)
% of the 'gap' mapping z = F(Z) = 2iZ / (Z^2 - 1)
% Note that Z will typically be vector valued

DF1 = -2i .* (Z.^2 + 1) ./ ((Z.^2 - 1).^2);