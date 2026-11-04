function DF3 = DF3(Z, a)

% DF3:
% returns the value of F'''(Z)
% of the 'gap' mapping z = F(Z) = 2iZ / (Z^2 - 1)
% Note that Z will typically be vector valued

DF3 = -12i .* (Z.^4 + 6.*(Z.^2) + 1) ./ ((Z.^2 - 1).^4);