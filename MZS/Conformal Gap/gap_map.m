function z = gap_map(Z, a)

% Conformal 'gap' mapping z = F(Z) = 2iZ / (Z^2 - 1)

z = 2i.*Z ./ (Z.^2 - 1);