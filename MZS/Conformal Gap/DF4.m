function DF4 = DF4(Z, a)
%DF4  Fourth derivative F''''(Z) of the 'gap' conformal map z = F(Z) = 2iZ / (Z^2 - 1).
%
%     F''''(Z) = 48i Z (Z^4 + 10 Z^2 + 5) / (Z^2 - 1)^5
%
%   Z  image-plane position(s), any array size (evaluated elementwise)
%   a  not used (gap half-width fixed at 1); see gap_map.m

DF4 = 48i.*Z .* (Z.^4 + 10.*(Z.^2) + 5) ./ ((Z.^2 - 1).^5);