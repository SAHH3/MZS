function DF1 = DF1(Z, a)
%DF1  First derivative F'(Z) of the 'gap' conformal map z = F(Z) = 2iZ / (Z^2 - 1).
%
%     F'(Z) = -2i (Z^2 + 1) / (Z^2 - 1)^2
%
%   Z  image-plane position(s), any array size (evaluated elementwise)
%   a  not used (gap half-width fixed at 1); see gap_map.m

DF1 = -2i .* (Z.^2 + 1) ./ ((Z.^2 - 1).^2);