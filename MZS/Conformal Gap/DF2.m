function DF2 = DF2(Z, a)
%DF2  Second derivative F''(Z) of the 'gap' conformal map z = F(Z) = 2iZ / (Z^2 - 1).
%
%     F''(Z) = 4i Z (Z^2 + 3) / (Z^2 - 1)^3
%
%   Z  image-plane position(s), any array size (evaluated elementwise)
%   a  not used (gap half-width fixed at 1); see gap_map.m

DF2 = 4i.*Z .* (Z.^2 + 3) ./ ((Z.^2 - 1).^3);