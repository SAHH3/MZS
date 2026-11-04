function DF3 = DF3(Z, a)
%DF3  Third derivative F'''(Z) of the 'gap' conformal map z = F(Z) = 2iZ / (Z^2 - 1).
%
%     F'''(Z) = -12i (Z^4 + 6 Z^2 + 1) / (Z^2 - 1)^4
%
%   Z  image-plane position(s), any array size (evaluated elementwise)
%   a  not used (gap half-width fixed at 1); see gap_map.m

DF3 = -12i .* (Z.^4 + 6.*(Z.^2) + 1) ./ ((Z.^2 - 1).^4);