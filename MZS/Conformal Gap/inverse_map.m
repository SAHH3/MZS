function Z = inverse_map(z, a)
%INVERSE_MAP  Inverse of the 'gap' conformal map (physical plane -> unit circle).
%
%   Z = inverse_map(z, a)
%
%   Solves z = 2iZ / (Z^2 - 1), i.e. z Z^2 - 2i Z - z = 0, whose roots are
%   Z = (i -/+ sqrt(z^2 - 1)) / z. They are written here in the equivalent form
%   (i - z -/+ s) / (i + z +/- s), s = sqrt(z^2 - 1), which is finite at z = 0.
%   The two roots multiply to -1, so exactly one lies inside the unit circle:
%   the first root is tried, and swapped for the second wherever |Z| > 1.
%
%   z  physical position(s), not on the walls; Z has the same size as z
%   a  not used (gap half-width fixed at 1); see gap_map.m

Z = (1i - z - sqrt(z.^2 - 1)) ./ (1i + z + sqrt(z.^2 - 1));
for j=1:1:length(z)
    if abs(Z(j)) > 1
        Z(j) = (1i - z(j) + sqrt(z(j).^2 - 1)) ./ (1i + z(j) - sqrt(z(j).^2 - 1));
    end
end