function Z = inverse_map(z, a)

% INVERSE_MAP:
% Computes the inverse of the gap conformal map (i.e. gap->circle)

Z = (1i - z - sqrt(z.^2 - 1)) ./ (1i + z + sqrt(z.^2 - 1));
for j=1:1:length(z)
    if abs(Z(j)) > 1
        Z(j) = (1i - z(j) + sqrt(z(j).^2 - 1)) ./ (1i + z(j) - sqrt(z(j).^2 - 1));
    end
end