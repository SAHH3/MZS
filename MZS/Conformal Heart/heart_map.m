function z = heart_map(Z, domain_var)
%HEART_MAP  The 'heart' conformal map from the unit circle to the physical domain.
%
%   z = heart_map(Z, domain_var)
%
%     z = F(Z) = a Z / ((1 - q^2 Z^2)(1 - i g Z)),   a = alpha_heart(q,g,pi)
%
%   The vortex positions are integrated in the image (unit-circle) plane Z;
%   this converts them to physical positions z, e.g. for plotting.
%
%   Z           image-plane position(s), |Z| < 1, any array size
%   domain_var  [q, g], the heart shape parameters (q = g = 0 is the unit circle)

q = domain_var(1); g = domain_var(2); a = alpha_heart(q, g, pi);

z = a.*Z ./ ((1 - (q.*Z).^2) .* (1 - 1i .* g.*Z));