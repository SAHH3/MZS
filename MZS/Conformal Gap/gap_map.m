function z = gap_map(Z, a)
%GAP_MAP  The 'gap' conformal map from the unit circle to the physical domain.
%
%   z = gap_map(Z, a)
%
%     z = F(Z) = 2 i Z / (Z^2 - 1)
%
%   The physical domain is the whole plane with two straight walls along the
%   real axis, x <= -1 and x >= 1, leaving a gap -1 < x < 1 between them.
%   The boundary |Z| = 1 maps onto the walls, Z = 0 maps to the centre of the gap,
%   and Z = +1, -1 map to infinity (one on each side of the walls).
%
%   The vortex positions are integrated in the image (unit-circle) plane Z;
%   this converts them to physical positions z, e.g. for plotting.
%
%   Z  image-plane position(s), |Z| < 1, any array size
%   a  not used (the gap half-width is fixed at 1); kept so that all the
%      conformal-map functions take the same (Z, domain_var) arguments

z = 2i.*Z ./ (Z.^2 - 1);