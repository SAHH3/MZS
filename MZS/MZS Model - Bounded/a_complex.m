function a = a_complex(z, zp)
%A_COMPLEX  The image-vortex factor |1 - Z conj(Z')|^2 for the unit circle, complex form.
%
%   a = a_complex(z, zp)
%
%     a = |1 - Z conj(Z')|^2 = 1 - 2 Re(Z conj(Z')) + |Z|^2 |Z'|^2
%
%   Same quantity as a.m, but taking complex image-plane positions.
%   z, zp  image-plane positions Z, Z'; a column and a row give the N x N
%          matrix for every pair (as used in uvUV_bounded_MZS.m)

a = 1 - 2.*real(z .* conj(zp)) + (abs(z) .* abs(zp)).^2;