function a = a_complex(z, zp)

a = 1 - 2.*real(z .* conj(zp)) + (abs(z) .* abs(zp)).^2;