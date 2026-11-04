function BB2 = BB2(z, zp, domain_var)

zbar = conj(z);
df1 = DF1(zp, domain_var); df2 = DF2(zp, domain_var);

BB2 = (zbar./(1 - zbar.*zp) - df2./df1) .* zbar ./ ((1 - zbar.*zp) .* (df1.^2));