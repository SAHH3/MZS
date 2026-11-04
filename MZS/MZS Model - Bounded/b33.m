function b33 = b33(z, zp, domain_var)

zbar = conj(z);
df1 = DF1(zp, domain_var); df2 = DF2(zp, domain_var);

b33 = (2.*zbar ./ (1 - zbar.*zp) - df2./df1) ./ (((1 - zbar.*zp) .* df1).^2);