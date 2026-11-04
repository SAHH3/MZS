function b44 = b44(z, zp, domain_var)

zpbar = conj(zp);
df1 = DF1(z, domain_var); df2 = DF2(z, domain_var); df3 = DF3(z, domain_var);

b44 = (2 .* (zpbar ./ (1 - z.*zpbar)).^2 - 3.*zpbar.*df2 ./ (df1 .* (1 - z.*zpbar)) + 3.*(df2./df1).^2 - df3./df1) .* zpbar ./ ((1 - z.*zpbar) .* (df1.^2));