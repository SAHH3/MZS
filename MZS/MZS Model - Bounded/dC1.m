function dC1 = dC1(z, domain_var)

df1 = DF1(z, domain_var); df2 = DF2(z, domain_var); df3 = DF3(z, domain_var); df4 = DF4(z, domain_var);

dC1 = -df4 ./ (3 .* (df1.^3)) + 5.*df3.*df2 ./ (2 .* (df1.^4)) - 3 .* (df2.^3) ./ (df1.^5);