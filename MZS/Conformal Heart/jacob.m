function [ jacb ] = jacob( bigz, domain_var )

% JACOB: 
% returns the value of the jacobian |F'(Z)| 
% of the `heart' mapping z=F(Z)=a(q,g) Z/(1-(qZ)^2)(1-I g Z)
% note that bigz will typically be vector valued

q=domain_var(1); g=domain_var(2); a=alpha_heart(q,g,pi);
jacb=abs(a*(1+q^2*bigz.^2-2*1i*q^2*g*bigz.^3)   ... 
     ./((1-q^2*bigz.^2).^2.*(1-1i*g*bigz).^2));

end

