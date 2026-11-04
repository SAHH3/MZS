
function [ bigz ] = inverse_map( z, domain_var )
%
% INVERSE_MAP:
% Computes the inverse of the heart conformal map (i.e. heart->circle)
% Input vector z=x+iy can include points outside the heart 
% (for use with set_ic_heart), 
% but output vector bigz consists of only those points inside the circle.
% 
q=domain_var(1); g=domain_var(2); a=alpha_heart(q,g,pi);

if (g>0) %get roots of cubic
    
    bigr=(-9*a*g + 2*1i*(9*g^2 + q^2)*z)./(54*g^3*q^2*z);
    bigq=(3*1i*a*g-3*g^2*z+q^2*z)./(9*g^2*q^2*z);
    bigs=(bigr+(bigq.^3+bigr.^2).^(1/2)).^(1/3);
    bigt=-bigq./bigs;
    bigzall=[bigs+bigt-1i/(3*g), ...
        -0.5*(bigs+bigt)-0.5*sqrt(3)*1i*(bigs-bigt)-1i/(3*g), ...
        -0.5*(bigs+bigt)+0.5*sqrt(3)*1i*(bigs-bigt)-1i/(3*g)];

else     %cubic replaced by quadratic for g=0
    
    bigzall=((a^2+4*q^2*z.^2).^0.5-a)./(2*q^2*z);
    
end
        
bigz=bigzall(abs(bigzall)<1);

end

