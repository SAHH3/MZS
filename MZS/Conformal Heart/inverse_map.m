function [ bigz ] = inverse_map( z, domain_var )
%INVERSE_MAP  Inverse of the 'heart' conformal map (physical heart -> unit disc).
%
%   bigz = inverse_map(z, domain_var)
%
%   Solves z = a Z / ((1 - q^2 Z^2)(1 - i g Z)) for Z. Rearranged, this is a
%   cubic in Z, which is solved with Cardano's formula; all three roots are
%   formed and only those with |Z| < 1 are kept.
%
%   The input z may contain points outside the heart (xyvors.m samples a box
%   around it): those have no root inside the disc and are simply dropped, so
%   bigz can be SHORTER than z, and its entries no longer line up with z.
%
%   Not valid for the unit disc itself (q = g = 0), where it divides by zero.
%
%   z           physical position(s), row vector
%   domain_var  [q, g], the heart shape parameters

q=domain_var(1); g=domain_var(2); a=alpha_heart(q,g,pi);

if (g>0) %get roots of cubic (Cardano: bigs, bigt are the two cube-root terms)

    bigr=(-9*a*g + 2*1i*(9*g^2 + q^2)*z)./(54*g^3*q^2*z);
    bigq=(3*1i*a*g-3*g^2*z+q^2*z)./(9*g^2*q^2*z);
    bigs=(bigr+(bigq.^3+bigr.^2).^(1/2)).^(1/3);
    bigt=-bigq./bigs;
    bigzall=[bigs+bigt-1i/(3*g), ...
        -0.5*(bigs+bigt)-0.5*sqrt(3)*1i*(bigs-bigt)-1i/(3*g), ...
        -0.5*(bigs+bigt)+0.5*sqrt(3)*1i*(bigs-bigt)-1i/(3*g)];

else     %cubic replaced by quadratic for g=0: q^2 z Z^2 + a Z - z = 0

    bigzall=((a^2+4*q^2*z.^2).^0.5-a)./(2*q^2*z);

end

% keep only the roots that lie inside the unit disc
bigz=bigzall(abs(bigzall)<1);

end