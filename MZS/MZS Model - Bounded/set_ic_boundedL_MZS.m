function [ xyvor ] = set_ic_boundedL_MZS(n, np, eval, eshell, lval, lshell, XY, eps)
%SET_IC_BOUNDEDL_MZS  Random initial positions in the UNIT CIRCLE with given
%                     energy and angular momentum.
%
%   xyvor = set_ic_boundedL_MZS(n, np, eval, eshell, lval, lshell, XY, eps)
%
%   For use with the circle only (q = g = 0 in main_bounded_MZS.m).
%   Finds positions satisfying
%       eval - eshell/2 < E < eval + eshell/2
%   AND lval - lshell/2 < L < lval + lshell/2, where the angular momentum is
%       L = sum_i k_i ( |Z_i|^2 + 2 eps_i (|X_i|^2 + 2) ).
%   It first draws random radii until L is in range, then keeps those radii
%   and draws random angles until E is in range.
%   NB. this means the positions are NOT drawn from a uniform distribution.
%
%   n, np           number of vortices; the first np have k = +1, the rest -1
%   eval, eshell    centre and width of the energy shell
%   lval, lshell    centre and width of the angular momentum shell
%   XY, eps         vortex shapes [X; Y] (2n x 1) and sizes (n x 1)
%   xyvor           [x; y], 2n x 1

    nm=n-np; l=-1e10; hvec=-1e10;
    kappas=[ones(1,np),-ones(1,nm)]';
    parvec=[kappas;0;0];                 % q = g = 0: unit circle

    while (l<lval-0.5*lshell || l>lval+0.5*lshell)    %first fix L
       rads = sqrt(rand(n,1));           % sqrt gives a uniform density over the unit circle
       l = sum(kappas .* (rads.^2 + 2.*eps .* (XY(1:n).^2 + XY(n+1:2*n).^2 + 2)));
    end

    while (hvec<eval-eshell/2 || hvec>eval+eshell/2) %next fix E
        angles = 2*pi*rand(n,1);
        xyvor = [rads.*cos(angles); rads.*sin(angles)];
        hvec = energy_bounded_MZS(0, xyvor, XY, eps, parvec);
    end

end