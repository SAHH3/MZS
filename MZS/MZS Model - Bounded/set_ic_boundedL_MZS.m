function [ xyvor ] = set_ic_boundedL_MZS(n, np, eval, eshell, lval, lshell, XY, eps)

% SET_IC_BOUNDEDL:
% function to generate random i.c.s for main_bounded_MZS
% for use with CIRCLE ONLY!!
% satisfying eval-eshell/2 < E <eval+eshell/2
% AND        lval-lshell/2 < L < lval+lshell/2
% n    - no. vortices, np - no. positive
% eval, eshell, lval, lshell - centres and thicknesses of (E,L) shell.
% first fixes radii to get angular momentum correct, 
% then perturbs angles to get energy correct
% nb. this means we are NOT drawing from uniform distribution 

    nm=n-np; l=-1e10; hvec=-1e10; 
    kappas=[ones(1,np),-ones(1,nm)]';
    parvec=[kappas;0;0];

    while (l<lval-0.5*lshell || l>lval+0.5*lshell)    %first fix L
       rads = sqrt(rand(n,1));
       l = sum(kappas .* (rads.^2 + 2.*eps .* (XY(1:n).^2 + XY(n+1:2*n).^2 + 2)));
    end       

    while (hvec<eval-eshell/2 || hvec>eval+eshell/2) %next fix E         
        angles = 2*pi*rand(n,1);
        xyvor = [rads.*cos(angles); rads.*sin(angles)];
        hvec = energy_bounded_MZS(0, xyvor, XY, eps, parvec);
    end

end
