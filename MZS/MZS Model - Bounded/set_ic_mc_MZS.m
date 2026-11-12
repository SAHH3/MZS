function [ xyvor ] = set_ic_mc_MZS(n, nvor, emin, XY, eps, parvec)
%SET_IC_MC_MZS  Random initial positions in the heart with energy above emin.
%
%   xyvor = set_ic_mc_MZS(n, nvor, emin, XY, eps, parvec)
%
%   An alternative to set_ic_bounded_MZS.m for high energies, where random
%   groups rarely land in the shell. Starts from n random candidate positions
%   (xyvors.m) and repeatedly tries replacing one vortex with another
%   candidate, keeping the change only if it raises the energy, until E > emin.
%
%   n        number of vortices
%   nvor     number of candidate positions to generate
%   emin     energy to exceed
%   XY, eps  vortex shapes [X; Y] (2n x 1) and sizes (n x 1)
%   parvec   [kappas; domain_var]
%   xyvor    [x; y], 2n x 1, image-plane positions
%   Needs boundary.m, so works with the Conformal Heart folder only.

    domain_vec=parvec(n+1:end)';
    bigzs=xyvors(nvor,domain_vec);       % candidate image-plane positions

    inds = randperm(nvor,n); % This is the running list or vortices
    xys = [real(bigzs(inds)), imag(bigzs(inds))]';
    hvec = energy_bounded_MZS(0, xys, XY, eps, parvec);

    while (hvec < emin)
        % Attempt to swap vortices i1 <--> i2
        i0 = randperm(n,1);
        i1 = inds(i0);
        i2 = randperm(nvor,1);
        if( ismember(i2,inds) )          % candidate already in use
            continue
        end

        inds(i0) = i2;
        xys = [real(bigzs( inds )) ...
               imag(bigzs( inds ))]';
        hnew = energy_bounded_MZS(0, xys, XY, eps, parvec);
        if (hnew>hvec) % Accept swap
            hvec = hnew;
        else
            inds(i0) = i1; % Reject swap
        end
    end
    xyvor = xys;
end