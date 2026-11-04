function [ xyvor ] = set_ic_mc_MZS(n, nvor, emin, XY, eps, parvec)

% SET_IC_BOUNDED:
% function to generate random i.c.s for main_bounded_MZS
% satisfying emin<E
% n    - no. vortices
% nvor - total number of possible vortex positons in record.

    domain_vec=parvec(n+1:end)';
    bigzs=xyvors(nvor,domain_vec);

    inds = randperm(nvor,n); % This is the running list or vortices
    xys = [real(bigzs(inds)), imag(bigzs(inds))]';
    hvec = energy_bounded_MZS(0, xys, XY, eps, parvec);

    while (hvec < emin)
        % Attempt to swap vortices i1 <--> i2
        i0 = randperm(n,1); 
        i1 = inds(i0); 
        i2 = randperm(nvor,1);
        if( ismember(i2,inds) )
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
