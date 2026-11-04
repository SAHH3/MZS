function [ xyvor ] = set_ic_bounded_MZS(n, nvor, eval, eshell, XY, eps, parvec)

% SET_IC_BOUNDED:
% function to generate random i.c.s for main_bounded_MZS
% satisfying eval-eshell/2<E<eval+eshell/2
% n    - no. vortices
% eval, ehell - centre and thickness of energy shell.
% nvor - total number of possible vortex positons in record.

    xyvor=[]; domain_vec=parvec(n+1:length(parvec))';
    bigzs=xyvors(nvor,domain_vec);
    nrec=floor(nvor/n); NGO=1;

    while (size(xyvor,1)==0)
        NGO
        i1=randperm(nvor); 

        for ii=1:nrec               

            xys = [real(bigzs(i1((ii-1)*n+1:ii*n))), ...
                 imag(bigzs(i1((ii-1)*n+1:ii*n)))]';
            hvec = energy_bounded_MZS(0, xys, XY, eps, parvec);

            if (hvec>eval-eshell/2 && hvec<eval+eshell/2)
                xyvor=xys; 
                return 
            end
        end

	NGO=NGO+1;
    end
end
