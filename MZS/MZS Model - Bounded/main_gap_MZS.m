tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% script for bounded MZS code
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% Main parameters %%%%%%%%%%%%%%

%rng('shuffle')                             %random number shuffle
rng(1)
addpath('../Conformal Gap/');              %path to select gap
dt = 0.5;                                  %time interval length
TEND = 300;                                %interval number  
N = 1; NP = 1; NM = N-NP;                  %vortex number, no. +ve, no.-ve  
kappas = [ones(1,NP), -ones(1,NM)]';       %vortex circulations
domain_var = 1;                            %half the size of the gap
parvec = [kappas; domain_var];             %parameters for ODEs
nvor=1e6;                                  %i.c. parameter
delta2=3e-10;                              %position tolerance 

options=odeset('RelTol',delta2, 'AbsTol',delta2, 'Events',@ElongationEvent2); %options for ODE solver

%%%% set up netcdf output %%%%%%%%%%%%%%%%

flname = ' - file.nc';
outname=strcat(pwd,flname);

% restart if file exists, create if not

if exist(outname,'file')==0

    nccreate(outname,'xyXYvar','Dimensions',{'xypos' 4*N 'list' Inf}, ...
           'Format','classic');
    nccreate(outname,'evals','Dimensions',{'list' Inf}, ...
           'Format','classic');

    %phi_0 = rand(N,1)*2*pi;
    phi_0 = ones(N,1)*0;
    %lambda_0 = rand(N,1)*3;
    lambda_0 = ones(N,1)*1;
    X = (lambda_0.^(1/2) - lambda_0.^(-1/2)) .* cos(2*phi_0);
    Y = (lambda_0.^(1/2) - lambda_0.^(-1/2)) .* sin(2*phi_0);
    XYvor0 = vertcat(X,Y);
    %eps = rand(N,1);
    eps = ones(N,1)*0.045;

    x0 = -4; y0 = 0.995; z0 = x0 + 1i*y0;
    Z0 = inverse_map(z0, domain_var); x0 = real(Z0); y0 = imag(Z0);
    xyvor = [x0; y0];
    xyXYvor = vertcat(xyvor, XYvor0);

    %h0 = energy_bounded_MZS(0, xyvor, XYvor0, eps, parvec);
    h0_list = energy_bounded_MZS(0, xyvor, XYvor0, eps, parvec, 1);
    h0 = h0_list(end);                               % initial energy
    errorh_list = zeros(TEND+1, 1); h_list = zeros(TEND+1, length(h0_list));
    h_list(1,:) = h0_list;

    TSTART=1;
    ncwrite(outname,'xyXYvar', xyXYvor, [1,1]);
    ncwrite(outname,'evals', h0, 1);
else
    vinfo = ncinfo(outname,'xyXYvar');
    xyXYvardim=vinfo.Size; T=xyXYvardim(2); TSTART=T;
    xyXYvor=ncread(outname,'xyXYvar',[1,T],[4*N,1],[1,1]);
end

%%%% initial conditions %%%%%%%%%%%%%%%%

xyvor = xyXYvor(1:2*N); XYvor = xyXYvor(2*N+1:4*N);
h0 = energy_bounded_MZS(0, xyvor, XYvor, eps, parvec), hprev = h0;

%%%% main loop %%%%%%%%%%%%%%%%%%%%%%%%%

for TT = TSTART:TEND

    if max(xyXYvor(2*N+1:3*N).^2 + xyXYvor(3*N+1:4*N).^2) >= 25
        break
    end
 
    [tout,xyXYnew] = ...
        ode113(@(t,x) uvUV_bounded_MZS(t, x, eps, parvec), [0,dt], xyXYvor, options);

    xyXYvor = xyXYnew(size(xyXYnew, 1), :)';
    xynew = xyXYvor(1:2*N); XYnew = xyXYvor(2*N+1:4*N);

    %hnew = energy_bounded_MZS(0, xynew, XYnew, eps, parvec);
    hnew_list = energy_bounded_MZS(0, xynew, XYnew, eps, parvec, 1);
    h_list(TT+1,:) = hnew_list;
    hnew = hnew_list(end);
    errorh = abs(hnew - hprev), hprev = hnew;
    errorh_list(TT+1) = errorh;
    TT

    ncwrite(outname,'xyXYvar', xyXYvor, [1,TT+1]);
    ncwrite(outname,'evals', hnew, TT+1);
end

xyXYvar = ncread(outname, 'xyXYvar');
evals = ncread(outname,'evals');

TEND = length(evals) - 1;

xyXYvar = xyXYvar(:, 1:TEND+1); evals = evals(1:TEND+1);

tnow = 0:dt:dt*TEND;

total_h_error = abs(evals(TEND+1) - evals(1))
cumulative_h_error = sum(abs(evals(2:TEND+1) - evals(1:TEND)))

% Plot of H against t
if 1==0
    subplot(1,2,1)
    plot(tnow, evals)
    xlabel('t')
    ylabel('H')
    title('H against t')
end

% Plot of error in H against t
if 1==0
    subplot(1,2,2)
    plot(tnow(2:TEND+1), evals(2:TEND+1) - evals(1:TEND))
    xlabel('t')
    ylabel('Error in H')
    title('Error in H against t')
end

% Plot how much of the energy is in each term of H
if 1==0
    subplot(1,1,1)
    hold on
    plot(tnow, h_list(1:TEND+1, 1:7), LineWidth=1)
    plot(tnow, h_list(1:TEND+1, 8), LineWidth=2)
    hold off
    xlabel('t')
    ylabel('Energy')
    title("Energy in each term of H")
    legend("h1", "h2", "h3", "h4", "h5", "h6", "h7", "H")
end

% Plot of lambda and phi against t
if 1==0
    X_0 = xyXYvar(2*N+1:3*N,:); Y_0 = xyXYvar(3*N+1:4*N,:);
    phi = atan2(Y_0, X_0)/2; absX2 = X_0.^2 + Y_0.^2;
    lambda = ((sqrt(absX2) + sqrt(absX2+4)).^2)/4;
    hold on
    plot(tnow, lambda)
    plot(tnow, phi)
    hold off
    xlabel('t')
    title('lambda, phi')
    legend("lambda", "phi")
end

% Plot of lambda against t
if 1==0
    X_0 = xyXYvar(2*N+1:3*N,:); Y_0 = xyXYvar(3*N+1:4*N,:);
    absX2 = X_0.^2 + Y_0.^2; lambda = ((sqrt(absX2) + sqrt(absX2+4)).^2)/4;
    plot(tnow, lambda)
    xlabel('t')
    title('lambda')
end

% Plot ellipses
if 1==1
    t = 0:0.01:2*pi;
    clear movie
    for k=1:1:TEND+1
        clf
        hold on
        plot([-10, -domain_var], [0, 0], Color='black')
        plot([domain_var, 10], [0, 0], Color='black')
        for n=1:1:N
            if kappas(n) > 0
                col = 'red';
            else
                col = 'blue';
            end
            %plot(xyXYvar(n,1:k), xyXYvar(N+n,1:k), col)
            eps_0 = eps(n); x_0 = xyXYvar(n,k); y_0 = xyXYvar(N+n,k);
            Z_0 = x_0 + 1i.*y_0; z_0 = gap_map(Z_0, domain_var);
            x_0 = real(z_0); y_0 = imag(z_0); X_0 = xyXYvar(2*N+n,k);
            Y_0 = xyXYvar(3*N+n,k); absX2 = X_0^2 + Y_0^2;
            lambda = ((sqrt(absX2) + sqrt(absX2+4))^2)/4; phi = atan2(Y_0,X_0)/2;
            if k==1
                lambda = lambda_0(n); phi = phi_0(n);
            end
            a=sqrt(8*eps_0*lambda);
            b=sqrt(8*eps_0/lambda);
            fill(a*cos(t)*cos(phi) - b*sin(t)*sin(phi) + x_0, a*cos(t)*sin(phi) + b*sin(t)*cos(phi) + y_0, col)
        end
        xlim([-5, 5])
        ylim([-2, 2])
        xlabel('x')
        ylabel('y')
        title('Evolution of Vortices')
        daspect([1,1,1])
        hold off
        movie(k) = getframe;
    end
    writer = VideoWriter('Vortices', 'MPEG-4');
    writer.FrameRate = 20;
    open(writer);
    writeVideo(writer, movie);
    close(writer);
end
toc