tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% script for bounded MZS code
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% Main parameters %%%%%%%%%%%%%%

%rng('shuffle')                             %random number shuffle
rng(0)
addpath('../Conformal Heart/');             %path to select heart
q = 0.8; g = 0; domain_var = [q; g];        %heart domain parameters
dt = 0.01;                                  %time interval length
TEND = 2500;                                %interval number  
N = 4; NP = N/2; NM = N-NP;                %vortex number, no. +ve, no.-ve  
kappas = [ones(1,NP), -ones(1,NM)]';        %vortex circulations
parvec = [kappas; domain_var];              %parameters for ODEs
nvor=1e6;                                   %i.c. parameter
delta2=1e-9;                               %position tolerance 

options=odeset('RelTol',delta2, 'AbsTol',delta2, 'Events',@ElongationEvent); %options for ODE solver

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
    eps = ones(N,1)*0.0003;

    %eval=(-0.1+0*0.1181249)*N; eshell=0.002*N;
    eval = 0; eshell = 20;

    if g + q == 0
        nccreate(outname,'lvals','Dimensions',{'list' Inf}, ...
               'Format','classic');   %if ang mom
        %lval=sqrt(N/48); lshell=0.01*sqrt(N); %sigma_l=sqrt(N/48)
        lval = 0; lshell = 20;
        xyvor = set_ic_boundedL_MZS(N, NP, eval, eshell, lval, lshell, XYvor0, eps);      % Circle

        modx2 = xyvor(1:N).^2 + xyvor(N+1:2*N).^2; modX2 = X.^2 + Y.^2;
        l1 = sum(kappas .* modx2); l2 = sum(kappas .* 2.*eps .* (modX2 + 2));
        l = l1 + l2; l_list = zeros(TEND+1, 3); l0_list = [l, l1, l2]; l_list(1,:) = l0_list;
        ncwrite(outname,'lvals', l, 1); %if ang mom
    else
        xyvor = set_ic_bounded_MZS(N, nvor, eval, eshell, XYvor0, eps, parvec);       % Other boundary
    end
    xyXYvor = vertcat(xyvor, XYvor0);

    %h0 = energy_bounded_MZS(0, xyvor, XYvor0, eps, parvec);
    h0_list = energy_bounded_MZS(0, xyvor, XYvor0, eps, parvec, 1)
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

if (g == 0) && (q == 0)
    modx2 = xyvor(1:N).^2 + xyvor(N+1:2*N).^2; modX2 = XYvor(1:N).^2 + XYvor(N+1:2*N).^2;
    l0 = sum(kappas .* (modx2 + 2.*eps .* (modX2 + 2))), lprev = l0;
end

%%%% main loop %%%%%%%%%%%%%%%%%%%%%%%%%

for TT = TSTART:TEND

    if max(xyXYvor(2*N+1:3*N).^2 + xyXYvor(3*N+1:4*N).^2) >= 4/3
        break
    end
 
    [tout,xyXYnew] = ...
        ode113(@(t,x) uvUV_bounded_MZS(t, x, eps, parvec), [0,dt], xyXYvor, options);
    xyXYvor = xyXYnew(size(xyXYnew, 1), :)';
    xynew = xyXYvor(1:2*N); XYnew = xyXYvor(2*N+1:4*N);

    if (g == 0) && (q == 0)
        modx2 = xynew(1:N).^2 + xynew(N+1:2*N).^2; modX2 = XYnew(1:N).^2 + XYnew(N+1:2*N).^2;
        l1new = sum(kappas .* modx2); l2new = sum(kappas .* 2.*eps .* (modX2 + 2));
        lnew = l1new + l2new; l_list(TT+1,:) = [lnew, l1new, l2new];
        errorl = abs(lnew - lprev), lprev = lnew;
        ncwrite(outname,'lvals', lnew, TT+1); %if ang mom
    end

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

TEND = length(evals) - 2;

xyXYvar = xyXYvar(:, 1:TEND+1); evals = evals(1:TEND+1);

tnow = 0:dt:dt*TEND;

total_h_error = abs(evals(TEND+1) - evals(1))
cumulative_h_error = sum(abs(evals(2:TEND+1) - evals(1:TEND)))

if g + q == 0
    lvals = ncread(outname,'lvals'); lvals = lvals(1:TEND+1);
    total_l_error = abs(lvals(TEND+1) - lvals(1))
    cumulative_l_error = sum(abs(lvals(2:TEND+1) - lvals(1:TEND)))
end

[h0_list(1), h0_list(2) + h0_list(3) + h0_list(6), h0_list(4) + h0_list(5) + h0_list(7), h0_list(8)]
h0
max(abs(evals-h0))
tnow(end)

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

% Plot of L against t
if (g + q == 0) && (1==0)
    subplot(2,2,3)
    plot(tnow, lvals)
    xlabel('t')
    ylabel('L')
    title('L against t')
end

% Plot of error in L against t
if (g + q == 0) && (1==0)
    subplot(2,2,4)
    plot(tnow(2:TEND+1), lvals(2:TEND+1) - lvals(1:TEND))
    xlabel('t')
    ylabel('Error in L')
    title('Error in L against t')
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

% Plot of energy change in each term of H
if 1==1
    hold on
    plot(tnow, h_list(1:TEND+1, 1) - h0_list(1), LineWidth=1)
    plot(tnow, h_list(1:TEND+1, 2) + h_list(1:TEND+1, 3) + h_list(1:TEND+1, 6) - h0_list(2) - h0_list(3) - h0_list(6), LineWidth=1)
    plot(tnow, h_list(1:TEND+1, 4) + h_list(1:TEND+1, 5) + h_list(1:TEND+1, 7) - h0_list(4) - h0_list(5) - h0_list(7), LineWidth=1)
    plot(tnow, h_list(1:TEND+1, 8) - h0_list(8), LineWidth=1)
    hold off
    set(gca, 'fontsize', 20)
    xlabel('Time t')
    ylabel('Energy')
    %title("Energy change in each term of H")
    lgd = legend("\Delta H_1: Self-Interaction", "\Delta H_2: Point Vortex", "\Delta H_3: Elliptic", "\Delta H: Total");
    fontsize(lgd, 20, 'points')
end

% Plot how much of the angular momentum is in each term of L
if (g + q == 0) && (1==0)
    subplot(1,2,1)
    plot(tnow, l_list(1:TEND+1,1:2), LineWidth=1)
    xlabel('t')
    ylabel('Angular Momentum')
    title("Angular Momentum in each term of L")
    legend("L", "Point Vortex")

    subplot(1,2,2)
    plot(tnow, l_list(1:TEND+1,3), LineWidth=1)
    xlabel('t')
    ylabel('Angular Momentum')
    title("Angular Momentum in each term of L")
    legend("Ellipticity")
end

% Plot of change in each term of L
if (g + q == 0) && (1==0)
    hold on
    plot(tnow, l_list(1:TEND+1, 1) - l0_list(1), LineWidth=1)
    plot(tnow, l_list(1:TEND+1, 2) - l0_list(2), LineWidth=1)
    plot(tnow, l_list(1:TEND+1, 3) - l0_list(3), LineWidth=1)
    hold off
    set(gca, 'fontsize', 15)
    xlabel('t')
    ylabel('Change in Angular Momentum')
    %title("Change in Angular Momentum in each term of L")
    lgd = legend("Change in L", "Change in Point Vortex Angular Momentum", "Change in Ellipticity Angular Momentum");
    fontsize(lgd, 15, 'points')
end

% Plot ellipses
if 1==0
    t = 0:0.01:2*pi;
    x_list = zeros(N, TEND+1); y_list = zeros(N, TEND+1);
    clear movie
    for k=1:1:TEND+1
        clf
        hold on
        if g + q == 0
            plot(cos(t), sin(t), 'black');
            xlim([-1, 1])
            ylim([-1, 1])
        else
            [xs, ys] = boundary(1000, domain_var);
            plot(xs, ys, 'black', LineWidth=2)
        end
        for n=1:1:N
            if kappas(n) > 0
                col = 'red';
            else
                col = 'blue';
            end
            x0 = xyXYvar(n,1); y0 = xyXYvar(N+n,1);
            eps_0 = eps(n); x_0 = xyXYvar(n,k); y_0 = xyXYvar(N+n,k);
            if g + q > 0
                Z0 = x0 + 1i.*y0; z0 = heart_map(Z0, domain_var); x0 = real(z0); y0 = imag(z0);
                Z_0 = x_0 + 1i.*y_0; z_0 = heart_map(Z_0, domain_var); x_0 = real(z_0); y_0 = imag(z_0);
            end
            scatter(x0, y0, col, "filled")

            x_list(n, k) = x_0; y_list(n, k) = y_0;
            plot(x_list(n, 1:k), y_list(n, 1:k), col)
            X_0 = xyXYvar(2*N+n,k); Y_0 = xyXYvar(3*N+n,k); absX2 = X_0^2 + Y_0^2;
            lambda = ((sqrt(absX2) + sqrt(absX2+4))^2)/4;
            phi = atan2(Y_0,X_0)/2;
            if k==1
                lambda = lambda_0(n); phi = phi_0(n);
            end
            a=sqrt(8*eps_0*lambda);
            b=sqrt(8*eps_0/lambda);
            fill(a*cos(t)*cos(phi) - b*sin(t)*sin(phi) + x_0, a*cos(t)*sin(phi) + b*sin(t)*cos(phi) + y_0, col, FaceAlpha=0.5)
        end
        %xlim([-1.2,1.2])
        %ylim([-1.2,1.2])
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