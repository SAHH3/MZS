%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MAIN_BOUNDED_MZS: script for bounded MZS code
%
% Integrates N elliptical vortices (MZS moment model) in a bounded domain
% given by a conformal map F of the unit circle: here the 'heart' map from the
% Conformal Heart folder, which with q = g = 0 is the unit circle itself.
%
% State vector (4N x 1): [x; y; X; Y]
%   x + i y  image-plane (unit-circle) position Z of each vortex centre;
%            the physical position is F(Z) (heart_map.m)
%   X + i Y  shape, (sqrt(lambda) - 1/sqrt(lambda)) exp(2 i phi), for an
%            ellipse of aspect ratio lambda and orientation phi, with
%            semi-axes sqrt(8 eps lambda), sqrt(8 eps / lambda) (area 8 pi eps)
% Equations of motion: uvUV_bounded_MZS.m. Energy: energy_bounded_MZS.m.
%
% The state after every step of length dt is appended to a NetCDF file
% (outname, saved in the folder ABOVE this one). If that file already exists,
% the run continues from its last saved state, so delete it to start afresh.
% The run stops early if any vortex reaches aspect ratio 3 (ElongationEvent.m).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% Main parameters %%%%%%%%%%%%%%

%rng('shuffle')                             %random number shuffle
rng(0)
addpath('../Conformal Heart/');             %path to select heart
q = 0.8; g = 0; domain_var = [q; g];        %heart domain parameters (q = g = 0: unit circle)
dt = 0.01;                                  %time interval length
TEND = 2500;                                %number of time steps (run length TEND*dt)
N = 4; NP = N/2; NM = N-NP;                %vortex number, no. +ve, no.-ve
kappas = [ones(1,NP), -ones(1,NM)]';        %vortex circulations
parvec = [kappas; domain_var];              %parameters for ODEs
nvor=1e6;                                   %no. of candidate positions for set_ic_bounded_MZS
delta2=1e-9;                               %ODE113 relative and absolute tolerance

% initial ellipse shapes and sizes (set here so they also exist on a restart)
% phi_0: orientations, lambda_0: aspect ratios, eps: sizes
%phi_0 = rand(N,1)*2*pi;
phi_0 = ones(N,1)*0;
%lambda_0 = rand(N,1)*3;
lambda_0 = ones(N,1)*1;
X = (lambda_0.^(1/2) - lambda_0.^(-1/2)) .* cos(2*phi_0);
Y = (lambda_0.^(1/2) - lambda_0.^(-1/2)) .* sin(2*phi_0);
XYvor0 = vertcat(X,Y);
%eps = rand(N,1);
eps = ones(N,1)*0.0003;

options=odeset('RelTol',delta2, 'AbsTol',delta2, 'Events',@ElongationEvent); %options for ODE solver

%%%% set up netcdf output %%%%%%%%%%%%%%%%

% variables in the file: xyXYvar (state after each step, 4N x steps),
% evals (energy after each step) and, for the unit circle, lvals (angular momentum)
flname = ' - file.nc';
outname=strcat(pwd,flname);                 %i.e. '<this folder's path> - file.nc'

% restart if file exists, create if not

if exist(outname,'file')==0

    nccreate(outname,'xyXYvar','Dimensions',{'xypos' 4*N 'list' Inf}, ...
           'Format','classic');
    nccreate(outname,'evals','Dimensions',{'list' Inf}, ...
           'Format','classic');

    % initial positions are drawn at random with energy in [eval - eshell/2, eval + eshell/2]
    %eval=(-0.1+0*0.1181249)*N; eshell=0.002*N;
    eval = 0; eshell = 20;

    if g + q == 0
        % unit circle: angular momentum L is also conserved, so fix it too
        nccreate(outname,'lvals','Dimensions',{'list' Inf}, ...
               'Format','classic');   %if ang mom
        %lval=sqrt(N/48); lshell=0.01*sqrt(N); %sigma_l=sqrt(N/48)
        lval = 0; lshell = 20;
        xyvor = set_ic_boundedL_MZS(N, NP, eval, eshell, lval, lshell, XYvor0, eps);      % Circle

        % L = l1 + l2: point-vortex part l1 = sum k |Z|^2, ellipse part l2 = sum 2 k eps (|X|^2 + 2)
        modx2 = xyvor(1:N).^2 + xyvor(N+1:2*N).^2; modX2 = X.^2 + Y.^2;
        l1 = sum(kappas .* modx2); l2 = sum(kappas .* 2.*eps .* (modX2 + 2));
        l = l1 + l2; l_list = zeros(TEND+1, 3); l0_list = [l, l1, l2]; l_list(1,:) = l0_list;
        ncwrite(outname,'lvals', l, 1); %if ang mom
    else
        xyvor = set_ic_bounded_MZS(N, nvor, eval, eshell, XYvor0, eps, parvec);       % Other boundary
    end
    xyXYvor = vertcat(xyvor, XYvor0);

    % energy, split into its terms [h1, ..., h7, H] (see energy_bounded_MZS.m)
    %h0 = energy_bounded_MZS(0, xyvor, XYvor0, eps, parvec);
    h0_list = energy_bounded_MZS(0, xyvor, XYvor0, eps, parvec, 1)
    h0 = h0_list(end);                               % initial energy
    errorh_list = zeros(TEND+1, 1); h_list = zeros(TEND+1, length(h0_list));
    h_list(1,:) = h0_list;

    % save the initial state as column 1
    TSTART=1;
    ncwrite(outname,'xyXYvar', xyXYvor, [1,1]);
    ncwrite(outname,'evals', h0, 1);
else
    % file exists: continue from its last saved state (column T)
    vinfo = ncinfo(outname,'xyXYvar');
    xyXYvardim=vinfo.Size; T=xyXYvardim(2); TSTART=T;
    xyXYvor=ncread(outname,'xyXYvar',[1,T],[4*N,1],[1,1]);

    % rebuild the energy (and angular momentum) diagnostics for the steps already saved
    xyXYsaved = ncread(outname,'xyXYvar',[1,1],[4*N,T]);
    h0_list = energy_bounded_MZS(0, xyXYsaved(1:2*N,1), xyXYsaved(2*N+1:4*N,1), eps, parvec, 1)
    h_list = zeros(max(TEND+1,T), length(h0_list)); errorh_list = zeros(max(TEND+1,T), 1);
    l_list = zeros(max(TEND+1,T), 3);
    for k = 1:T
        xy_k = xyXYsaved(1:2*N,k); XY_k = xyXYsaved(2*N+1:4*N,k);
        h_list(k,:) = energy_bounded_MZS(0, xy_k, XY_k, eps, parvec, 1);
        if (g == 0) && (q == 0)
            l1 = sum(kappas .* (xy_k(1:N).^2 + xy_k(N+1:2*N).^2));
            l2 = sum(kappas .* 2.*eps .* (XY_k(1:N).^2 + XY_k(N+1:2*N).^2 + 2));
            l_list(k,:) = [l1 + l2, l1, l2];
        end
    end
    errorh_list(2:T) = abs(diff(h_list(1:T,end)));
    l0_list = l_list(1,:);
end

%%%% initial conditions %%%%%%%%%%%%%%%%

% the state this run starts from, and its energy (and angular momentum)
xyvor = xyXYvor(1:2*N); XYvor = xyXYvor(2*N+1:4*N);
h0 = energy_bounded_MZS(0, xyvor, XYvor, eps, parvec), hprev = h0;

if (g == 0) && (q == 0)
    modx2 = xyvor(1:N).^2 + xyvor(N+1:2*N).^2; modX2 = XYvor(1:N).^2 + XYvor(N+1:2*N).^2;
    l0 = sum(kappas .* (modx2 + 2.*eps .* (modX2 + 2))), lprev = l0;
end

%%%% main loop %%%%%%%%%%%%%%%%%%%%%%%%%

% step TT takes the state in column TT of the file to column TT+1
for TT = TSTART:TEND

    % stop once any vortex has reached aspect ratio 3 (|X|^2 = 4/3)
    if max(xyXYvor(2*N+1:3*N).^2 + xyXYvor(3*N+1:4*N).^2) >= 4/3
        break
    end

    % integrate one step of length dt; keep only the final state
    [tout,xyXYnew] = ...
        ode113(@(t,x) uvUV_bounded_MZS(t, x, eps, parvec), [0,dt], xyXYvor, options);
    xyXYvor = xyXYnew(size(xyXYnew, 1), :)';
    xynew = xyXYvor(1:2*N); XYnew = xyXYvor(2*N+1:4*N);

    % conserved quantities, and how much they changed over this step (printed)
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

    % append the new state and energy to the file
    ncwrite(outname,'xyXYvar', xyXYvor, [1,TT+1]);
    ncwrite(outname,'evals', hnew, TT+1);
end

%%%% read back the whole run %%%%%%%%%%%%

xyXYvar = ncread(outname, 'xyXYvar');
evals = ncread(outname,'evals');

TEND = length(evals) - 1;                   % number of steps actually taken

xyXYvar = xyXYvar(:, 1:TEND+1); evals = evals(1:TEND+1);

tnow = 0:dt:dt*TEND;

% angular momentum errors: total drift, and sum of the step-to-step changes
if g + q == 0
    lvals = ncread(outname,'lvals'); lvals = lvals(1:TEND+1);
    total_l_error = abs(lvals(TEND+1) - lvals(1))
    cumulative_l_error = sum(abs(lvals(2:TEND+1) - lvals(1:TEND)))
end

% largest energy change from the start of this run
max(abs(evals-h0))

%%%% plots (switch each on or off with 1==1 / 1==0) %%%%

% Plot of energy change in each term of H
% grouped as self-interaction h1, point vortex h2+h3+h6, elliptic h4+h5+h7, total H
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

% Plot ellipses: movie of the vortices in the physical domain, saved as Vortices.mp4
% (red: positive circulation, blue: negative)
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
            % starting point (x0, y0) and current position (x_0, y_0), mapped to the physical plane
            x0 = xyXYvar(n,1); y0 = xyXYvar(N+n,1);
            eps_0 = eps(n); x_0 = xyXYvar(n,k); y_0 = xyXYvar(N+n,k);
            if g + q > 0
                Z0 = x0 + 1i.*y0; z0 = heart_map(Z0, domain_var); x0 = real(z0); y0 = imag(z0);
                Z_0 = x_0 + 1i.*y_0; z_0 = heart_map(Z_0, domain_var); x_0 = real(z_0); y_0 = imag(z_0);
            end
            scatter(x0, y0, col, "filled")

            x_list(n, k) = x_0; y_list(n, k) = y_0;
            plot(x_list(n, 1:k), y_list(n, 1:k), col)
            % aspect ratio and orientation from (X, Y), then the ellipse semi-axes a, b
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