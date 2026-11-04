%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MAIN_GAP_MZS: script for MZS code with the 'gap' domain
%
% Integrates elliptical vortices (MZS moment model) in the plane with two
% walls along the real axis, x <= -1 and x >= 1, leaving a gap between them
% (Conformal Gap folder, gap_map.m). Set up for a single vortex starting
% above the wall at x = -4. Uses the same bounded-domain equations as
% main_bounded_MZS.m (uvUV_bounded_MZS.m, energy_bounded_MZS.m); see that
% script for the state vector [x; y; X; Y].
%
% The state after every step of length dt is appended to a NetCDF file
% (outname, saved in the folder ABOVE this one). If that file already exists,
% the run continues from its last saved state, so delete it to start afresh.
% The run stops early if any vortex reaches |X|^2 = 25 (ElongationEvent2.m).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% Main parameters %%%%%%%%%%%%%%

%rng('shuffle')                             %random number shuffle
rng(1)
addpath('../Conformal Gap/');              %path to select gap
dt = 0.5;                                  %time interval length
TEND = 300;                                %number of time steps (run length TEND*dt)
N = 1; NP = 1; NM = N-NP;                  %vortex number, no. +ve, no.-ve
kappas = [ones(1,NP), -ones(1,NM)]';       %vortex circulations
domain_var = 1;                            %half the size of the gap (the map assumes 1; only used for plotting)
parvec = [kappas; domain_var];             %parameters for ODEs
nvor=1e6;                                  %i.c. parameter (not used here)
delta2=3e-10;                              %ODE113 relative and absolute tolerance

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
eps = ones(N,1)*0.045;

options=odeset('RelTol',delta2, 'AbsTol',delta2, 'Events',@ElongationEvent2); %options for ODE solver

%%%% set up netcdf output %%%%%%%%%%%%%%%%

% variables in the file: xyXYvar (state after each step, 4N x steps) and
% evals (energy after each step)
flname = ' - Gap - file.nc';
outname=strcat(pwd,flname);                %i.e. '<this folder's path> - Gap - file.nc'

% restart if file exists, create if not

if exist(outname,'file')==0

    nccreate(outname,'xyXYvar','Dimensions',{'xypos' 4*N 'list' Inf}, ...
           'Format','classic');
    nccreate(outname,'evals','Dimensions',{'list' Inf}, ...
           'Format','classic');

    % starting position, given in the physical plane and mapped to the unit disc
    x0 = -4; y0 = 0.995; z0 = x0 + 1i*y0;
    Z0 = inverse_map(z0, domain_var); x0 = real(Z0); y0 = imag(Z0);
    xyvor = [x0; y0];
    xyXYvor = vertcat(xyvor, XYvor0);

    % energy, split into its terms [h1, ..., h7, H] (see energy_bounded_MZS.m)
    %h0 = energy_bounded_MZS(0, xyvor, XYvor0, eps, parvec);
    h0_list = energy_bounded_MZS(0, xyvor, XYvor0, eps, parvec, 1);
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

    % rebuild the energy diagnostics for the steps already saved
    xyXYsaved = ncread(outname,'xyXYvar',[1,1],[4*N,T]);
    h0_list = energy_bounded_MZS(0, xyXYsaved(1:2*N,1), xyXYsaved(2*N+1:4*N,1), eps, parvec, 1);
    h_list = zeros(max(TEND+1,T), length(h0_list)); errorh_list = zeros(max(TEND+1,T), 1);
    for k = 1:T
        h_list(k,:) = energy_bounded_MZS(0, xyXYsaved(1:2*N,k), xyXYsaved(2*N+1:4*N,k), eps, parvec, 1);
    end
    errorh_list(2:T) = abs(diff(h_list(1:T,end)));
end

%%%% initial conditions %%%%%%%%%%%%%%%%

% the state this run starts from, and its energy (printed)
xyvor = xyXYvor(1:2*N); XYvor = xyXYvor(2*N+1:4*N);
h0 = energy_bounded_MZS(0, xyvor, XYvor, eps, parvec), hprev = h0;

%%%% main loop %%%%%%%%%%%%%%%%%%%%%%%%%

% step TT takes the state in column TT of the file to column TT+1
for TT = TSTART:TEND

    % stop once any vortex has reached |X|^2 = 25 (aspect ratio about 27)
    if max(xyXYvor(2*N+1:3*N).^2 + xyXYvor(3*N+1:4*N).^2) >= 25
        break
    end

    % integrate one step of length dt; keep only the final state
    [tout,xyXYnew] = ...
        ode113(@(t,x) uvUV_bounded_MZS(t, x, eps, parvec), [0,dt], xyXYvor, options);

    xyXYvor = xyXYnew(size(xyXYnew, 1), :)';
    xynew = xyXYvor(1:2*N); XYnew = xyXYvor(2*N+1:4*N);

    % energy, and how much it changed over this step (printed)
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

% largest energy change from the start of this run
max(abs(evals-h0))

%%%% plots (switch on or off with 1==1 / 1==0) %%%%

% Plot ellipses: movie of the vortices in the physical plane, saved as Vortices.mp4
% (black lines: the walls; red: positive circulation, blue: negative)
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
            % current position mapped to the physical plane; aspect ratio and
            % orientation from (X, Y); then the ellipse semi-axes a, b
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