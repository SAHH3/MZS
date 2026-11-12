%
% MAIN_DOUBLY_MZS: script for doubly-periodic MZS code
%
% Integrates N elliptical vortices (MZS moment model) in the doubly
% periodic domain [0, 2pi) x [0, 2pi).
%
% State vector (4N x 1): [x; y; X; Y]
%   (x, y)   position of each vortex centre
%   X + i Y  shape, (sqrt(lambda) - 1/sqrt(lambda)) exp(2 i phi), for an
%            ellipse of aspect ratio lambda and orientation phi, with
%            semi-axes sqrt(8 eps lambda), sqrt(8 eps / lambda) (area 8 pi eps)
% Equations of motion: uvUV_doubly_MZS.m. Energy: energy_doubly_MZS.m.
%
% The state and time after every step of length dt are appended to a NetCDF
% file (outname, saved in the folder ABOVE this one). If that file already
% exists, the run continues from its last saved state, so delete it to start
% afresh. The run stops early if any vortex reaches aspect ratio 3
% (ElongationEvent.m).
%

%%%%%%%%%%% Main parameters %%%%%%%%%%%%%%

%rng('shuffle')      %random number shuffle
rng(3)
N = 10;              %number of vortices (even: half +ve, half -ve)
dt = 0.05;           %time interval for integration steps
delta1 = 10^(-9);    %tolerances for ODE113, 1 RelTol, 2 AbsTol
delta2 = 10^(-9);
%delta2 = vertcat(ones(2*N,1)*10^(-9), ones(2*N,1)*10^(-7));
NI = 6;              %number of terms in Weiss-McWilliams summation
TEND = 800;        % number of time steps (run length TEND*dt)

options = odeset('RelTol',delta1, 'AbsTol',delta2, 'Events',@ElongationEvent);   %reset events

kappas = [ones(1,N/2), -ones(1,N/2)]';         %neutral vortex gas
%kappas = 1;
parvec = [NI; kappas];                    %parameters for integrator
dt1 = 0;                                  %start time of each ODE113 call

% initial ellipse shapes and sizes (set here so they also exist on a restart)
% phi_0: orientations, lambda_0: aspect ratios, eps: sizes
%phi_0 = rand(N,1)*2*pi;
phi_0 = zeros(N,1);
%lambda_0 = rand(N,1)*3;
lambda_0 = ones(N,1)*1;
X = (lambda_0.^(1/2) - lambda_0.^(-1/2)) .* cos(2*phi_0);
Y = (lambda_0.^(1/2) - lambda_0.^(-1/2)) .* sin(2*phi_0);
XYvor0 = vertcat(X,Y);
%eps = rand(N,1);
eps = ones(N,1)*0.007;

%%%% set up netcdf output %%%%%%%%%%%%%%%%

% variables in the file: xyXYvar (state after each step, 4N x steps),
% evals (energy after each step) and tnow (time of each step)
flname=' - file.nc';
outname=strcat(pwd,flname);              %i.e. '<this folder's path> - file.nc'

% restart if file exists, create if not

if exist(outname,'file')==0

    nccreate(outname,'xyXYvar','Dimensions',{'xypos' 4*N 'list' Inf}, ...
              'Format','classic');
    nccreate(outname,'evals','Dimensions',{'list' Inf}, ...
              'Format','classic');
    nccreate(outname,'tnow','Dimensions',{'list' Inf}, ...
              'Format','classic');

    % random initial positions with energy in [e0 - eshell/2, e0 + eshell/2],
    % no two vortices closer than 3 vortex radii
    e0=0; eshell=10;
    [ xyvor0 ] = set_ic_MZS_chain_punc(N, e0, eshell, XYvor0, eps, 3);
    %xyvor0 = [pi; pi];
    % (switched off) alternative: take the initial positions of the
    % fastest (or slowest) run to elongate from a saved batch of runs
    if 1==0
        M = 200; exc = 3;
        outname0 = sprintf('Doubly_Periodic - eps=0.0005, exclusion=%d, N=%d, %d runs - file.nc', exc, N, M);
        hit_times = ncread(outname0, 'hit_times');
        xy = ncread(outname0, 'xyvor0'); XY = ncread(outname0, 'XY_end');
        if 1==1
            hit_time = min(hit_times);
            min_index = find(hit_times == min(hit_times)); index = min_index(1);
        else
            hit_time = max(hit_times);
            max_index = find(hit_times == max(hit_times)); index = max_index(1);
        end
        xyvor0 = xy(:,index); XY_end = XY(:,index);
    end
    xyfull=xyvor0;

    % energy, split into its terms [h1, h2 + h5, h3, h4, H] (see energy_doubly_MZS.m)
    %h0 = energy_doubly_MZS(0, xyvor0, XYvor0, eps, parvec);
    h0_list = energy_doubly_MZS(0, xyvor0, XYvor0, eps, parvec, 1)
    h0 = h0_list(5);                               % initial energy
    errorh_list = zeros(TEND+1,1); h_list = zeros(TEND+1,5);
    h_list(1,:) = h0_list;

    % all vortices are active (kept for selecting a subset of xyfull)
    tnow=0; active=ones(1,N); inds=(1:N);
    xyvor=xyfull([active==1,active==1]);
    XYvor = XYvor0;
    xyXYvor = vertcat(xyvor, XYvor);

    % save the initial state as column 1
    TSTART=1;
    ncwrite(outname,'tnow', tnow, 1);
    ncwrite(outname,'xyXYvar', xyXYvor, [1,1]);
    ncwrite(outname,'evals', h0, 1);

else
    % file exists: continue from its last saved state (column T)
    vinfo = ncinfo(outname,'xyXYvar');
    xyXYvardim=vinfo.Size; T=xyXYvardim(2); TSTART=T;
    xyXYvor=ncread(outname,'xyXYvar',[1,T],[4*N,1],[1,1]);
    tnow = ncread(outname,'tnow',T,1);        %time of the last saved state

    % rebuild the energy diagnostics for the steps already saved
    xyXYsaved = ncread(outname,'xyXYvar',[1,1],[4*N,T]);
    xyvor0 = xyXYsaved(1:2*N,1);
    h0_list = energy_doubly_MZS(0, xyvor0, xyXYsaved(2*N+1:4*N,1), eps, parvec, 1)
    h_list = zeros(max(TEND+1,T), 5); errorh_list = zeros(max(TEND+1,T), 1);
    for k = 1:T
        h_list(k,:) = energy_doubly_MZS(0, xyXYsaved(1:2*N,k), xyXYsaved(2*N+1:4*N,k), eps, parvec, 1);
    end
    errorh_list(2:T) = abs(diff(h_list(1:T,5)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% set up initial conditions (fixing E) %%%%%

% the state this run starts from, and its energy (printed)
xyvor = xyXYvor(1:2*N); XYvor = xyXYvor(2*N+1:4*N);
h0 = energy_doubly_MZS(0, xyvor, XYvor, eps, parvec), hprev = h0;

%%%%%%%%% Main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%

% step TT takes the state in column TT of the file to column TT+1
for TT=TSTART:TEND

    % stop once any vortex has reached aspect ratio 3 (|X|^2 = 4/3)
    if max(xyXYvor(2*N+1:3*N).^2 + xyXYvor(3*N+1:4*N).^2) >= 4/3
        break
    end

    % integrate one step of length dt; keep only the final state
    [tout,xyXYnew] = ode113(@(t,x) ...
         uvUV_doubly_MZS(t,x,eps,parvec),[dt1,dt],xyXYvor,options);

    tout(length(tout));
    xyXYnew = xyXYnew(size(xyXYnew, 1), :)';
    xynew = xyXYnew(1:2*N); XYnew = xyXYnew(2*N+1:4*N);
    xynew = mod(xynew, 2*pi);                  % wrap positions back into [0, 2pi)

    % energy, and how much it changed over this step (printed)
    %hnew = energy_doubly_MZS(0, xynew, XYnew, eps, parvec);
    hnew_list = energy_doubly_MZS(0, xynew, XYnew, eps, parvec, 1);
    h_list(TT+1,:) = hnew_list;
    hnew = hnew_list(5);
    errorh=abs(hnew-hprev)
    errorh_list(TT+1) = errorh;
    tnow = tnow+dt
    hprev=hnew;
    xyXYvor = vertcat(xynew, XYnew);

    % append the new time, state and energy to the file
    ncwrite(outname,'tnow', tnow, TT+1);
    ncwrite(outname,'xyXYvar', xyXYvor, [1,TT+1]);
    ncwrite(outname,'evals', hnew, TT+1);

end

%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%

% read back the whole run
tnow = ncread(outname,'tnow');
xyXYvar = ncread(outname,'xyXYvar');
evals = ncread(outname,'evals');

TEND = length(tnow) - 1;                     % number of steps actually taken

% largest energy change from the start of this run
max(abs(evals-h0))

% (switch each plot on or off with 1==1 / 1==0)

% Plot of energy change in each term of H
% columns of h_list: self-interaction h1, point vortex h2 + h5, elliptic h3, total H
% (h4 is left out: it is zero for a neutral gas)
if 1==1
    hold on
    plot(tnow, h_list(1:TEND+1, 1) - h0_list(1), LineWidth=1)
    plot(tnow, h_list(1:TEND+1, 2) - h0_list(2), LineWidth=1)
    plot(tnow, h_list(1:TEND+1, 3) - h0_list(3), LineWidth=1)
    %plot(tnow, h_list(1:TEND+1, 4) - h0_list(4), LineWidth=1)
    plot(tnow, h_list(1:TEND+1, 5) - h0_list(5), LineWidth=1)
    hold off
    set(gca, 'fontsize', 22)
    xlabel('Time t')
    ylabel('Energy')
    %title("Energy change in each term of H")
    %lgd = legend("\Delta H_1: Self-Interaction", "\Delta H_2: Point Vortex", "\Delta H_3: Elliptic", "\Delta H_4: ", "\Delta H: Total");
    lgd = legend("\Delta H_1: Self-Interaction", "\Delta H_2: Point Vortex", "\Delta H_3: Elliptic", "\Delta H: Total");
    fontsize(lgd, 22, 'points')
end

% Plot ellipses: movie of the vortices, saved as Vortices.mp4
% (red: positive circulation, blue: negative; tracks drawn in the 3 x 3
% neighbouring periodic copies so they stay visible when they wrap around)
if 1==0
    t = 0:0.01:2*pi;
    clear movie
    for k=1:1:TEND+1
        clf
        hold on
        for n=1:1:N
            if kappas(n) > 0
                col = 'red';
            else
                col = 'blue';
            end
            M=1;
            for i=-M:1:M
                for j=-M:1:M
                    plot(unwrap(xyXYvar(n,1:k)) + i*2*pi, unwrap(xyXYvar(N+n,1:k)) + j*2*pi, col)
                end
            end
            x0 = xyXYvar(n,1); y0 = xyXYvar(N+n,1);
            scatter(x0, y0, col, "filled")

            eps_0 = eps(n); x_0 = xyXYvar(n,k); y_0 = xyXYvar(N+n,k); X_0 = xyXYvar(2*N+n,k); Y_0 = xyXYvar(3*N+n,k);
            absX2 = X_0^2 + Y_0^2;
            lambda = ((sqrt(absX2) + sqrt(absX2+4))^2)/4;
            phi = atan2(Y_0,X_0)/2;
            if k==1
                lambda = lambda_0(n); phi = phi_0(n);
            end
            a=sqrt(8*eps_0*lambda);
            b=sqrt(8*eps_0/lambda);
            fill(a*cos(t)*cos(phi) - b*sin(t)*sin(phi) + x_0, a*cos(t)*sin(phi) + b*sin(t)*cos(phi) + y_0, col, FaceAlpha=0.5)
        end
        xlim([0,2*pi])
        ylim([0,2*pi])
        xlabel('x')
        ylabel('y')
        title('Evolution of Vortices')
        daspect([1,1,1])
        hold off
        movie(k) = getframe;
    end
    writer = VideoWriter('Vortices', 'MPEG-4');
    writer.FrameRate = 30;
    open(writer);
    writeVideo(writer, movie);
    close(writer);
end

% Plot ellipses at the end: final ellipses, full tracks and starting points,
% with the periodic domain outlined by dashed lines
if 1==1
    t = 0:0.01:2*pi;
    k=TEND+1;
    clf
    hold on
    for n=1:1:N
        if kappas(n) > 0
            col = 'red';
        else
            col = 'blue';
        end
        M=1;
        for i=-M:1:M
            for j=-M:1:M
                plot(unwrap(xyXYvar(n,1:k)) + i*2*pi, unwrap(xyXYvar(N+n,1:k)) + j*2*pi, col)
            end
        end
        x0 = xyXYvar(n,1); y0 = xyXYvar(N+n,1);
        scatter(x0, y0, col, "filled")

        eps_0 = eps(n); x_0 = xyXYvar(n,k); y_0 = xyXYvar(N+n,k); X_0 = xyXYvar(2*N+n,k); Y_0 = xyXYvar(3*N+n,k);
        absX2 = X_0^2 + Y_0^2;
        lambda = ((sqrt(absX2) + sqrt(absX2+4))^2)/4;
        phi = atan2(Y_0,X_0)/2;
        a=sqrt(8*eps_0*lambda);
        b=sqrt(8*eps_0/lambda);
        fill(a*cos(t)*cos(phi) - b*sin(t)*sin(phi) + x_0, a*cos(t)*sin(phi) + b*sin(t)*cos(phi) + y_0, col, FaceAlpha=0.5)
    end
    plot([0, 2*pi], [0, 0], LineStyle="--", Color='black')
    plot([0, 2*pi], [2*pi, 2*pi], LineStyle="--", Color='black')
    plot([0, 0], [0, 2*pi], LineStyle="--", Color='black')
    plot([2*pi, 2*pi], [0, 2*pi], LineStyle="--", Color='black')

    e0 = 0.5;
    xlim([0-e0, 2*pi+e0])
    ylim([0-e0, 2*pi+e0])
    xlabel('x')
    ylabel('y')
    daspect([1,1,1])
    hold off
end