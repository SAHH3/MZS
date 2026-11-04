tic
%
% script for doubly-periodic MZS code
%

%%%%%%%%%%% Main parameters %%%%%%%%%%%%%%

%rng('shuffle')      %random number shuffle
rng(3)
N = 10;              %number of vortices
dt = 0.05;           %time interval for integration steps
delta1 = 10^(-9);    %tolerances for ODE113, 1 RelTol, 2 AbsTol
delta2 = 10^(-9);
%delta2 = vertcat(ones(2*N,1)*10^(-9), ones(2*N,1)*10^(-7));
NI = 6;              %number of terms in Weiss-McWilliams summation 
TEND = 800;        % run end time

options = odeset('RelTol',delta1, 'AbsTol',delta2, 'Events',@ElongationEvent);   %reset events

kappas = [ones(1,N/2), -ones(1,N/2)]';         %neutral vortex gas
%kappas = 1;
parvec = [NI; kappas];                    %parameters for integrator

%%%% set up netcdf output %%%%%%%%%%%%%%%%

flname=' - file.nc';
outname=strcat(pwd,flname);

% restart if file exists, create if not

if exist(outname,'file')==0

    nccreate(outname,'xyXYvar','Dimensions',{'xypos' 4*N 'list' Inf}, ...
              'Format','classic');
    nccreate(outname,'evals','Dimensions',{'list' Inf}, ...
              'Format','classic');
    nccreate(outname,'tnow','Dimensions',{'list' Inf}, ...
              'Format','classic');

    %phi_0 = rand(N,1)*2*pi;
    phi_0 = zeros(N,1);
    %lambda_0 = rand(N,1)*3;
    lambda_0 = ones(N,1)*1;
    X = (lambda_0.^(1/2) - lambda_0.^(-1/2)) .* cos(2*phi_0);
    Y = (lambda_0.^(1/2) - lambda_0.^(-1/2)) .* sin(2*phi_0);
    XYvor0 = vertcat(X,Y);
    %eps = rand(N,1);
    eps = ones(N,1)*0.007;

    e0=0; eshell=10;
    [ xyvor0 ] = set_ic_MZS_chain_punc(N, e0, eshell, XYvor0, eps, 3);
    %xyvor0 = [pi; pi];
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

    %h0 = energy_doubly_MZS(0, xyvor0, XYvor0, eps, parvec);
    h0_list = energy_doubly_MZS(0, xyvor0, XYvor0, eps, parvec, 1)
    h0 = h0_list(5);                               % initial energy
    errorh_list = zeros(TEND+1,1); h_list = zeros(TEND+1,5);
    h_list(1,:) = h0_list;

    tnow=0; dt1=0; active=ones(1,N); inds=(1:N);
    xyvor=xyfull([active==1,active==1]);
    XYvor = XYvor0;
    xyXYvor = vertcat(xyvor, XYvor);

    TSTART=1;
    ncwrite(outname,'tnow', tnow, 1);
    ncwrite(outname,'xyXYvar', xyXYvor, [1,1]);
    ncwrite(outname,'evals', h0, 1);

else
    vinfo = ncinfo(outname,'xyXYvar');
    xyXYvardim=vinfo.Size; T=xyXYvardim(2); TSTART=T+1;
    xyXYvor=ncread(outname,'xyXYvar',[1,T],[4*N,1],[1,1]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% set up initial conditions (fixing E) %%%%%

xyvor = xyXYvor(1:2*N); XYvor = xyXYvor(2*N+1:4*N);
h0 = energy_doubly_MZS(0, xyvor, XYvor, eps, parvec), hprev = h0;

%%%%%%%%% Main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%

for TT=TSTART:TEND

    if max(xyXYvor(2*N+1:3*N).^2 + xyXYvor(3*N+1:4*N).^2) >= 4/3
        break
    end

    [tout,xyXYnew] = ode113(@(t,x) ...
         uvUV_doubly_MZS(t,x,eps,parvec),[dt1,dt],xyXYvor,options);

    tout(length(tout));
    xyXYnew = xyXYnew(size(xyXYnew, 1), :)';
    xynew = xyXYnew(1:2*N); XYnew = xyXYnew(2*N+1:4*N);
    xynew = mod(xynew, 2*pi);

    %hnew = energy_doubly_MZS(0, xynew, XYnew, eps, parvec);
    hnew_list = energy_doubly_MZS(0, xynew, XYnew, eps, parvec, 1);
    h_list(TT+1,:) = hnew_list;
    hnew = hnew_list(5);
    errorh=abs(hnew-hprev)
    errorh_list(TT+1) = errorh;
    tnow = tnow+dt
    hprev=hnew;
    xyXYvor = vertcat(xynew, XYnew);

    ncwrite(outname,'tnow', tnow, TT+1);
    ncwrite(outname,'xyXYvar', xyXYvor, [1,TT+1]);
    ncwrite(outname,'evals', hnew, TT+1);

end

%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%

%herror_list

tnow = ncread(outname,'tnow');
xyXYvar = ncread(outname,'xyXYvar');
evals = ncread(outname,'evals');

TEND = length(tnow) - 1;

total_h_error = abs(evals(TEND+1) - evals(1))

h0_list
max(abs(evals-h0))

%hit_time
%XY_diff = abs(xyXYvar(2*N+1:4*N,end) - XY_end)'

% Plot of locations
if 1==0
    %subplot(2,3,1)
    hold on
    for k=1:1:N
        plot(xyXYvar(k,:), xyXYvar(N+k,:))
        scatter(xyvor0(k),xyvor0(N+k))
    end
    hold off
    xlim([0,2*pi])
    ylim([0,2*pi])
    xlabel('x')
    ylabel('y')
    title('Vortex Positions')
end

% Plot of x and y against t
if 1==0
    %subplot(2,3,2)
    hold on
    for k=1:1:N
        plot(tnow, xyXYvar(k,:), LineStyle="-")
        plot(tnow, xyXYvar(N+k,:), LineStyle="--")
    end
    hold off
    xlabel('t')
    title('X and Y')
end

% Plot of X and Y against t
if 1==0
    %subplot(2,3,2)
    hold on
    for k=1:1:N
        plot(tnow, xyXYvar(2*N+k,:), LineStyle="-")
        plot(tnow, xyXYvar(3*N+k,:), LineStyle="--")
    end
    hold off
    xlabel('t')
    title('X and Y')
end

% Plot of |X| against t
if 1==0
    %subplot(2,3,3)
    hold on
    for k=1:1:N
        plot(tnow, sqrt(xyXYvar(2*N+k,:).^2 + xyXYvar(3*N+k,:).^2))
    end
    plot([0, TEND*dt], [2/sqrt(3), 2/sqrt(3)], LineStyle="--")
    hold off
    xlabel('t')
    title('|X|')
end

% Histogram of |X|
if 1==0
    for k=0:1:8
        subplot(3,3,k+1)
        histogram(sqrt(xyXYvar(2*N+1:3*N, ceil(k*TEND/8)+1).^2 + xyXYvar(3*N+1:4*N, ceil(k*TEND/8)+1).^2), 100)
        xlabel('|X|')
        ylabel("Number of Vortices")
        title(sprintf("Histogram of |X| at t = %d", ceil(k*TEND/8)*dt))
    end
end

% Plot of H against t
if 1==0
    subplot(1,2,1)
    %subplot(2,3,4)
    plot(tnow, evals)
    xlabel('t')
    ylabel('H')
    title('H against t')
end

% Plot of error in H against t
if 1==0
    subplot(1,2,2)
    %subplot(2,3,5)
    plot(tnow, errorh_list(1:TEND+1))
    xlabel('t')
    ylabel('Error in H')
    title('Error in H against t')
end

% Plot how much of the energy is in each term of H
if 1==0
    %subplot(2,3,6)
    subplot(1,2,1)
    hold on
    plot(tnow, h_list(1:TEND+1, 1), LineWidth=1)
    plot(tnow, h_list(1:TEND+1, 3), LineWidth=1)
    xlabel('t')
    ylabel('Energy')
    title("Energy in each term of H")
    legend("h1", "h3")

    subplot(1,2,2)
    hold on
    plot(tnow, h_list(1:TEND+1, 2), LineWidth=1)
    plot(tnow, h_list(1:TEND+1, 5), LineWidth=1)
    xlabel('t')
    ylabel('Energy')
    title("Energy in each term of H")
    legend("h2", "H")
end

% Plot of energy change in each term of H
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

% Plot ellipses
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

% Plot ellipses at the end
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

toc