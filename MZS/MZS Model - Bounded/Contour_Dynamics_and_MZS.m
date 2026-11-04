%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% script for bounded MZS code
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% Main parameters %%%%%%%%%%%%%%

%rng('shuffle')                             %random number shuffle
rng(1)
addpath('../Conformal Gap/');              %path to select gap
dt = 1;                                  %time interval length
TEND = 800;                                %interval number  
N = 1; NP = 1; NM = N-NP;                  %vortex number, no. +ve, no.-ve
eps = 0.00125;
kappas = 16 * (pi^2) * eps;                %vortex circulations
domain_var = 1;                            %half the size of the gap
parvec = [kappas; domain_var];             %parameters for ODEs
nvor=1e6;                                  %i.c. parameter
delta2=3e-14;                              %position tolerance 

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
    %eps = ones(N,1)*0.00125;

    x0 = -4; y0 = 0.997; z0 = x0 + 1i*y0;
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
    xyXYvardim=vinfo.Size; T=xyXYvardim(2); TSTART=T+1;
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
    hnew_list = energy_bounded_MZS(0, xynew, XYnew, eps, parvec, 1)
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
T = 0:0.01:2*pi;
x_list = zeros(1,TEND); y_list = zeros(1,TEND);



%   Plot cd movie and also save as matlab file
%   July 2011: no printing - much faster
%              sensible defaults
%              calculation of area, centroid, orientation and aspect ratio
    
datafile = input('CD data file (default = fort.11)  ','s');
if isempty(datafile)
      datafile = 'fort.11'; 
end
numhead = 3;  % ignore these header lines
cdata  = readmatrix ('fort.11', 'FileType', 'text', 'NumHeaderLines', numhead);

    
fprintf('Maximum number of frames is some submultiple of %u \n', cdata(1,1));
nframes = input  ('How many frames to plot (default = all)?  ');
if isempty(nframes)
      nframes = 1000; 
end

xmin = input('what xmin (default = -6)?  ');
if isempty(xmin)
      xmin = -6; 
end

xmax = input('what xmax (default = 2)?  ');
if isempty(xmax)
      xmax = 2; 
end

ymin = input('what ymin (default = -2)?  ');
if isempty(ymin)
      ymin = -2; 
end

ymax = input('what ymax (default = 2)?  ');
if isempty(ymax)
      ymax = 2; 
end

fdirout = 'figfiles';
if exist(fdirout, 'dir')
      disp('directory figfiles already exists');
      reply = input('Overwrite files it? y/n [n]: ', 's');
      if isempty(reply) reply = 'n'; end;
      if reply ~= 'y', return, end;
else
      mkdir figfiles 
end

    disp ('You have the option to save the movie as a matlab movie file');
    moviefile = input('Save CD film data to (default = do not save)  ','s');
    savemovie = true;
    if isempty(moviefile)
          savemovie = false; 
    end
    
%   Start drawing frames one-by-one
    row = 1;
    aarea = zeros(1, nframes);
    centre = zeros(1, nframes);
    phi = zeros(1, nframes);
    lam = zeros(1, nframes);
    frame_root = 'CD_';

    frames = 0;
    %abc = '(c)';
    %frames = [1, 1 + 15/dt, 1 + 25/dt, 1 + 40/dt, 1 + 105/dt, nframes];
    clear cdmovie
    for k = 1:nframes
%	     Get the number, n, of contours in this frame
         n   = cdata(row,1);
         npt = cdata(row,2);
         t   = cdata(row,3);
         time(k) = t;
%	     Get the lengths & starting points of contours         
         np(1:n) = cdata(row+1:row+n,1);
         mb(1:n) = cdata(row+1:row+n,2);
         
         row = row + 1;
         
% 	     Get the position of the nodes for all contours in this frame
	     zall = complex( cdata(row+n+1:row+n+npt,1), cdata(row+n+1:row+n+npt,2));
	     
%	     Now plot individual loops one-by-one
         fig1 = figure (1);
         for j = 1,1 ; 
%        Only need first loop plotted in this problem (change to n for all loops)

%        First calculate various moments of the patch
              z = zall([ mb(j)+1:mb(j)+np(j) mb(j)+1]);
              nz = length(z);
              dz = z(2:nz) - z(1:nz-1);
              zm  = (z(2:nz) + z(1:nz-1))/2;
              zmc_dz = conj(zm).*dz;
              
              area = 0.5*imag(sum(zmc_dz));
              aarea(k)=area;
              
              zc = sum(zm.*zmc_dz)/(2i*area);
              
              zz = sum( (zm-zc).^2 .* zmc_dz)/(2i);
              zmod = sum( (zm-zc).*conj(zm-zc).^2 .* dz)/(4i);
              
              ixx = 0.5*real(zmod+zz);
              iyy = 0.5*real(zmod-zz);
              ixy = 0.5*imag(zz);
              cov = [ ixx ixy; ixy iyy];
              [V,D] = eig(cov);              
              smaxes = 2.*sqrt(diag(D)/area);
              smajor = max(smaxes);
              sminor = min(smaxes);             
              phi(k) = 0.5*angle(zz);
              lam(k) = smajor/sminor;           
              
              centre(k) = zc;
              aa =      smajor*exp(1i*phi(k));
              bb = 1i * sminor*exp(1i*phi(k));

              if isscalar(frames)
                  plot (z, Color='blue')
                  hold on
                  plot (zc,'+', Color='red')
                  plot (centre(1:k),'--', Color='cyan')
                  plot([zc-aa zc+aa], Color='magenta')
                  plot([zc-bb zc+bb], Color='green')
              else
                  hold on
                  if k == frames(1)
                      plot (z, Color='blue')
                      plot (zc,'+', Color='red')
                      plot (centre(1:k),'--', Color='cyan')
                      plot([zc-aa zc+aa], Color='magenta')
                      plot([zc-bb zc+bb], Color='green')
                  elseif k == frames(2)
                      plot (z, Color='blue')
                      plot (zc,'+', Color='red')
                      plot (centre(1:k),'--', Color='cyan')
                      plot([zc-aa zc+aa], Color='magenta')
                      plot([zc-bb zc+bb], Color='green')
                  elseif k == frames(3)
                      plot (z, Color='blue')
                      plot (zc,'+', Color='red')
                      plot (centre(1:k),'--', Color='cyan')
                      plot([zc-aa zc+aa], Color='magenta')
                      plot([zc-bb zc+bb], Color='green')
                  elseif k == frames(4)
                      plot (z, Color='blue')
                      plot (zc,'+', Color='red')
                      plot (centre(1:k),'--', Color='cyan')
                      plot([zc-aa zc+aa], Color='magenta')
                      plot([zc-bb zc+bb], Color='green')
                  elseif k == frames(5)
                      plot (z, Color='blue')
                      plot (zc,'+', Color='red')
                      plot (centre(1:k),'--', Color='cyan')
                      plot([zc-aa zc+aa], Color='magenta')
                      plot([zc-bb zc+bb], Color='green')
                      time_now = (frames(5) - 1)*dt;
                      txt = ['t = ' num2str(time_now)];
                      text(real(centre(k)), -1.5, txt, HorizontalAlignment="center", FontSize=15)
                  elseif k == frames(6)
                      plot (z, Color='blue')
                      plot (zc,'+', Color='red')
                      plot (centre(1:k),'--', Color='cyan')
                      plot([zc-aa zc+aa], Color='magenta')
                      plot([zc-bb zc+bb], Color='green')
                      time_now = (frames(6) - 1)*dt;
                      txt = ['t = ' num2str(time_now)];
                      text(real(centre(k)), -1.5, txt, HorizontalAlignment="center", FontSize=15)
%                  elseif k == frames(7)
 %                     plot (z, Color='blue')
  %                    plot (zc,'+', Color='red')
   %                   plot (centre(1:k),'--', Color='cyan')
    %                  plot([zc-aa zc+aa], Color='magenta')
     %                 plot([zc-bb zc+bb], Color='green')
                  end
              end
         end

         % Plot MZS
         %clf
         %plot([-10, -domain_var], [0, 0], Color='black')
         %plot([domain_var, 10], [0, 0], Color='black')
         if kappas > 0
             col = 'red';
         else
             col = 'blue';
         end
         eps_0 = eps; x_0 = xyXYvar(1,k); y_0 = xyXYvar(2,k);
         Z_0 = x_0 + 1i.*y_0; z_0 = gap_map(Z_0, domain_var);
         x_0 = real(z_0); y_0 = imag(z_0); X_0 = xyXYvar(3,k);
         Y_0 = xyXYvar(4,k); absX2 = X_0^2 + Y_0^2;
         x_list(k) = x_0; y_list(k) = y_0;
         lambda = ((sqrt(absX2) + sqrt(absX2+4))^2)/4; Phi = atan2(Y_0,X_0)/2;
         if k==1
             lambda = lambda_0; Phi = phi_0;
         end
         a=sqrt(8*eps_0*lambda);
         b=sqrt(8*eps_0/lambda);
         if isscalar(frames)
             plot(x_list(1:k), y_list(1:k), col, Linestyle=":")
             fill(a*cos(T)*cos(Phi) - b*sin(T)*sin(Phi) + x_0, a*cos(T)*sin(Phi) + b*sin(T)*cos(Phi) + y_0, col, FaceAlpha=0.3)
         else
             if k == frames(1)
                 text(-5.5, 1.8, abc, HorizontalAlignment="center", FontSize=15)
                 plot(x_list(1:k), y_list(1:k), col, Linestyle=":")
                 fill(a*cos(T)*cos(Phi) - b*sin(T)*sin(Phi) + x_0, a*cos(T)*sin(Phi) + b*sin(T)*cos(Phi) + y_0, col, FaceAlpha=0.3)
                 time_now = (frames(1) - 1)*dt;
                 txt = ['t = ' num2str(time_now)];
                 text(x_0, 1.5, txt, HorizontalAlignment="center", FontSize=15)
             elseif k == frames(2)
                 plot(x_list(1:k), y_list(1:k), col, Linestyle=":")
                 fill(a*cos(T)*cos(Phi) - b*sin(T)*sin(Phi) + x_0, a*cos(T)*sin(Phi) + b*sin(T)*cos(Phi) + y_0, col, FaceAlpha=0.3)
                 time_now = (frames(2) - 1)*dt;
                 txt = ['t = ' num2str(time_now)];
                 text(x_0, 1.5, txt, HorizontalAlignment="center", FontSize=15)
             elseif k == frames(3)
                 plot(x_list(1:k), y_list(1:k), col, Linestyle=":")
                 fill(a*cos(T)*cos(Phi) - b*sin(T)*sin(Phi) + x_0, a*cos(T)*sin(Phi) + b*sin(T)*cos(Phi) + y_0, col, FaceAlpha=0.3)
                 time_now = (frames(3) - 1)*dt;
                 txt = ['t = ' num2str(time_now)];
                 text(x_0, 1.3, txt, HorizontalAlignment="center", FontSize=15)
             elseif k == frames(4)
                 plot(x_list(1:k), y_list(1:k), col, Linestyle=":")
                 fill(a*cos(T)*cos(Phi) - b*sin(T)*sin(Phi) + x_0, a*cos(T)*sin(Phi) + b*sin(T)*cos(Phi) + y_0, col, FaceAlpha=0.3)
                 time_now = (frames(4) - 1)*dt;
                 txt = ['t = ' num2str(time_now)];
                 text(0.5, -0.5, txt, HorizontalAlignment="center", FontSize=15)
             elseif k == frames(5)
                 plot(x_list(1:k), y_list(1:k), col, Linestyle=":")
                 fill(a*cos(T)*cos(Phi) - b*sin(T)*sin(Phi) + x_0, a*cos(T)*sin(Phi) + b*sin(T)*cos(Phi) + y_0, col, FaceAlpha=0.3)
                 time_now = (frames(5) - 1)*dt;
                 txt = ['t = ' num2str(time_now)];
                 text(x_0, 1.3, txt, HorizontalAlignment="center", FontSize=15)
             elseif k == frames(6)
                 plot(x_list(1:k), y_list(1:k), col, Linestyle=":")
                 fill(a*cos(T)*cos(Phi) - b*sin(T)*sin(Phi) + x_0, a*cos(T)*sin(Phi) + b*sin(T)*cos(Phi) + y_0, col, FaceAlpha=0.3)
                 time_now = (frames(6) - 1)*dt;
                 txt = ['t = ' num2str(time_now)];
                 text(x_0, 1.5, txt, HorizontalAlignment="center", FontSize=15)
%             elseif k == frames(7)
 %                plot(x_list(1:k), y_list(1:k), col, Linestyle=":")
  %               fill(a*cos(T)*cos(Phi) - b*sin(T)*sin(Phi) + x_0, a*cos(T)*sin(Phi) + b*sin(T)*cos(Phi) + y_0, col, FaceAlpha=0.3)
   %              time_now = (frames(7) - 1)*dt;
    %             txt = ['t = ' num2str(time_now)];
     %            text(x_0, -1.5, txt, HorizontalAlignment="center", FontSize=15)
             end
         end

         % Finished this frame
         if isscalar(frames)
             %xlim([-5, 5])
             %ylim([-2, 2])
             %title('Evolution of Vortex')
             %daspect([1,1,1])
             xlabel('x')
             ylabel('y')
             axis equal
             axis([xmin, xmax, ymin, ymax])
             plot([xmin -1],[0 0],'k','LineWidth', 1)
             plot([ 1 xmax],[0 0],'k','LineWidth', 1)
             hold off
             cdmovie(k) = getframe;
             filename = [fdirout, '/', frame_root, num2str(k,'%03i'), '.png'];
             print (fig1, filename, '-dpng');
         else
             if ismember(k, frames)
                 %xlim([-5, 5])
                 %ylim([-2, 2])
                 %title('Evolution of Vortex')
                 %daspect([1,1,1])
                 xlabel('x')
                 ylabel('y')
                 axis equal
                 axis([xmin, xmax, ymin, ymax])
                 plot([xmin -1],[0 0],'k','LineWidth', 1)
                 plot([ 1 xmax],[0 0],'k','LineWidth', 1)
                 hold off
                 cdmovie(k) = getframe;
                 filename = [fdirout, '/', frame_root, num2str(k,'%03i'), '.png'];
                 print (fig1, filename, '-dpng');
             end
         end
         row = row + n + npt + 1;
         if row > length(cdata) 
              break 
         end
    end
    if savemovie 
         %save(moviefile, 'cdmovie', 'xmax', 'ymax') 
         writer = VideoWriter('CD & MZS', 'MPEG-4');
         writer.FrameRate = 20;
         open(writer);
         writeVideo(writer, cdmovie);
         close(writer);
    end
    
    %   Plot time series of moments 
    fig2 = figure (2);
    %dt = cdata(2,1); % time step in integrations
    ns = cdata(1,2); % number of time steps per save 
    dtime=diff(time);
    dphi = diff(unwrap(2*phi(1:k)))./(2.*dtime);
    dxdt = diff(real(centre(1:k)))./dtime;
    % time = (0:k-1) * dt * ns;
    hold on
    plot(time(1:k), real(centre(1:k)), Displayname='x_c', Linestyle='--', Color='black', LineWidth=1)
    plot(time(1:k), imag(centre(1:k)), Displayname='y_c', Linestyle='--', Color='red', LineWidth=1)
    %plot(time(1:k), phi(1:k), Displayname='\phi_c', Linestyle='--', Color='yellow', LineWidth=1)
    plot(time(1:k), lam(1:k), Displayname='\lambda_c', Linestyle='--', Color='blue', LineWidth=1)
    %plot(time(1:k), aarea(1:k), Displayname='area_c', Linestyle='--', Color='magenta', LineWidth=1)
    plot(0.5*(time(1:k-1)+time(2:k)), dphi, Displayname='\dot{\phi}_c', Linestyle='--', Color='green', LineWidth=1)
    %plot(0.5*(time(1:k-1)+time(2:k)), dxdt,'k')

    % Plot x, y, phi, lambda, area, dphi with MZS
    x_0 = xyXYvar(1,:); y_0 = xyXYvar(2,:); X = xyXYvar(3,:); Y = xyXYvar(4,:);
    z_0 = x_0 + 1i*y_0; z = gap_map(z_0, domain_var); x = real(z); y = imag(z);
    Phi = atan2(Y, X)/2; absX2 = X.^2 + Y.^2;
    lambda = ((sqrt(absX2) + sqrt(absX2+4)).^2)/4;
    dPhi = diff(unwrap(2*Phi)) ./ (2.*dt);
    plot(tnow, x, DisplayName='x_{MZS}', Linestyle=':', Color='black', LineWidth=1.5)
    plot(tnow, y, DisplayName='y_{MZS}', Linestyle=':', Color='red', LineWidth=1.5)
    %plot(tnow, Phi, DisplayName='\phi_{MZS}', Linestyle=':', Color='yellow', LineWidth=1.5)
    plot(tnow, lambda, DisplayName='\lambda_{MZS}', Linestyle=':', Color='blue', LineWidth=1.5)
    %plot(tnow, eps*8*pi .* ones(length(tnow),1), DisplayName='area_{MZS}', Linestyle=':', Color='magenta', LineWidth=1.5)
    plot(0.5 * (tnow(1:end-1) + tnow(2:end)), dPhi, DisplayName='\dot{\phi}_{MZS}', Linestyle=':', Color='green', LineWidth=1.5)

    %PhiDot = kappas .* ((1 ./ (eps .* sqrt(absX2 + 4))) - (X .* (absX2 + 2) ./ (y.^2 .* absX2))) ./ (8.*pi .* sqrt(absX2 + 4));
    %plot(tnow, PhiDot, DisplayName='\dot{\phi}_{test}', Linestyle='-', Color='blue', LineWidth=1.5)
    %xlim([0, 50])
    %ylim([1.25, 1.6])

    hold off
    set(gca, 'fontsize', 22)
    lgd = legend('Location','east', 'NumColumns',2);
    fontsize(lgd, 22, 'points')
    xlabel('Time t')
    filename = [fdirout, '/', frame_root, 'evolution.png'];
    print (fig2, filename, '-dpng');

    time(frames)
    tnow(frames)