%   Plot cd movie and also save as matlab file
%   July 2011: no printing - much faster
%              sensible defaults
%              calculation of area, centroid, orientation and aspect ratio

    clear all
    close all

    
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
              
              plot (z)
              hold on
              centre(k) = zc;
              plot (zc,'+')
              plot (centre(1:k),'--')
              aa =      smajor*exp(1i*phi(k));
              bb = 1i * sminor*exp(1i*phi(k));
              plot([zc-aa zc+aa])
              plot([zc-bb zc+bb])
              axis equal
              axis([xmin, xmax, ymin, ymax])
         end
         
%	     Finished this frame
         plot([xmin -1],[0 0],'k','LineWidth', 1)
         plot([ 1 xmax],[0 0],'k','LineWidth', 1)
         hold off
         cdmovie(k) = getframe;
         filename = [fdirout, '/', frame_root, num2str(k,'%03i'), '.pdf'];
         print (fig1, filename, '-dpdf');      
	     row = row + n + npt + 1;
         if row > length(cdata) 
              break 
         end
    end
    if savemovie 
         %save(moviefile, 'cdmovie', 'xmax', 'ymax') 
         writer = VideoWriter('Contour', 'MPEG-4');
         writer.FrameRate = 50;
         open(writer);
         writeVideo(writer, cdmovie);
         close(writer);
    end
    
%   Plot time series of moments 
    fig2 = figure (2);
    dt = cdata(2,1); % time step in integrations
    ns = cdata(1,2); % number of time steps per save 
    dtime=diff(time);
    dphi = diff(unwrap(2*phi(1:k)))./dtime;
    dxdt = diff(real(centre(1:k)))./dtime;
   % time = (0:k-1) * dt * ns;
    hx = plot (time(1:k), real(centre(1:k)),':k');
    hold on
    hy = plot (time(1:k), imag(centre(1:k)),'--k');
    h1 = plot (time(1:k), phi(1:k),'g');
    h2 = plot (time(1:k), lam(1:k),'r');
    h3 = plot (time(1:k), aarea(1:k),'b');
    h4 = plot (0.5*(time(1:k-1)+time(2:k)), dphi,'m');   
    h5 = plot (0.5*(time(1:k-1)+time(2:k)), dxdt,'k');   
    set(hx,'Displayname','$x_c$')
    set(hy,'Displayname','$y_c$')
    set(h1,'Displayname','$\phi$')
    set(h2,'Displayname','$\lambda$')
    set(h3,'Displayname','area')
    set(h4,'Displayname','dphi')
    legend('Location','east')
    title('Centroid evolution', 'FontWeight','bold')
    filename = [fdirout, '/', frame_root, 'evolution.pdf'];
    print (fig2, filename, '-dpdf'); 

