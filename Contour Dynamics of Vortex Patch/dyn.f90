!       Planar Contour Surgery    
!       20/12/93 / DGD / Toulouse    
!        July 96 / SRC,ERJ / UCL    
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::    
    
!     This program is designed to follow the dynamics of    
! piecewise-constant vorticity distributions for long times.    
     
!     Previously, contour dynamics programs could not cope with    
! the formation of corners, the splitting of a vortex into two,    
! or the merger of vortices.  The current algorithm surmounts these    
! difficulties by introducing a cutoff distance scale.  Basically,    
! if the thinness of filaments that are produced in the calculation    
! becomes less than this cutoff, the part of the filament which is    
! too thin is eliminated.  This applies (with certain restrictions)    
! to cusps, filament "bridges" between larger regions of vorticity,    
! and merger situations.  Of course, true inviscid dynamics would    
! not allow contours to ever break; however, we expect that the    
! elimination of extremely thin filaments will be dynamically    
! insignificant and enable a further look at bizzare nl behavior.    
    
! Time integration scheme: Fourth order Runge-Kutta.    
!  Spatial discretization: Nodes are redistributed at each time    
!                          step according to local curvature but    
!                          in such a way as to leave enough    
!                          resolution in low curvature regions    
!                          to allow contours to bend.    
!     Spatial integration: Explicit formula to first order in the    
!                          contour's departure from a straight line    
!                          segment.  This departure is kept small    
!                          at all times and under all circumstances    
!                          by the re-noding according to curvature.    
                
!*********************************************************************    
! n, npt:      Number of vorticity levels, contours and nodes    
! t:            Beginning time    
! dt,nloops,ns: Time step, # of time saves, # of time steps/save    
! amu:          Dimensionless node-spacing PARAMETER    
! ell:          Large-scale length of the flow    
! rdi:          Inverse Rossby radius    
! dm:           Surgical scale (NB: dm must satisfy     
!               dm >= amu*amu*ell/4.  This is forced below.)    
! np(j):        Number of nodes on contour j.    
! mb(j)+1:      Starting node index for contour j    
! ind(j),om(j): Vorticity level and vorticity jump    
!               for contour j; only same levels can reconnect.    
! corn:         LOGICAL array telling whether a node is a corner or not.    
!               The contour must make an acute angle to form one.    
!               This is needed for surgery and node redistribution.    
!    
!x,y:          Coordinates of contour j    
!    
!NOTE    
!xtop,ytop,toppts,nptop etc have the same meaning but SOLELY for the    
!topographic curves, which advect the contours but remain undeflected    
!themselves. The cubic coefficients of the topographic alptop,bettop,    
!gamtop, along with dxtop, dytop etc, are calculated once and for all    
!at the start of a run, and are only used in uv_     
!***********************************************************************    
    
PROGRAM dyn    
    
    USE dyn_com    
    USE tabgen_    
    USE advect_    
    USE dump_    
    USE renode_    
    USE surgery_    
    USE cubic_    
    
    IMPLICIT NONE    
    REAL(kr) :: amu    
    INTEGER  :: loop,nloops,ns,istep    
    CHARACTER(LEN=1) :: ch_dum    
    
     
  
    
    ! Read INITIAL data    
      
    OPEN (11,FILE='fort.11',STATUS='OLD')    
    REWIND (11)    
    READ (11,*) nloops,ns    
    READ (11,*) dt,amu,ell,dm    
    READ (11,*) rdi,ng,flux,stream,beta   
    WRITE (*,*) nloops,ns    
    WRITE (*,*) dt,amu,ell,dm    
    WRITE (*,*) rdi,ng,flux,stream,beta    
    
    ! Position data set just before EOF record so that more can be written    
    DO    
        READ  (11,*,END=88,ERR=101) n,npt,t    
        WRITE (*,*) n,npt,t    
        READ (11,*) (np(j),mb(j),ind(j),om(j), j=1,n)    
        read (11,*) xloc,yloc    
        READ (11,*) (xw(j),yw(j),ch_dum, j=1,npt)    
    END DO    
88  BACKSPACE (11)    
    CLOSE (11)    
    
    ! Initialize constants    
    
    ! Small parameters for CUBIC, UV    
    small = EPSILON(one)    
    dm = (half*amu)**2*ell    
    ! Parameter for RENODE    
    densf = one/(amu*SQRT(ell))    
    rdsq = one/rdi**2  
 
       
    
    WRITE (*,*) 'dm',dm    
    WRITE (*,*) 'densf',densf    
    WRITE (*,*) 'small',small    
    
    ALLOCATE (corn(npt))    
    ALLOCATE(x(npt),y(npt),o(npt),u(npt),v(npt),&
             dx(npt),dy(npt),ugap(npt),vgap(npt),xzi(npt),eta(npt))    
    ALLOCATE (xi(npt),yi(npt),xf(npt),yf(npt))
    allocate (z(npt),zeta(npt),ucomp(npt))    
    allocate (bg(ng),ug(ng),vg(ng),xg(ng),yg(ng),etag(ng),fg(ng),ag(ng), &
               vgneg(ng),ugneg(ng))    
    x = xw(:npt)    
    y=  yw(:npt)   
   
    
    DO j = 1,n    
        o(mb(j)+1:mb(j)+np(j)) = om(j)    
    END DO    
            
    CALL tabgen    
    CALL renode    
    WRITE(*,*) 'Renoded contour:',npt    
    
    DO loop = 1,nloops    
      DO istep = 1,ns    
        ! Step from time t to time t + dt:    
        CALL advect    
        CALL surgery    
        CALL renode    
      END DO    
      CALL dump    
   END DO    
   
101 CLOSE (11)   
   
END PROGRAM dyn    
