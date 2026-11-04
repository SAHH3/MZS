 
 
Program patch  
  
! Data for planecs -- infinitely long escarpment along y=0  
  
    IMPLICIT NONE 
 
 
    INTEGER, PARAMETER     :: kr=kind(.0d0) 
    REAL(kr), PARAMETER    :: zero=0._kr,one=1._kr,two=2.d0,half=one/two  
    REAL(kr), PARAMETER    :: four=4._kr,ten=10._kr  
    REAL(kr), PARAMETER    :: pi=3.141592653589793238_kr,twopi=two*pi  
    REAL(kr), PARAMETER    :: pif=pi/180._kr  
    REAL(kr), DIMENSION(1) :: om,rad  
    INTEGER,  DIMENSION(1) :: np,mb,ind  
    REAL(kr), DIMENSION(:), ALLOCATABLE :: x,y  
    INTEGER  :: npt,i,ns,nloops,ng  
    REAL(kr) :: t,curv,ff,rdi,dm,dt,amu,ell,pradx,prady,xploc, & 
                yploc,theta,flux,stream,beta
        
  
    PRINT *,'The initial flow consists of a vortex patch near a wall with a gap.'  
  
    PRINT *, ' '   
    PRINT *,  'Enter the reciprocal of the Rossby radius (small!!!!)'  
    READ *, rdi 

    PRINT *,  'Enter the patch (ellipse) half-dimensions in the x and y directions'  
    READ *, pradx,prady
    rad(1)=pradx 
    PRINT *,  'Enter the patch vorticity (in units of 2*pi)'  
    READ *, om(1) 
!    om(1)=om(1)/twopi 
    PRINT *,  'Enter the location of the centre of the patch '  
    READ *, xploc,yploc   
 
     
    PRINT *,  'Enter the uniform stream velocity parallel to wall '  
    READ *, stream
    PRINT *,  'Enter the  stream velocity far from wall (i.e. beta) '  
    READ *, beta  
    PRINT *,  'Enter the flux of fluid through gap '  
    READ *, flux   
 
    ! Ask for algorithm parameters:  
    PRINT *, ' '   
    PRINT *,  'Next, set up the numerical parameters.'  
    PRINT *, ' '   
    PRINT *,  ' Enter the number of Fourier modes'  
    READ *,  ng
    PRINT *,  ' Enter the spatial resolution parameter, mu:'  
    PRINT *,  '(a value <= 0.15 is advisable)'  
    READ *, amu  
    ell = one  
    dm = amu**2*ell/four  
    PRINT *,  'Surgery will be performed at the scale ',dm  
    PRINT *, ' '   
    dt = 0.1_kr  
    PRINT *,  ' The time step dt should satisfy dt <= ',dt  
    PRINT *,  '(for high accuracy, half this value should be used)'  
    PRINT *, ' '   
    PRINT *,  ' Enter dt, the number of time steps per data save '  
    PRINT *,  ' and the total number of SAVES'  
    READ *, dt,ns,nloops  
  
    ! Write data for input to cs calculation:  
    OPEN (11,FILE='fort.11',STATUS='REPLACE')  
    WRITE (11,*) nloops,ns  
    WRITE (11,*) dt,amu,ell,dm  
    WRITE (11,*) rdi,ng,flux,stream,beta  
  
    ! npt: the total number of nodes  
    npt = 0  
        curv = SQRT(one/ell**2+one/rad(1)**2)  
        ff = SQRT(curv/ell)/amu+half*curv  
        np(1) = 1+NINT (twopi*rad(1)*ff)  
        npt = npt+np(1) 
	  
 
  
    ALLOCATE (x(npt),y(npt))  
    t = zero  
    WRITE (11,'(2i8,f20.10)') 1,npt,t  
    mb(1)=0  
    ind(1) = 1  
      
    WRITE (11,'(3i8,f20.10,L8)') (np(i),mb(i),ind(i),om(i), i=1,1) 
      
        DO i = 1,np(1) 
          theta=twopi*float(i-1)/float(np(1)) 
           x(i) =pradx*cos(theta)+xploc 
           y(i) =prady*sin(theta)+yploc    
        END DO 
  
       
  
    !translate the points so that vortex is at x=0 initially  
  
    WRITE (11,'(2f20.10)') xploc,yploc  
    WRITE (11,'(2f20.10,L8)') (x(i),y(i),.FALSE.,i=1,npt)  
       
    CLOSE (11)  
  
 
END Program patch 
 
