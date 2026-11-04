PROGRAM DIPOLE

! Data for planecs -- two vortex patches

    IMPLICIT NONE
    INTEGER, PARAMETER     :: kr=kind(.0d0)
    REAL(kr), PARAMETER    :: zero=0._kr,one=1._kr,two=2.d0,half=one/two
    REAL(kr), PARAMETER    :: four=4._kr,ten=10._kr
    REAL(kr), PARAMETER    :: pi=3.141592653589793238_kr,twopi=two*pi
    REAL(kr), PARAMETER    :: pif=pi/180._kr
    REAL(kr), DIMENSION(2) :: om,rad,yc,asp
    INTEGER,  DIMENSION(2) :: np,mb,ind
    REAL(kr), DIMENSION(:), ALLOCATABLE :: x,y
    INTEGER  :: npt,i,ns,nloops,k
    REAL(kr) :: t,cth,sth,curv,ff,ang,xx,yy,bb,aa,dang,theta,rdi,        &
                vrat,arat,dis,dm,dt,amu,ell

    PRINT *,'The initial flow consists of two elliptical vortex patches.'
    PRINT *,'The vorticity of the lower vortex is 2*pi and its area is pi.'
    PRINT *, ' ' 
    PRINT *,' Enter the vorticity ratio of the two vortices:'
    READ *, vrat
    om(1) = one
    om(2) = vrat
    PRINT *, ' ' 
    DO
        PRINT *,  ' Enter the vortex area ratio:'
        READ *, arat
        IF (arat > zero) EXIT
        PRINT *, ' The area ratio must be positive.'
    END DO
    rad(2) = SQRT(arat)
    rad(1) = one
    PRINT *, ' ' 
    PRINT *,'The centres of the two vortices are initially on the y axis.'
    PRINT *, ' ' 
    PRINT *,  ' Enter the distance between vortex centers:'
    READ *, dis
    yc(2) = dis/(one+arat)
    yc(1) = yc(2)-dis
    PRINT *, ' ' 
    PRINT *,  ' Enter the aspect ratio (height/width) of both vortices.'
    READ *, asp(1),asp(2)
    PRINT *, ' ' 
    PRINT *,  'Enter the reciprocal of the Rossby radius'
    READ *, rdi

    ! Ask for algorithm parameters:
    PRINT *, ' ' 
    PRINT *,  'Next, set up the numerical parameters.'
    PRINT *, ' ' 
    PRINT *,  ' Enter the spatial resolution parameter, mu:'
    PRINT *,  '(a value <= 0.15 is advisable)'
    READ *, amu
    ell = one
    dm = amu**2*ell/four
    PRINT *,  'Surgery will be performed at the scale ',dm
    PRINT *, ' ' 
    dt = 0.1_kr/max(one,abs(vrat))
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
    WRITE (11,*) rdi

    ! npt: the total number of nodes
    npt = 0
    DO k = 1,2
        curv = SQRT(one/ell**2+one/rad(k)**2)
        ff = SQRT(curv/ell)/amu+half*curv
        np(k) = 1+NINT (twopi*rad(k)*ff)
        npt = npt+np(k)
    END DO

    ALLOCATE (x(npt),y(npt))
    t = zero
    WRITE (11,'(2i8,f20.10)') 2,npt,t
    mb(1) = 0
    mb(2) = np(1) 
    ind(1) = 1
    ind(2) = -1
    IF ( ABS(vrat-one) < ten*epsilon(one) ) THEN
        ind(2) = 1
        om(2) = one
    END IF
    WRITE (11,'(3i8,f20.10)') (np(i),mb(i),ind(i),om(i), i=1,2)
    theta = pif*theta
    cth = COS (theta)
    sth = SIN (theta)
    DO k = 1,2
        bb = rad(k)*SQRT(asp(k))
        aa = bb/asp(k)
        dang = twopi/float(np(k))
        DO i = 1,np(k)
            ang = dang*(i-1)
            xx = aa*COS(ang)
            yy = yc(k)+bb*SIN(ang)
            x(mb(k)+i) = xx*cth-yy*sth
            y(mb(k)+i) = yy*cth+xx*sth
           
        END DO
    END DO
    WRITE (11,'(2f20.10,L8)') (x(i),y(i),.FALSE.,i=1,npt)

END PROGRAM dipole
