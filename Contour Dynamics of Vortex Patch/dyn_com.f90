MODULE dyn_com    
    
    ! npm: maximum total number of nodes on all contours    
    !  nm: maximum number of contours allowed    
    INTEGER,  PARAMETER :: npm = 5000, nm = 200    
    INTEGER,  PARAMETER :: kr = KIND(.0d0)    
    REAL(kr), PARAMETER :: pi=3.1415926535897932385_kr,                &    
                           zero=0.0_kr, half=0.5_kr, one=1.0_kr,       &    
                           two=2.0_kr, three=3.0_kr, twopi=two*pi,     &    
                           barg=13.0_kr, bargx=barg+two, dr=0.0001_kr  
    complex(kr), parameter :: ci=(0._kr,1._kr)		       
                              
    INTEGER, PARAMETER  :: ntab=bargx/dr         
    ! barg is chosen so that barg*K_1(barg) can be neglected    
    ! compared to 1.  bargx is chosen for safety reasons.      
    ! dr is the increment in the table.    
    ! The error in table lookup is less than dr**2.    
    
    LOGICAL,  DIMENSION(npm)    :: cornw    
    REAL(kr), DIMENSION(npm)    :: xw,yw    
    REAL(kr), DIMENSION(nm)     :: om        
    INTEGER,  DIMENSION(nm)     :: np,mb,ind    
    REAL(kr), DIMENSION(0:ntab) :: agf, bgf    
    
    ! Allocatable work arrays    
    LOGICAL,  DIMENSION(:), ALLOCATABLE :: corn    
    REAL(kr), DIMENSION(:), ALLOCATABLE :: x,y,o,u,v,dx,dy,xi,yi,xf, &    
                                           yf,xb1,yb1,xb2,yb2,theta,xzi,eta, &
					   bg,ug,vg,xg,yg,etag,fg,ag,ugap, &
					   vgap,ugneg,vgneg
					   
    complex(kr), dimension(:),allocatable::zeta,z,ucomp	  
    complex(kr) :: sumcomp
    
    INTEGER  :: n,npt,i,j,m,jjj,ng   
    REAL(kr) :: t,dt,small,dm,densf,ell,rdi,drdi,rdsq,bdsq,  &    
                rad1,rad2,cent,xloc,yloc,flux,stream,beta 
               
 
    
END MODULE dyn_com    
    
