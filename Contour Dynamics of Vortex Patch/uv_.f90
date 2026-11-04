MODULE uv_   
   
! Calculates the velocity field at each of the nodal   
! positions on all of the contours. The velocity is calculated    
! to first order in the departure from a straight line between nodes.   
   
    USE dyn_com   
    USE cubic_   
    IMPLICIT NONE   
    PRIVATE   
    PUBLIC uv   
   
CONTAINS   
   
SUBROUTINE uv   
      
   
    REAL(kr), PARAMETER :: p1=0.21132486540519_kr, p2=0.78867513459481_kr, &   
                dp1=two*p1, dp2=two*p2, tp1=three/two*p1, tp2=three/two*p2   
    REAL(kr) :: rr1,qy1,qy2,qx1,qx2,dd1,gf1,gg,rr,ff,dd2,rr2,qy,bb,qx,gf2, &   
                aa,ee,hh,bs,cc,ac,cs,dl   
    INTEGER  :: jj,k1,k2,nn  
    REAL(kr), DIMENSION(:), ALLOCATABLE :: xx,yy,uu,uuu,vv,vvv,odx,ody,edx,edy,    &   
                                           a,b,c,d,e,f,g,h,eu,alp,bet,gam,     &   
                                           uplus,vplus,uplusim,vplusim, &
					   uneg,vneg,unegim,vnegim, &  
                                           xxx,yyy,ddx,ddy  
 					     
    INTEGER,  DIMENSION(:), ALLOCATABLE :: loc 
      
   
    ALLOCATE (odx(npt),ody(npt),edx(npt),edy(npt),a(npt),b(npt),   &   
              c(npt),d(npt),e(npt),f(npt),g(npt),h(npt),eu(npt),   &   
              alp(npt),bet(npt),gam(npt),uplus(npt),ddx(npt),ddy(npt),      &   
              vplus(npt),uplusim(npt),vplusim(npt),uneg(npt),vneg(npt), &
	      unegim(npt),vnegim(npt))
    
    CALL cubic (eu,alp,bet,gam,c,e,g)  
    odx = o*dx  
    ody = o*dy  
    edx = eu*dx  
    edy = eu*dy            
    vplus= zero   
    uplus = zero
    uplusim=zero
    vplusim=zero
    uneg=zero
    vneg=zero
    unegim=zero
    vnegim=zero
    
    ug=zero
    vg=zero
    ugneg=zero
    vgneg=zero 
    
    v = zero   
    u = zero   
    a = x + p1*dx - dy*p1*(alp+p1*(bet+p1*gam))  
    b = y + p1*dy + dx*p1*(alp+p1*(bet+p1*gam))  
    c = odx - ody*(alp+dp1*(bet+tp1*gam))  
    d = ody + odx*(alp+dp1*(bet+tp1*gam))  
    e = x + p2*dx - dy*p2*(alp+p2*(bet+p2*gam))  
    f = y + p2*dy + dx*p2*(alp+p2*(bet+p2*gam))  
    g = odx - ody*(alp+dp2*(bet+tp2*gam))  
    h = ody + odx*(alp+dp2*(bet+tp2*gam))  
    ALLOCATE (xx(npt),yy(npt),yyy(npt),xxx(npt),uu(npt),vv(npt),uuu(npt),    &  
              loc(npt),vvv(npt))
  
    ddy=dy  
    ddx=dx  
    yyy=y  
    xxx=x
    
    do i=1,npt
    if (y(i) .lt. zero)then     
    y(i)=zero
    endif  
    enddo           
 
    CALL cubic (eu,alp,bet,gam,c,e,g)   
 
    odx = o*dx   
    ody = o*dy   
    edx = eu*dx   
    edy = eu*dy  
 
    
   
    a = x + p1*dx - dy*p1*(alp+p1*(bet+p1*gam))   
    b = y + p1*dy + dx*p1*(alp+p1*(bet+p1*gam))   
    c = odx - ody*(alp+dp1*(bet+tp1*gam))   
    d = ody + odx*(alp+dp1*(bet+tp1*gam))   
    e = x + p2*dx - dy*p2*(alp+p2*(bet+p2*gam))   
    f = y + p2*dy + dx*p2*(alp+p2*(bet+p2*gam))   
    g = odx - ody*(alp+dp2*(bet+tp2*gam))   
    h = ody + odx*(alp+dp2*(bet+tp2*gam))   
 
	       
   
   
    !Find the velocity at all nodes DUE to i-th node   
u1: DO i = 1,npt   
    
        nn = 0   
   
   
        !Find only the nodes that will contribute to the integral   
        DO j = 1,npt   
            IF ((xxx(j)-x(i))**2+(yyy(j)-y(i))**2 > bdsq) CYCLE	      
            nn = nn+1   
            loc(nn) = j   
        END DO   
   
        DO jj = 1,nn   
            xx(jj) = xxx(loc(jj))   
            yy(jj) = yyy(loc(jj))   
        END DO   
   
        DO jj = 1,nn   
   
            ! Non-singular part (table look-up)   
            qx1 = xx(jj)-a(i)   
            qy1 = yy(jj)-b(i)   
            dd1 = qx1**2+qy1**2   
            rr1 = drdi*SQRT(dd1)   
            k1 = INT(rr1)   
               
            gf1 = (agf(k1)+bgf(k1)*rr1)*(qx1*c(i)+qy1*d(i))/dd1   
               
            qx2 = xx(jj)-e(i)   
            qy2 = yy(jj)-f(i)   
            dd2 = qx2**2+qy2**2   
            rr2 = drdi*SQRT(dd2)   
            k2 = INT(rr2)   
            gf2 = (agf(k2)+bgf(k2)*rr2)*(qx2*g(i)+qy2*h(i))/dd2   
	       
            ! Singular part (explicit)    
            qx = xx(jj)-x(i)   
            qy = yy(jj)-y(i)   
            bb = edx(i)*qx+edy(i)*qy   
            cc = edx(i)*qy-edy(i)*qx   
            bs = bb*bb   
            cs = cc*cc+small   
            ac = two*bb   
            dl = half*LOG(one+(one-ac)/(bs+cs))   
            aa = bs-cs   
            rr = ATAN2(cc,cs-bb*(one-bb))   
            ff = cc*(alp(i)+ac*bet(i)+(three*bs-cs)*gam(i))   
            gg = bb*(alp(i)+(bs-three*cs)*gam(i))+aa*bet(i)-cc   
            hh = one+gg*rr+ff*dl+cc*(bet(i)+(half+ac)*gam(i))   
            ee = gg*dl-ff*rr+(bb-half)*bet(i)+(half*bb-two/three+aa)*gam(i)   
            uu(jj) = odx(i)*hh-ody(i)*ee+gf1*qx1+gf2*qx2   
            vv(jj) = ody(i)*hh+odx(i)*ee+gf1*qy1+gf2*qy2   
        END DO   
        DO jj = 1,nn   
            j = loc(jj)   
            uplus(j) = uplus(j)+uu(jj)   
            vplus(j) = vplus(j)+vv(jj)   
        END DO   
    END DO u1   

 

do j=1,ng
   etag(j)=pi*(float(j)/float(ng+1))
   xg(j)=cos(etag(j))
   yg(j)=zero  
u2: DO i = 1,npt   
    
        nn = 0   
   
   
        !Find only the nodes that will contribute to the integral   
          
            IF ((xg(j)-x(i))**2+(yg(j)-y(i))**2 > bdsq) CYCLE	      
       
 
   
  
   
            ! Non-singular part (table look-up)   
            qx1 = xg(j)-a(i)   
            qy1 = yg(j)-b(i)   
            dd1 = qx1**2+qy1**2   
            rr1 = drdi*SQRT(dd1)   
            k1 = INT(rr1)   
               
            gf1 = (agf(k1)+bgf(k1)*rr1)*(qx1*c(i)+qy1*d(i))/dd1   
               
            qx2 = xg(jj)-e(i)   
            qy2 = yg(jj)-f(i)   
            dd2 = qx2**2+qy2**2   
            rr2 = drdi*SQRT(dd2)   
            k2 = INT(rr2)   
            gf2 = (agf(k2)+bgf(k2)*rr2)*(qx2*g(i)+qy2*h(i))/dd2   
	       
            ! Singular part (explicit)    
            qx = xg(j)-x(i)   
            qy = yg(j)-y(i)   
            bb = edx(i)*qx+edy(i)*qy   
            cc = edx(i)*qy-edy(i)*qx   
            bs = bb*bb   
            cs = cc*cc+small   
            ac = two*bb   
            dl = half*LOG(one+(one-ac)/(bs+cs))   
            aa = bs-cs   
            rr = ATAN2(cc,cs-bb*(one-bb))   
            ff = cc*(alp(i)+ac*bet(i)+(three*bs-cs)*gam(i))   
            gg = bb*(alp(i)+(bs-three*cs)*gam(i))+aa*bet(i)-cc   
            hh = one+gg*rr+ff*dl+cc*(bet(i)+(half+ac)*gam(i))   
            ee = gg*dl-ff*rr+(bb-half)*bet(i)+(half*bb-two/three+aa)*gam(i)   
            ug(j) = ug(j)+odx(i)*hh-ody(i)*ee+gf1*qx1+gf2*qx2   
            vg(j) = vg(j)+ody(i)*hh+odx(i)*ee+gf1*qy1+gf2*qy2   

    END DO u2
enddo 

 
     do i=1,npt 
    y(i)=yyy(i) 
    dy(i)=ddy(i)  
    x(i)=xxx(i)  
    dx(i)=ddx(i)  
    enddo 
 
 
    
    ddy=dy  
    ddx=dx  
    yyy=y  
    xxx=x 
    
    do i=1,npt
    if (y(i) .lt. zero)then     
    y(i)=zero
    endif  
    enddo       
    
    do i=1,npt
      y(i)=-y(i)
    enddo     
  
     CALL cubic (eu,alp,bet,gam,c,e,g)   
 
    odx = o*dx   
    ody = o*dy   
    edx = eu*dx   
    edy = eu*dy  
 

   
    v = zero   
    u = zero   
    a = x + p1*dx - dy*p1*(alp+p1*(bet+p1*gam))   
    b = y + p1*dy + dx*p1*(alp+p1*(bet+p1*gam))   
    c = odx - ody*(alp+dp1*(bet+tp1*gam))   
    d = ody + odx*(alp+dp1*(bet+tp1*gam))   
    e = x + p2*dx - dy*p2*(alp+p2*(bet+p2*gam))   
    f = y + p2*dy + dx*p2*(alp+p2*(bet+p2*gam))   
    g = odx - ody*(alp+dp2*(bet+tp2*gam))   
    h = ody + odx*(alp+dp2*(bet+tp2*gam))   
 
	       
   
    
    !Find the velocity at all nodes DUE to i-th node   
u3: DO i = 1,npt   
    
        nn = 0   
   
   
        !Find only the nodes that will contribute to the integral   
        DO j = 1,npt   
            IF ((xxx(j)-x(i))**2+(yyy(j)-y(i))**2 > bdsq) CYCLE	      
            nn = nn+1   
            loc(nn) = j   
        END DO   
   
        DO jj = 1,nn   
            xx(jj) = xxx(loc(jj))   
            yy(jj) = yyy(loc(jj))   
        END DO   
   
        DO jj = 1,nn   
   
            ! Non-singular part (table look-up)   
            qx1 = xx(jj)-a(i)   
            qy1 = yy(jj)-b(i)   
            dd1 = qx1**2+qy1**2   
            rr1 = drdi*SQRT(dd1)   
            k1 = INT(rr1)   
               
            gf1 = (agf(k1)+bgf(k1)*rr1)*(qx1*c(i)+qy1*d(i))/dd1   
               
            qx2 = xx(jj)-e(i)   
            qy2 = yy(jj)-f(i)   
            dd2 = qx2**2+qy2**2   
            rr2 = drdi*SQRT(dd2)   
            k2 = INT(rr2)   
            gf2 = (agf(k2)+bgf(k2)*rr2)*(qx2*g(i)+qy2*h(i))/dd2   
	       
            ! Singular part (explicit)    
            qx = xx(jj)-x(i)   
            qy = yy(jj)-y(i)   
            bb = edx(i)*qx+edy(i)*qy   
            cc = edx(i)*qy-edy(i)*qx   
            bs = bb*bb   
            cs = cc*cc+small   
            ac = two*bb   
            dl = half*LOG(one+(one-ac)/(bs+cs))   
            aa = bs-cs   
            rr = ATAN2(cc,cs-bb*(one-bb))   
            ff = cc*(alp(i)+ac*bet(i)+(three*bs-cs)*gam(i))   
            gg = bb*(alp(i)+(bs-three*cs)*gam(i))+aa*bet(i)-cc   
            hh = one+gg*rr+ff*dl+cc*(bet(i)+(half+ac)*gam(i))   
            ee = gg*dl-ff*rr+(bb-half)*bet(i)+(half*bb-two/three+aa)*gam(i)   
            uuu(jj) = odx(i)*hh-ody(i)*ee+gf1*qx1+gf2*qx2   
            vvv(jj) = ody(i)*hh+odx(i)*ee+gf1*qy1+gf2*qy2   
        END DO   
        DO jj = 1,nn   
            j = loc(jj)   
            uplusim(j) = uplusim(j)+uuu(jj)   
            vplusim(j) = vplusim(j)+vvv(jj)   
        END DO   
    END DO u3
     
    do i=1,npt 
    y(i)=yyy(i) 
    dy(i)=ddy(i)  
    x(i)=xxx(i)  
    dx(i)=ddx(i)  
    enddo 




    ddy=dy  
    ddx=dx  
    yyy=y  
    xxx=x
    
    do i=1,npt
    if (y(i) .gt. zero)then     
    y(i)=zero
    endif  
    enddo           
 
    CALL cubic (eu,alp,bet,gam,c,e,g)   
 
    odx = o*dx   
    ody = o*dy   
    edx = eu*dx   
    edy = eu*dy  
 
    
   
    a = x + p1*dx - dy*p1*(alp+p1*(bet+p1*gam))   
    b = y + p1*dy + dx*p1*(alp+p1*(bet+p1*gam))   
    c = odx - ody*(alp+dp1*(bet+tp1*gam))   
    d = ody + odx*(alp+dp1*(bet+tp1*gam))   
    e = x + p2*dx - dy*p2*(alp+p2*(bet+p2*gam))   
    f = y + p2*dy + dx*p2*(alp+p2*(bet+p2*gam))   
    g = odx - ody*(alp+dp2*(bet+tp2*gam))   
    h = ody + odx*(alp+dp2*(bet+tp2*gam))   

	       
   
    
    !Find the velocity at all nodes DUE to i-th node   
u4: DO i = 1,npt   
    
        nn = 0   
   
   
        !Find only the nodes that will contribute to the integral   
        DO j = 1,npt   
            IF ((xxx(j)-x(i))**2+(yyy(j)-y(i))**2 > bdsq) CYCLE	      
            nn = nn+1   
            loc(nn) = j   
        END DO   
   
        DO jj = 1,nn   
            xx(jj) = xxx(loc(jj))   
            yy(jj) = yyy(loc(jj))   
        END DO   
   
        DO jj = 1,nn   
   
            ! Non-singular part (table look-up)   
            qx1 = xx(jj)-a(i)   
            qy1 = yy(jj)-b(i)   
            dd1 = qx1**2+qy1**2   
            rr1 = drdi*SQRT(dd1)   
            k1 = INT(rr1)   
               
            gf1 = (agf(k1)+bgf(k1)*rr1)*(qx1*c(i)+qy1*d(i))/dd1   
               
            qx2 = xx(jj)-e(i)   
            qy2 = yy(jj)-f(i)   
            dd2 = qx2**2+qy2**2   
            rr2 = drdi*SQRT(dd2)   
            k2 = INT(rr2)   
            gf2 = (agf(k2)+bgf(k2)*rr2)*(qx2*g(i)+qy2*h(i))/dd2   
	       
            ! Singular part (explicit)    
            qx = xx(jj)-x(i)   
            qy = yy(jj)-y(i)   
            bb = edx(i)*qx+edy(i)*qy   
            cc = edx(i)*qy-edy(i)*qx   
            bs = bb*bb   
            cs = cc*cc+small   
            ac = two*bb   
            dl = half*LOG(one+(one-ac)/(bs+cs))   
            aa = bs-cs   
            rr = ATAN2(cc,cs-bb*(one-bb))   
            ff = cc*(alp(i)+ac*bet(i)+(three*bs-cs)*gam(i))   
            gg = bb*(alp(i)+(bs-three*cs)*gam(i))+aa*bet(i)-cc   
            hh = one+gg*rr+ff*dl+cc*(bet(i)+(half+ac)*gam(i))   
            ee = gg*dl-ff*rr+(bb-half)*bet(i)+(half*bb-two/three+aa)*gam(i)   
            uu(jj) = odx(i)*hh-ody(i)*ee+gf1*qx1+gf2*qx2   
            vv(jj) = ody(i)*hh+odx(i)*ee+gf1*qy1+gf2*qy2   
        END DO   
        DO jj = 1,nn   
            j = loc(jj)   
            uneg(j) = uneg(j)+uu(jj)   
            vneg(j) = vneg(j)+vv(jj)   
        END DO   
    END DO u4   

  
do j=1,ng
   etag(j)=pi*(float(j)/float(ng+1))
   xg(j)=cos(etag(j))
   yg(j)=zero  
u5: DO i = 1,npt   
    
        nn = 0   
   
   
        !Find only the nodes that will contribute to the integral   
          
            IF ((xg(j)-x(i))**2+(yg(j)-y(i))**2 > bdsq) CYCLE	      
       
 
   
  
   
            ! Non-singular part (table look-up)   
            qx1 = xg(j)-a(i)   
            qy1 = yg(j)-b(i)   
            dd1 = qx1**2+qy1**2   
            rr1 = drdi*SQRT(dd1)   
            k1 = INT(rr1)   
               
            gf1 = (agf(k1)+bgf(k1)*rr1)*(qx1*c(i)+qy1*d(i))/dd1   
               
            qx2 = xg(jj)-e(i)   
            qy2 = yg(jj)-f(i)   
            dd2 = qx2**2+qy2**2   
            rr2 = drdi*SQRT(dd2)   
            k2 = INT(rr2)   
            gf2 = (agf(k2)+bgf(k2)*rr2)*(qx2*g(i)+qy2*h(i))/dd2   
	       
            ! Singular part (explicit)    
            qx = xg(j)-x(i)   
            qy = yg(j)-y(i)   
            bb = edx(i)*qx+edy(i)*qy   
            cc = edx(i)*qy-edy(i)*qx   
            bs = bb*bb   
            cs = cc*cc+small   
            ac = two*bb   
            dl = half*LOG(one+(one-ac)/(bs+cs))   
            aa = bs-cs   
            rr = ATAN2(cc,cs-bb*(one-bb))   
            ff = cc*(alp(i)+ac*bet(i)+(three*bs-cs)*gam(i))   
            gg = bb*(alp(i)+(bs-three*cs)*gam(i))+aa*bet(i)-cc   
            hh = one+gg*rr+ff*dl+cc*(bet(i)+(half+ac)*gam(i))   
            ee = gg*dl-ff*rr+(bb-half)*bet(i)+(half*bb-two/three+aa)*gam(i)   
            ugneg(j) = ugneg(j)+odx(i)*hh-ody(i)*ee+gf1*qx1+gf2*qx2   
            vgneg(j) = vgneg(j)+ody(i)*hh+odx(i)*ee+gf1*qy1+gf2*qy2   

    END DO u5
enddo 

    
    
    do i=1,npt 
    y(i)=yyy(i) 
    dy(i)=ddy(i)  
    x(i)=xxx(i)  
    dx(i)=ddx(i)  
    enddo 
    
    ddy=dy  
    ddx=dx  
    yyy=y  
    xxx=x 
    
    do i=1,npt
    if (y(i) .gt. zero)then     
    y(i)=zero
    endif  
    enddo       
    
    do i=1,npt
      y(i)=-y(i)
    enddo     
  
     CALL cubic (eu,alp,bet,gam,c,e,g)   
 
    odx = o*dx   
    ody = o*dy   
    edx = eu*dx   
    edy = eu*dy  
 

   
    v = zero   
    u = zero   
    a = x + p1*dx - dy*p1*(alp+p1*(bet+p1*gam))   
    b = y + p1*dy + dx*p1*(alp+p1*(bet+p1*gam))   
    c = odx - ody*(alp+dp1*(bet+tp1*gam))   
    d = ody + odx*(alp+dp1*(bet+tp1*gam))   
    e = x + p2*dx - dy*p2*(alp+p2*(bet+p2*gam))   
    f = y + p2*dy + dx*p2*(alp+p2*(bet+p2*gam))   
    g = odx - ody*(alp+dp2*(bet+tp2*gam))   
    h = ody + odx*(alp+dp2*(bet+tp2*gam))   
 
	       
   
    
    !Find the velocity at all nodes DUE to i-th node   
u6: DO i = 1,npt   
    
        nn = 0   
   
   
        !Find only the nodes that will contribute to the integral   
        DO j = 1,npt   
            IF ((xxx(j)-x(i))**2+(yyy(j)-y(i))**2 > bdsq) CYCLE	      
            nn = nn+1   
            loc(nn) = j   
        END DO   
   
        DO jj = 1,nn   
            xx(jj) = xxx(loc(jj))   
            yy(jj) = yyy(loc(jj))   
        END DO   
   
        DO jj = 1,nn   
   
            ! Non-singular part (table look-up)   
            qx1 = xx(jj)-a(i)   
            qy1 = yy(jj)-b(i)   
            dd1 = qx1**2+qy1**2   
            rr1 = drdi*SQRT(dd1)   
            k1 = INT(rr1)   
               
            gf1 = (agf(k1)+bgf(k1)*rr1)*(qx1*c(i)+qy1*d(i))/dd1   
               
            qx2 = xx(jj)-e(i)   
            qy2 = yy(jj)-f(i)   
            dd2 = qx2**2+qy2**2   
            rr2 = drdi*SQRT(dd2)   
            k2 = INT(rr2)   
            gf2 = (agf(k2)+bgf(k2)*rr2)*(qx2*g(i)+qy2*h(i))/dd2   
	       
            ! Singular part (explicit)    
            qx = xx(jj)-x(i)   
            qy = yy(jj)-y(i)   
            bb = edx(i)*qx+edy(i)*qy   
            cc = edx(i)*qy-edy(i)*qx   
            bs = bb*bb   
            cs = cc*cc+small   
            ac = two*bb   
            dl = half*LOG(one+(one-ac)/(bs+cs))   
            aa = bs-cs   
            rr = ATAN2(cc,cs-bb*(one-bb))   
            ff = cc*(alp(i)+ac*bet(i)+(three*bs-cs)*gam(i))   
            gg = bb*(alp(i)+(bs-three*cs)*gam(i))+aa*bet(i)-cc   
            hh = one+gg*rr+ff*dl+cc*(bet(i)+(half+ac)*gam(i))   
            ee = gg*dl-ff*rr+(bb-half)*bet(i)+(half*bb-two/three+aa)*gam(i)   
            uuu(jj) = odx(i)*hh-ody(i)*ee+gf1*qx1+gf2*qx2   
            vvv(jj) = ody(i)*hh+odx(i)*ee+gf1*qy1+gf2*qy2   
        END DO   
        DO jj = 1,nn   
            j = loc(jj)   
            unegim(j) = unegim(j)+uuu(jj)   
            vnegim(j) = vnegim(j)+vvv(jj)   
        END DO   
    END DO u6
   
    do i=1,npt 
    y(i)=yyy(i) 
    dy(i)=ddy(i)  
    x(i)=xxx(i)  
    dx(i)=ddx(i)  
    enddo 





   

    
    
    do i=1,ng
      fg(i)=two*(ug(i)-ugneg(i))*sin(etag(i))       
    enddo
    
    do i=1,ng
      bg(i)=zero
       do j=1,ng
         bg(i)=bg(i)+fg(j)*sin(float(i)*etag(j))
       enddo
     bg(i)=two*bg(i)/float(ng+1)
      
     
    enddo
     
    
    do i=1,npt
     z(i)=x(i)+ci*y(i)
     zeta(i)=z(i)+ci*sqrt(-z(i)**2+one)
     sumcomp=zero+ci*zero
!     do j=1,ng
      do j=ng,1,-1
       if  (y(i) .ge. zero)then
       sumcomp=sumcomp+bg(j)*zeta(i)**(-j-1)
       else
       sumcomp=sumcomp+bg(j)*zeta(i)**(j-1)
       endif
     enddo
       
      ucomp(i)=(one/(one-one/zeta(i)**2))*sumcomp
      ugap(i)=real(ucomp(i))
      vgap(i)=-aimag(ucomp(i))
     enddo 
    
       	   


   
    !calculate total velocity at each node   
    DO i=1,npt
      
     if (y(i) .gt. zero) then      
      u(i)=uplus(i)+uplusim(i)+ugap(i)+ &
          (flux/pi)*real(one/(ci*sqrt(one-z(i)**2)))+stream+&
	   (beta/two)*real(one+z(i)/(ci*sqrt(one-z(i)**2))) 	     	        
      v(i)=vplus(i)+vplusim(i)+vgap(i)- &
           (flux/pi)*aimag(one/(ci*sqrt(one-z(i)**2)))-&
           (beta/two)*aimag(one+z(i)/(ci*sqrt(one-z(i)**2)))
     else
      u(i)=uneg(i)+unegim(i)+ugap(i)+ &
         (flux/pi)*real(one/(ci*sqrt(one-z(i)**2)))+stream+&
	  (beta/two)*real(one+z(i)/(ci*sqrt(one-z(i)**2)))	        
      v(i)=vneg(i)+vnegim(i)+vgap(i)-&
         (flux/pi)*aimag(one/(ci*sqrt(one-z(i)**2)))-&
	 (beta/two)*aimag(one+z(i)/(ci*sqrt(one-z(i)**2)))      
     endif 
         
    END DO  
  
  
     


    DEALLOCATE(xx,yy,uu,uuu,vv,vvv,loc,odx)
 
    DEALLOCATE (ody,edx,edy,a,b,c,d)

      DEALLOCATE(e,f,g,h)

     DEALLOCATE( &
     eu,alp,bet,gam,uplus,vplus,uplusim,vplusim,uneg,vneg,unegim,vnegim) 

      DEALLOCATE(xxx,yyy,ddy,ddx) 
		  
 
      
          
END SUBROUTINE uv   

   
END MODULE uv_   
   
   
   
   
