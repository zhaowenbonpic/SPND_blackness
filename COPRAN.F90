SUBROUTINE COPRAN(MMAX,IG,R,SIG,P,GAM)
! collision probabilities in annular geometry
! by Ingvar Carvik
! MMAX = number of annular region
! IG   = number of points in integration
! R()  = radii in increasing order
! p(i,j) = vol(i)*sig(i)*(colprob from region I to J)
! GAM()  = partial blacknesses
! TAU()  = optcal distance
! SIGV() = sigma*volume
REAL KI3 !function need write
dimension R(MMAX),SIG(MMAX),P(MMAX,MMAX),GAM(mmax)
dimension GJC(36),GJP(6,6),GJW(6,6)
dimension R2(mmax),SIGV(mmax),TAU(20),TX(21)
EQUIVALENCE (GJC(1),GJP(1,1),GJW(1,1)),(TX(2),TAU(1))
!  Gauss-Jacobi integration
!  x(i) = 1 - x(i)**2
!  w(i) = 4*w(i)
DATA GJC &
	/.00000000,.55555556,.87393877,.95491150,.98046718,.99029084, &
	 2.0000000,.00000000,.28606124,.65127016,.82660307,.90725799, &
	 .72783448,1.27216522,.0000000,.16932809,.47704397,.68412769, &
	 .27930792,.91696442,.80372766,.00000000,.11094751,.35681753, &
	 .12472388,.51939018,.81385828,.54202764,.00000000,.07803490, &
	 .06299166,.29563548,.58554794,.66869856,.38712636,.00000000/

DATA PI / 3.1415927/

N = MMAX
y = 0.0d0
do i=1,n
	R2(i) = R(i)*R(i)
	SIGV(i) = sig(i)*pi*(r2(i)-y)
	y = r2(i)
	do j = 1,n
          p(i,j) = 0.0d0
	enddo
enddo

do i = 1,n
	rstart = 0.0
	if(i.gt.1) rstart = r(i-1)
	dr = r(i)-rstart
	do k = 1,ig
	  y = rstart+dr*GJP(IG+1,K)
	  fac = dr*GJW(k,IG+1)
	  TAU(i-1) = 0.0
	  tplus = 0.0
	  Y2 = Y*Y
	  do j = i,n
	    tminus = tplus
	    tplus = sqrt(r2(j)-y2)
	    tau(j) = tau(j-1)+sig(j)*(tplus-tminus)
	    do l=i,j
	      p(l,j) = p(l,j) + FAC*(KI3(TAU(J)+TAU(L))-KI3(TAU(J)-TAU(L)))
	    enddo
	  enddo
	enddo
enddo
!print*, P

do i = 1,n
	j = n-i+1
	do k=j,n
	  l = j+n-k
	  if(l.eq.1) goto 30
	  if(l.eq.j) p(j,l-1)=p(l-1,j)
	  p(j,l)=p(j,l)-p(j,l-1)
	  if(j.eq.1) goto 30
	  p(j,l) = p(j,l) -p(j-1,l)+p(j-1,l-1)
 30	    p(l,j) = p(j,l)
	enddo
	p(j,j) = p(j,j)+sigv(j)
enddo

s = 2.0/(PI*r(N))
DO I = 1,N
	SUM = SIGV(I)
	DO J = 1,N
	  SUM = SUM - p(I,J)
	ENDDO
	gam(i) = s*sum
ENDDO
write(*,*) gam
RETURN
END

real FUNCTION KI3(XX)
! 
! 3 and 4 term power series corresponding to chebyshev expansion
! intevals of 0.05 in range from 0.0 to 0.5, coeffcients A
!             0.10 in range from 0.5 to 1.0, A()
!             0.40 in range from 1.0 to 5.0,B()
!             0.80 in range from 5.0 to 9.0, B()
dimension a(45),b(51),indexa(20),indexb(20)
data a / 0.7853961,-0.9990226,0.7266088, &
        0.7852024,-0.9912340,0.6466375, &
        0.7845986,-0.9791293,0.5856605, &
        0.7834577,-0.9638914,0.5346648, &
        0.7817094,-0.9463843,0.4907827, &
        0.7793031,-0.9271152,0.4521752, &
        0.7762107,-0.9064822,0.4177388, &
        0.7724519,-0.8849865,0.3869945, &
        0.7679903,-0.8626685,0.3590753, &
        0.7628988,-0.8400133,0.3338676, &
        0.7540982,-0.8054172,0.2998569, &
        0.7401279,-0.7587821,0.2609154, &
        0.7239594,-0.7125290,0.2278226, &
        0.7058777,-0.6672761,0.1994999, &
        0.6861762,-0.6234536,0.1751248/
data b / 0.7247294,-0.7538355,0.3203223,-5.337485E-2, &
        0.6663720,-0.6279752,0.2295280,-3.146833E-2, &
        0.5956163,-0.5094124,0.1631667,-1.906198E-2, &
        0.5191031,-0.4046007,0.1152418,-1.174752E-2, &
        0.4425954,-0.3159648,8.097913E-2,-7.328415E-3, &
        0.3703178,-0.2434341,5.669960E-2,-4.617254E-3, &
        0.1684022,-7.158569E-2,7.923547E-3, &
        0.1278307,-5.016344E-2,5.095111E-3, &
        9.611422E-2,-3.501524E-2,3.286040E-3, &
        7.170491E-2,-2.437465E-2,2.126242E-3, &
        4.616317E-2,-1.425519E-2,1.123687E-3, &
        2.475115E-2,-6.810124E-3,4.762937E-4, &
        1.302864E-2,-3.232035E-3,2.031843E-4, &
        6.749972E-3,-1.524126E-3,8.701440E-5, &
        3.454768E-3,-7.157367E-4,3.742673E-5 /

data  indexa / 3,6,9,12,15,18,21,24,27,30, &
              33,33,36,36,39,39,42,42,45,45/
data  indexb / 4,8,12,16,20,24,27,30,33,36, &
             39,39,42,42,45,45,48,48,51,51/

!
!  inteval 0.0 - 1.0
x = ABS(xx)
if (x.ge.1.0d0) goto 10
I = ifix(20.0*x) + 1
i = indexa(i)
KI3 = x*(x*A(i)+a(i-1))+a(i-2)
return
!  inteval 1.0 - 3.4
10    i = ifix(2.5*(x-1.0)) + 1
if (x.ge.3.4) goto 20
i = indexb(i)
ki3 = x*(x*(x*b(i)+b(i-1))+b(i-2))+b(i-3)
return
20    if(x.ge.9.0) goto 30
i = indexb(i)
ki3 = x*(x*b(i)+b(i-1))+b(i-2)
return
30    ki3 = 0.0
return
end

   	   




