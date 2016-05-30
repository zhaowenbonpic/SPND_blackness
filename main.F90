program main
common /CPR000/ a(1000)
common /CPR001/ nin,NOT,nx1(2)
common /CPR002/ mmax,ig,nx2(4)
1010 FORMAT(2i6)
nin=15
not=16
open(nin,file='inp')
read(nin,*) mmax,ig
!
! ADdressing
!
lrad = 1
lvol = lrad+mmax
lsig = lvol+mmax
lsigs = lsig +mmax
lsigr = lsigs +mmax
lsig1 = lsigr +mmax
lsec = lsig1 +mmax
lsigv = lsec +mmax
!
lp = lsigv +mmax
lgam = lp + mmax*mmax
lyi = lgam  !?
lxik = lyi +mmax
lqi = lxik +mmax*mmax
lvoly = lqi+mmax
lfluq = lvoly
lfluj = lfluq +mmax
lflux = lfluj +mmax
lremo = lflux +mmax
last = lremo +mmax
!
call prepar(A(lrad),A(LVOL),A(LSIG),A(LSIGS),A(LSIGR),A(LSIG1),A(LSEC),A(LSIGV))
print*, "prepar end"
CALL COPRAN(MMAX,IG,A(LRAD),A(LSIG),A(LP),A(LGAM))
print*, "COPRAN end"
CALL GAUSCP(MMAX,A(LSEC),A(LSIG1),A(LSIGV),A(LP),A(LGAM),A(LXIK))
print*, "GAUSCP end"
CALL CALFLU(A(LVOL),A(LSIGR),A(LYI),A(LXIK),A(LQI),A(LVOLY), &
&           A(LFLUQ),A(LFLUJ),A(LFLUX),A(LREMO))
print*, "CALFLU end"
STOP
END



 




 





