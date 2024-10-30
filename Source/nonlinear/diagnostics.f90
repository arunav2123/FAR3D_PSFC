!
!     Routine for high frequency (every time step) output of a few quantities
!
      subroutine hifreq
      use var_para
      use cotrol
      use domain
      use equil
      use dynamo
      use scratch
      use param
      real(IDP) :: bthta, btht1, btht2, btht3, btht4, phi00,psi00,pr00,vprl00
      integer :: rd2, rd4, rd6, rd8, rda
      character(len=23) :: formath='(1pe13.6,9(a1,1pe15.8))'
      character(len=1) :: t
		interface
			subroutine dbydr_par(d,a,c1,c2,k)
				use param
				use var_para
				implicit none
				integer :: k
				real(IDP) :: c1,c2
				real(IDP), dimension(mj_start:,0:) :: a,d
			end subroutine dbydr_par
			subroutine mult(f,g,itypeg,h,itypeh,c1,c2)
				use param
				use var_para
				implicit none
				integer :: itypeg,itypeh
				real(IDP) :: c1,c2
				real(IDP), dimension(mj_start:,0:) :: f,g,h
			end subroutine mult
			subroutine eqtodyn_par(adyn,aeq,c1,c2)
				use param
				use var_para
				implicit none
				real(IDP) :: c1,c2
				real(IDP), dimension(mj_start:,0:) :: adyn
				real(IDP), dimension(0:,0:) :: aeq
			end subroutine eqtodyn_par
			subroutine trnsfr0e(f,ftype)
				use param
				use var_para
				implicit none
				integer :: ftype
				real(IDP), dimension(mj_start:,0:) :: f
			end subroutine trnsfr0e
			subroutine quadq(n,x,y,x0,b,c,d,result)
				use param
				implicit none
				integer :: n
				real(IDP) :: x0,result
				real(IDP), dimension(:) :: x,y,b,c,d
			end subroutine quadq
	       end interface
      call dbydr_par(sc1,psi,0.0_IDP,1.0_IDP,0)
      call eqtodyn_par(sc2,bmod,0.0_IDP,1.0_IDP)
      call mult(sc3,sc2,2,sc1,1,0.0_IDP,1.0_IDP)
      do 13 l=1,lmax
      sc2(:,l)=uzt(mj_start:mj_end,l)
      uzt(mj_start:mj_end,l)=sc3(:,l)
  13  continue
      call trnsfr0e(uzt,1)
      if (myPE == 0) then
      t=char(9)
      bthta = 0.0_IDP; btht1 = 0.0_IDP; btht2 = 0.0_IDP
      btht3 = 0.0_IDP; btht4 = 0.0_IDP
      rda=mj;rd2=mj/5;rd4=2*mj/5;rd6=3*mj/5;rd8=4*mj/5
      ltype = +1
      do 10 l=1,lmax
      lsign=signl(l)
      if(ltype*lsign .gt. 0) bthta=bthta+uzt(rda,l)*eps/(feq(rda)-qqinv(rda)*cureq(rda))
      if(ltype*lsign .gt. 0) btht1=btht1+uzt(rd2,l)*eps/(feq(rd2)-qqinv(rd2)*cureq(rd2))
      if(ltype*lsign .gt. 0) btht2=btht2+uzt(rd4,l)*eps/(feq(rd4)-qqinv(rd4)*cureq(rd4))
      if(ltype*lsign .gt. 0) btht3=btht3+uzt(rd6,l)*eps/(feq(rd6)-qqinv(rd6)*cureq(rd6))
      if(ltype*lsign .gt. 0) btht4=btht4+uzt(rd8,l)*eps/(feq(rd8)-qqinv(rd8)*cureq(rd8))
  10  continue
      do 11 j=0,mj
      sd1(j) = phi(j,l0)*r(j)
      sd2(j) = psi(j,l0)*r(j)
      sd3(j) = nf(j,l0)*r(j)
      sd4(j) = vprlf(j,l0)*r(j)
  11  continue
      mjp = mj + 1
      call quadq(mjp,r,sd1,1._IDP,sd5,sd6,sd7,phi00)
      call quadq(mjp,r,sd2,1._IDP,sd5,sd6,sd7,psi00)
      call quadq(mjp,r,sd3,1._IDP,sd5,sd6,sd7,pr00)
      call quadq(mjp,r,sd4,1._IDP,sd5,sd6,sd7,vprl00)
      write(82,formath) time,t,btht1,t,btht2,t,btht3,t,btht4,  &
		  t,bthta,t,phi00,t,psi00,t,pr00,t,vprl00
      end if
      do 15 l=1,lmax
      uzt(mj_start:mj_end,l)=sc2(:,l)
  15  continue
      return
      end subroutine hifreq

