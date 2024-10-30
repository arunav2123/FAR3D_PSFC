	subroutine bigf(f,g,itypeg,h,itypeh,sx1,sx2,c1,c2)

		use mpi
		use param
		use processor
		use var_para
		use domain
		implicit none

		integer :: itypeg,itypeh,l
		real(IDP) :: c1,c2
		real(IDP), dimension(mj_start:,0:) :: f,g,h,sx1,sx2
		
		interface
			subroutine dbydth_par(d,a,ltype,c1,c2,k)
				use param
				use var_para
				implicit none
				integer :: k,ltype
				real(IDP) :: c1,c2
				real(IDP), dimension(mj_start:,0:) :: a,d
			end subroutine dbydth_par
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
		end interface

		call dbydr_par(sx2,g,0.0_IDP,-1.0_IDP,0)
		call mult(sx1,h,itypeh,sx2,itypeg,0.0_IDP,1.0_IDP)
		call dbydth_par(sx2,sx1,itypeg*itypeh,0.0_IDP,1.0_IDP,1)
		if (myPE == 0) then
			do l=1,lmax
				sx2(0,l)=0.
			end do
		end if

!$OMP PARALLEL DO
		do l=1,lmax
			f(mj_start:mj_end,l)=c1*f(mj_start:mj_end,l)+c2*sx2(mj_start:mj_end,l)
		end do
!$OMP END PARALLEL DO

		call dbydth_par(sx1,g,itypeg,0.0_IDP,1.0_IDP,0)
		call mult(sx2,h,itypeh,sx1,-itypeg,0.0_IDP,1.0_IDP)
		if (myPE == 0) then
			do l=1,lmax
				if (mm(l) == 0) f(0,l)=f(0,l)+c2*2.*rinv(1)*sx2(1,l)
			end do
		end if
!$OMP PARALLEL DO
		do l=1,lmax
			sx2(mj_start:mj_end,l)=r(mj_start:mj_end)*sx2(mj_start:mj_end,l)
		end do
!$OMP END PARALLEL DO
		call dbydr_par(sx1,sx2,0.0_IDP,1.0_IDP,0)
!$OMP PARALLEL DO
		do l=1,lmax
			f(mj_start:mj_end,l)=f(mj_start:mj_end,l)+c2*rinv(mj_start:mj_end)*sx1(mj_start:mj_end,l)
		end do
!$OMP END PARALLEL DO

	end subroutine bigf

	subroutine divv(f,c1,c2)

		use param
		use var_para
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		integer :: l
		real(IDP) :: c1,c2
		real(IDP), dimension(mj_start:,0:) :: f
		
		interface
			subroutine dbydth_par(d,a,ltype,c1,c2,k)
				use param
				use var_para
				implicit none
				integer :: k,ltype
				real(IDP) :: c1,c2
				real(IDP), dimension(mj_start:,0:) :: a,d
			end subroutine dbydth_par
			subroutine dbydr_par(d,a,c1,c2,k)
				use param
				use var_para
				implicit none
				integer :: k
				real(IDP) :: c1,c2
				real(IDP), dimension(mj_start:,0:) :: a,d
			end subroutine dbydr_par
			subroutine dbydzt_par(d,a,ltype,c1,c2)
				use param
				use var_para
				implicit none
				integer :: ltype
				real(IDP) :: c1,c2
				real(IDP), dimension(mj_start:,0:) :: a,d
			end subroutine dbydzt_par
			subroutine grdpar(d,a,ltype,c1,c2)
				use param
				use var_para
				implicit none
				integer :: ltype
				real(IDP) :: c1,c2
				real(IDP), dimension(mj_start:,0:) :: a,d
			end subroutine grdpar
			subroutine grpareq(d,a,ltype,c1,c2)
				use param
				implicit none
				integer :: ltype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine grpareq
			subroutine dbydr0(d,a,c1,c2,k)
				use param
				implicit none
				integer :: k
				real(IDP) :: c1,c2
				real(IDP), dimension(0:) :: a,d
			end subroutine dbydr0
			subroutine eqtodyn_par(adyn,aeq,c1,c2)
				use param
				use var_para
				implicit none
				real(IDP) :: c1,c2
				real(IDP), dimension(mj_start:,0:) :: adyn
				real(IDP), dimension(0:,0:) :: aeq
			end subroutine eqtodyn_par
			subroutine mult(f,g,itypeg,h,itypeh,c1,c2)
				use param
				use var_para
				implicit none
				integer :: itypeg,itypeh
				real(IDP) :: c1,c2
				real(IDP), dimension(mj_start:,0:) :: f,g,h
			end subroutine mult
		end interface

		sd2=feq/(feq-qqinv*cureq)
		sd3=cureq/(feq-qqinv*cureq)
		do l=1,leqmax
			sceq1(:,l)=-(sd2*djtoj(:,l)-rinv*sd3*djzoj(:,l))
		end do
		call eqtodyn_par(sc2,sceq1,0.0_IDP,1.0_IDP)
		call dbydr_par(sc1,phi,0.0_IDP,1.0_IDP,0)
		call mult(sc3,sc1,-1,sc2,-2,0.0_IDP,1.0_IDP)
		do l=1,leqmax
			sceq1(:,l)=-(r*dbsjzoj(:,l)-sd2*djroj(:,l))
		end do
		call eqtodyn_par(sc2,sceq1,0.0_IDP,1.0_IDP)
		call dbydth_par(sc1,phi,-1,0.0_IDP,1.0_IDP,0)
		call mult(sc3,sc1,1,sc2,2,1.0_IDP,1.0_IDP)
		call dbydr0(sd4,sd2,0.0_IDP,1.0_IDP,0)
!$OMP PARALLEL DO
		do l=1,lmax
			sc3(:,l)=sc3(:,l)+sd4(mj_start:mj_end)*sc1(:,l)
		end do
!$OMP END PARALLEL DO
		do l=1,leqmax
			sceq1(:,l)=-(rinv*sd3*djroj(:,l)-dbsjtoj(:,l))
		end do
		call eqtodyn_par(sc2,sceq1,0.0_IDP,1.0_IDP)
		call dbydzt_par(sc1,phi,-1,0.0_IDP,1.0_IDP)
		call mult(sc3,sc1,1,sc2,2,1.0_IDP,1.0_IDP)
		call dbydr0(sd4,sd3,0.0_IDP,1.0_IDP,0)
!$OMP PARALLEL DO
		do l=1,lmax
			sc3(:,l)=sc3(:,l)-rinv(mj_start:mj_end)*sd4(mj_start:mj_end)*sc1(:,l)
		end do
!$OMP END PARALLEL DO

!  parallel thermal velocity term		

		do l=1,leqmax
			sceq1(:,l)=bmod(:,l)/(feq-qqinv*cureq)
		end do
		call eqtodyn_par(sc2,sceq1,0.0_IDP,1.0_IDP)
		call grdpar(sc1,vthprlf,-1,0.0_IDP,-1.0_IDP)
		call mult(sc3,sc1,1,sc2,2,1.0_IDP,1.0_IDP)
		call grpareq(sceq2,bmod,1,0.0_IDP,1.0_IDP)
		do l=1,leqmax
			sceq1(:,l)=sceq2(:,l)/(feq-qqinv*cureq)
		end do
		call eqtodyn_par(sc2,sceq1,0.0_IDP,1.0_IDP)
		call mult(sc3,vthprlf,-1,sc2,-2,1.0_IDP,1.0_IDP)

		call mult(f,pr,1,sc3,1,c1,c2)

	end subroutine divv

	subroutine chtype(a,b)
	
		use param
		use var_para
		use domain
		implicit none

		integer :: l
		real(IDP), dimension(mj_start:,0:) :: a,b
		
		a=0.0_IDP
!$OMP PARALLEL DO
		do l=1,lmaxn
			a(:,lln(lo(l)))=b(:,lln(l))
		end do
!$OMP END PARALLEL DO
		a(:,0)=0.0_IDP
		
	end subroutine chtype