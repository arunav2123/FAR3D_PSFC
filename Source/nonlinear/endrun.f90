	subroutine endrun

		use param
		use processor
		use var_para
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		integer :: i,j,l,lwrt,l2t,mnum
		real(IDP) :: scmn,scmx,sclmn,sclmx
		real(IDP), dimension(:), allocatable :: scnorm
		character(len=1) :: t
		character(len=12) :: confil
		character(len=32) :: formatt='("r",1200(a1,i4,"/",i4))'
		character(len=32) :: formatv='(1pe13.6,1200(a1,1pe15.8))'

		interface
			subroutine dbydth(d,a,ltype,c1,c2,k)
				use param
				implicit none
				integer :: k,ltype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine dbydth
			subroutine dbydr(d,a,c1,c2,k)
				use param
				implicit none
				integer :: k
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine dbydr
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
			subroutine trnsfr0e(f,ftype)
				use param
				use var_para
				implicit none
				integer :: ftype
				real(IDP), dimension(mj_start:,0:) :: f
			end subroutine trnsfr0e
			subroutine delstar(ss,ff,itypf,wk1,wk2,wk3,c1,c2)
				use param
				use var_para
				implicit none
				integer :: itypf
				real(IDP) :: c1,c2
				real(IDP), dimension(mj_start:,0:) :: ss,ff,wk1,wk2,wk3
			end subroutine delstar			
		end interface

		t=char(9)

!  store uzt-values into sc4

		do l=1,lmax
			sc4(:,l)=uzt(mj_start:mj_end,l)
		end do

		allocate (scnorm(nnum))
		scnorm=1.0_IDP

		if (myPE == 0) then

!  vr up
			call dbydth(uzt,phi,-1,0.0_IDP,-eps,0)

			do i=1,nnum
				if (i < nst) then
					mnum=mnumn(i)
				else
					mnum=2*mnumn(i)
				end if
				if (nonlin == 0) then
					scmn=0.
					scmx=0.
					do l2t=1,mnum
						l=l2t+lnumn(i-1)
						sclmn=minval(uzt(:,lln(l)))
						sclmx=maxval(uzt(:,lln(l)))
						scmn=min(scmn,sclmn)
						scmx=max(scmx,sclmx)
					end do
					if (scmx > abs(scmn)) then
						scnorm(i)=scmx
					else
						scnorm(i)=scmn
					endif
				end if
				if(scnorm(i) /= 0.0_IDP) then
					do l2t=1,mnum
						l=l2t+lnumn(i-1)
						psi(:,lln(l))=psi(:,lln(l))/scnorm(i)
						where (abs(psi(:,lln(l))) < 1.e-50_IDP) psi(:,lln(l))=0
						phi(:,lln(l))=phi(:,lln(l))/scnorm(i)
						where (abs(phi(:,lln(l))) < 1.e-50_IDP) phi(:,lln(l))=0
						pr(:,lln(l))=pr(:,lln(l))/scnorm(i)
						where (abs(pr(:,lln(l))) < 1.e-50_IDP) pr(:,lln(l))=0
						nf(:,lln(l))=nf(:,lln(l))/scnorm(i)
						where (abs(nf(:,lln(l))) < 1.e-50_IDP) nf(:,lln(l))=0
						vprlf(:,lln(l))=vprlf(:,lln(l))/scnorm(i)
						where (abs(vprlf(:,lln(l))) < 1.e-50_IDP) vprlf(:,lln(l))=0
						vthprlf(:,lln(l))=vthprlf(:,lln(l))/scnorm(i)
						where (abs(vthprlf(:,lln(l))) < 1.e-50_IDP) vthprlf(:,lln(l))=0
						uzt(:,lln(l))=uzt(:,lln(l))/scnorm(i)
						where (abs(uzt(:,lln(l))) < 1.e-50_IDP) uzt(:,lln(l))=0
					end do
				end if
			end do
			lwrt=abs(lplots)
			if (lwrt > lmax) lwrt=lmax
			write(confil,'("phi_",2a2)') numrun(1),numrun(2)
			open(unit=92,file=confil,recl=19384)
			write(92,formatt) (t,mm(l),nn(l),l=1,lwrt)
			do j=0,mj
				write(92,formatv) r(j),(t,phi(j,l),l=1,lwrt)
			end do
			close(92)
			write(confil,'("psi_",2a2)') numrun(1),numrun(2)
			open(unit=92,file=confil,recl=19384)
			write(92,formatt) (t,mm(l),nn(l),l=1,lwrt)
			do j=0,mj
				write(92,formatv) r(j),(t,psi(j,l),l=1,lwrt)
			end do
			close(92)
			write(confil,'("pr_",2a2)') numrun(1),numrun(2)
			open(unit=92,file=confil,recl=19384)
			write(92,formatt) (t,mm(l),nn(l),l=1,lwrt)
			do j=0,mj
				write(92,formatv) r(j),(t,pr(j,l),l=1,lwrt)
			end do
			close(92)
			write(confil,'("vr_",2a2)') numrun(1),numrun(2)
			open(unit=92,file=confil,recl=19384)
			write(92,formatt) (t,mm(l),nn(l),l=1,lwrt)
			do j=0,mj
				write(92,formatv) r(j),(t,uzt(j,l),l=1,lwrt)
			end do
			close(92)
			write(confil,'("nf_",2a2)') numrun(1),numrun(2)
			open(unit=92,file=confil,recl=19384)
			write(92,formatt) (t,mm(l),nn(l),l=1,lwrt)
			do j=0,mj
				write(92,formatv) r(j),(t,nf(j,l),l=1,lwrt)
			end do
			close(92)
			write(confil,'("vprlf_",2a2)') numrun(1),numrun(2)
			open(unit=92,file=confil,recl=19384)
			write(92,formatt) (t,mm(l),nn(l),l=1,lwrt)
			do j=0,mj
				write(92,formatv) r(j),(t,vprlf(j,l),l=1,lwrt)
			end do
			close(92)
			write(confil,'("vthprlf_",2a2)') numrun(1),numrun(2)
			open(unit=92,file=confil,recl=19384)
			write(92,formatt) (t,mm(l),nn(l),l=1,lwrt)
			do j=0,mj
				write(92,formatv) r(j),(t,vthprlf(j,l),l=1,lwrt)
			end do
			close(92)

!  vth up
			call dbydr(uzt,phi,0.0_IDP,eps,0)
			write(confil,'("vth_",2a2)') numrun(1),numrun(2)
			open(unit=92,file=confil,recl=19384)
			write(92,formatt) (t,mm(l),nn(l),l=1,lwrt)
			do j=0,mj
				write(92,formatv) r(j),(t,uzt(j,l),l=1,lwrt)
			end do
			close(92)
		end if

!  br up
		call dbydth_par(sc3,psi,1,0.0_IDP,-eps,0)
!  bth up
		call dbydr_par(sc5,psi,0.0_IDP,eps,0)

		call eqtodyn_par(sc8,sqgi,0.0_IDP,1.0_IDP)

		call mult(uzt,sc8,1,sc3,-1,0.0_IDP,1.0_IDP)
		call trnsfr0e(uzt,-1)
		if (myPE == 0) then
			write(confil,'("br_",2a2)') numrun(1),numrun(2)
			open(unit=92,file=confil,recl=19384)
			write(92,formatt) (t,mm(l),nn(l),l=1,lwrt)
			do j=0,mj
				write(92,formatv) r(j),(t,uzt(j,l),l=1,lwrt)
			end do
			close(92)
		end if

		call mult(uzt,sc8,1,sc5,1,0.0_IDP,1.0_IDP)
		call trnsfr0e(uzt,1)
		if (myPE == 0) then
			write(confil,'("bth_",2a2)') numrun(1),numrun(2)
			open(unit=92,file=confil,recl=19384)
			write(92,formatt) (t,mm(l),nn(l),l=1,lwrt)
			do j=0,mj
				write(92,formatv) r(j),(t,uzt(j,l),l=1,lwrt)
			end do
			close(92)
		end if

		call delstar(uzt,psi,1,sc2,sc3,sc5,0.0_IDP,1.0_IDP)
		uzt=eps*uzt
		call trnsfr0e(uzt,1)
		if (myPE == 0) then
			write(confil,'("curzt_",2a2)') numrun(1),numrun(2)
			open(unit=92,file=confil,recl=19384)
			write(92,formatt) (t,mm(l),nn(l),l=1,lwrt)
			do j=0,mj
				write(92,formatv) r(j),(t,uzt(j,l),l=1,lwrt)
			end do
			close(92)
		end if

!  restore uzt-values

		do l=1,lmax
			uzt(mj_start:mj_end,l)=sc4(:,l)
		end do
		call trnsfr0e(uzt,-1)
		if (myPE == 0) then
			do i=1,nnum
				if (i < nst) then
					mnum=mnumn(i)
				else
					mnum=2*mnumn(i)
				end if
!				if (nonlin == 0) then
!					scmn=0.
!					scmx=0.
!					do l2t=1,mnum
!						l=l2t+lnumn(i-1)
!						sclmn=minval(uzt(:,lln(l)))
!						sclmx=maxval(uzt(:,lln(l)))
!						scmn=min(scmn,sclmn)
!						scmx=max(scmx,sclmx)
!					end do
!					if (scmx > abs(scmn)) then
!						scnorm(i)=scmx
!					else
!						scnorm(i)=scmn
!					endif
!				end if
				if(scnorm(i) /= 0.0_IDP) then
					do l2t=1,mnum
						l=l2t+lnumn(i-1)
						uzt(:,lln(l))=uzt(:,lln(l))/scnorm(i)
						where (abs(uzt(:,lln(l))) < 1.e-50_IDP) uzt(:,lln(l))=0
					end do
				end if
			end do
			write(confil,'("uzt_",2a2)') numrun(1),numrun(2)
			open(unit=92,file=confil,recl=19384)
			write(92,formatt) (t,mm(l),nn(l),l=1,lwrt)
			do j=0,mj
				write(92,formatv) r(j),(t,uzt(j,l),l=1,lwrt)
			end do
			close(92)
		end if

	end subroutine endrun