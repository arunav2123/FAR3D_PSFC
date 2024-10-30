	subroutine delstar(ss,ff,itypf,wk1,wk2,wk3,c1,c2)

		use mpi
		use param
		use processor
		use var_para
		use cotrol
		use domain
		use equil
		implicit none

		integer :: itypf,l,j,tag,ierr
		real(IDP) :: c1,c2,xm
		real(IDP), dimension(mj_start:,0:) :: ss,ff,wk1,wk2,wk3
		real(IDP), dimension(lmax) :: ssend,srecv
		integer, dimension(MPI_STATUS_SIZE) :: status

		interface
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
		end interface

		ff(mj_start:mj_end,0)=0. 
		ss(mj_start:mj_end,0)=0. 
		wk1=0. 
		wk2=0. 
		wk3=0. 

		if (myPE == 0) then
!$OMP PARALLEL DO PRIVATE(j)
			do l=1,lmax
				do j=1,mj_end-1
					wk2(j,l)=del2cp(j)*(ff(j+1,l)-ff(j,l))+del2cm(j)*(ff(j-1,l)-ff(j,l))
				end do
			end do
!$OMP END PARALLEL DO
		else if (myPE == numPEsm1) then
!$OMP PARALLEL DO PRIVATE(j)
			do l=1,lmax
				do j=mj_start+1,mjm1
					wk2(j,l)=del2cp(j)*(ff(j+1,l)-ff(j,l))+del2cm(j)*(ff(j-1,l)-ff(j,l))
				end do
			end do
!$OMP END PARALLEL DO
		else
!$OMP PARALLEL DO PRIVATE(j)
			do l=1,lmax
				do j=mj_start+1,mj_end-1
					wk2(j,l)=del2cp(j)*(ff(j+1,l)-ff(j,l))+del2cm(j)*(ff(j-1,l)-ff(j,l))
				end do
			end do
!$OMP END PARALLEL DO
		end if
		call eqtodyn_par(wk3,gttoj,0.0_IDP,1.0_IDP)
		call mult(ss,wk2,itypf,wk3,2,c1,c2)
		call dbydr_par(wk2,wk3,0.0_IDP,1.0_IDP,2)
		call eqtodyn_par(wk3,grtoj,0.0_IDP,-1.0_IDP)
		call dbydth_par(wk1,wk3,-1,0.0_IDP,1.0_IDP,2)
!$OMP PARALLEL DO
		do l=1,lmax
			wk1(:,l)=wk1(:,l)+wk2(:,l)
		end do
!$OMP END PARALLEL DO
		call dbydr_par(wk2,ff,0.0_IDP,1.0_IDP,0)
		if (myPE == 0) then
			do l=1,lmax
				wk2(0,l)=0.
			end do
		end if
		call mult(ss,wk2,itypf,wk1,2,1.0_IDP,c2)
		call dbydth_par(wk1,wk2,itypf,0.0_IDP,1.0_IDP,0)
		if (myPE == 0) then
			do l=1,lmax
				wk1(0,l)=0.
			end do
		end if
		call mult(ss,wk1,-itypf,wk3,-2,1.0_IDP,2.*c2)
		call dbydr_par(wk1,wk3,0.0_IDP,1.0_IDP,0)
		call eqtodyn_par(wk3,grroj,0.0_IDP,1.0_IDP)
		call dbydth_par(wk2,wk3,1,0.0_IDP,1.0_IDP,2)
!$OMP PARALLEL DO
		do l=1,lmax
			wk1(:,l)=wk1(:,l)+wk2(:,l)
		end do
!$OMP END PARALLEL DO
		call dbydth_par(wk2,ff,itypf,0.0_IDP,1.0_IDP,0)
		if (myPE == 0) then
			do l=1,lmax
				wk2(0,l)=0.
			end do
		end if
		call mult(ss,wk2,-itypf,wk1,-2,1.0_IDP,c2)

		if (myPE == 0) then
!$OMP PARALLEL DO
			do l=1,lmax
				xm=mm(l)
				wk2(1:mj_end,l)=-(rinv(1:mj_end)*xm)**2*ff(1:mj_end,l)
			end do
!$OMP END PARALLEL DO
		else
!$OMP PARALLEL DO
			do l=1,lmax
				xm=mm(l)
				wk2(mj_start:mj_end,l)=-(rinv(mj_start:mj_end)*xm)**2*ff(mj_start:mj_end,l)
			end do
!$OMP END PARALLEL DO
		end if
		call mult(ss,wk3,2,wk2,itypf,1.0_IDP,c2)

		if (myPE == numPEsm1) then
!$OMP PARALLEL DO
			do l=1,lmax
				ss(mj,l)=(ss(mjm1,l)*(r(mj)-r(mjm2))-ss(mjm2,l)*(r(mj)-r(mjm1)))/(r(mjm1)-r(mjm2))
			end do
!$OMP END PARALLEL DO
		else if (myPE == 0) then
!$OMP PARALLEL DO
			do l=1,lmax
				ss(0,l)=0.0_IDP
				if (mm(l) == 0) ss(0,l)=(r(2)**2*ss(1,l)-r(1)**2*ss(2,l))/(r(2)**2-r(1)**2)
			end do
!$OMP END PARALLEL DO
		end if

		if (myPE < numPEsm1) then
			do l=1,lmax
				ssend(l)=ss(mj_end-1,l)
			end do
			tag=myPE+1
			call MPI_SEND(ssend,lmax,MPI_DOUBLE_PRECISION,myPE+1,tag,MPI_COMM_WORLD,ierr)
		end if
		if (myPE > 0) then
			tag=myPE
			call MPI_RECV(srecv,lmax,MPI_DOUBLE_PRECISION,myPE-1,tag,MPI_COMM_WORLD,status,ierr)
			do l=1,lmax
				ss(mj_start,l)=srecv(l)
			end do
		end if
		if (myPE > 0) then
			do l=1,lmax
				ssend(l)=ss(mj_start+1,l)
			end do
			tag=numPEs+myPE
			call MPI_SEND(ssend,lmax,MPI_DOUBLE_PRECISION,myPE-1,tag,MPI_COMM_WORLD,ierr)
		end if
		if (myPE < numPEsm1) then
			tag=numPEs+myPE+1
			call MPI_RECV(srecv,lmax,MPI_DOUBLE_PRECISION,myPE+1,tag,MPI_COMM_WORLD,status,ierr)
			do l=1,lmax
				ss(mj_end,l)=srecv(l)
			end do
		end if

		wk1(mj_start:mj_end,0)=0. 
		wk2(mj_start:mj_end,0)=0. 
		wk3(mj_start:mj_end,0)=0. 
		ss(mj_start:mj_end,0)=0. 
		ff(mj_start:mj_end,0)=0. 

	end subroutine delstar

	subroutine dlstar(ss,ff,itypf,wk1,wk2,wk3,c1,c2)

		use mpi
		use param
		use processor
		use var_para
		use domain
		use equil
		implicit none

		integer :: itypf,l,j,tag,ierr
		real(IDP) :: c1,c2,xm
		real(IDP), dimension(mj_start:,0:) :: ss,ff,wk1,wk2,wk3
		real(IDP), dimension(lmax) :: ssend,srecv
		integer, dimension(MPI_STATUS_SIZE) :: status

		interface
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
		end interface

		ff(mj_start:mj_end,0)=0. 
		ss(mj_start:mj_end,0)=0. 
		wk1=0. 
		wk2=0. 
		wk3=0. 

!                write(0,'(" mj_start = ",i8)') mj_start
!                write(0,'(" mj_end = ",i8)') mj_end
!                write(0,'(" myPE = ",i8)') myPE
!                write(0,'(" lmax = ",i8)') lmax

		if (myPE == 0) then
			do l=1,lmax
				do j=1,mj_end-1
					wk2(j,l)=del2cp(j)*(ff(j+1,l)-ff(j,l))+del2cm(j)*(ff(j-1,l)-ff(j,l))
				end do
			end do
		else if (myPE == numPEsm1) then
			do l=1,lmax
				do j=mj_start+1,mjm1
					wk2(j,l)=del2cp(j)*(ff(j+1,l)-ff(j,l))+del2cm(j)*(ff(j-1,l)-ff(j,l))
				end do
			end do
		else
			do l=1,lmax
				do j=mj_start+1,mj_end-1
					wk2(j,l)=del2cp(j)*(ff(j+1,l)-ff(j,l))+del2cm(j)*(ff(j-1,l)-ff(j,l))
				end do
			end do
		end if

!		do l=1,lmax
!                        write(0,'(" wk2 = ",1pe10.3)') wk2(0,l)
!                        write(0,'(" wk2 = ",1pe10.3)') wk2(mj_end-1,l)
!		end do

!                write(0,'(" leqmax = ",i8)') leqmax
!		do l=1,leqmax
!			j=ll(mmeq(l),nneq(l))
!                        write(0,'(" l = ",i8)') l
!                       write(0,'(" j = ",i8)') j
!                        write(0,'(" mmeq = ",i8)') mmeq(l)
!                        write(0,'(" nneq = ",i8)') nneq(l)
!		end do
  
		call eqtodyn_par(wk3,jbgtt,0.0_IDP,1.0_IDP)
		do l=1,lmax
			wk3(:,l)=denseq(mj_start:mj_end)*wk3(:,l)
		end do
		call mult(ss,wk2,itypf,wk3,2,c1,c2)
		call dbydr_par(wk2,wk3,0.0_IDP,1.0_IDP,2)
		call eqtodyn_par(wk3,jbgrt,0.0_IDP,-1.0_IDP)
		do l=1,lmax
			wk3(:,l)=denseq(mj_start:mj_end)*wk3(:,l)
		end do
		call dbydth_par(wk1,wk3,-1,0.0_IDP,1.0_IDP,2)
		do l=1,lmax
			wk1(:,l)=wk1(:,l)+wk2(:,l)
		end do
		call dbydr_par(wk2,ff,0.0_IDP,1.0_IDP,0)
		call mult(ss,wk2,itypf,wk1,2,1.0_IDP,c2)
		call dbydth_par(wk1,wk2,itypf,0.0_IDP,1.0_IDP,0)
		call mult(ss,wk1,-itypf,wk3,-2,1.0_IDP,2.*c2)
		call dbydr_par(wk1,wk3,0.0_IDP,1.0_IDP,0)
		call eqtodyn_par(wk3,jbgrr,0.0_IDP,1.0_IDP)
		do l=1,lmax
			wk3(:,l)=denseq(mj_start:mj_end)*wk3(:,l)
		end do
		call dbydth_par(wk2,wk3,1,0.0_IDP,1.0_IDP,2)
		do l=1,lmax
			wk1(:,l)=wk1(:,l)+wk2(:,l)
		end do
		call dbydth_par(wk2,ff,itypf,0.0_IDP,1.0_IDP,0)
		call mult(ss,wk2,-itypf,wk1,-2,1.0_IDP,c2)
		if (myPE == 0) then
			do l=1,lmax
				xm=mm(l)
				wk2(1:mj_end,l)=-(rinv(1:mj_end)*xm)**2*ff(1:mj_end,l)
			end do
		else
			do l=1,lmax
				xm=mm(l)
				wk2(mj_start:mj_end,l)=-(rinv(mj_start:mj_end)*xm)**2*ff(mj_start:mj_end,l)
			end do
		end if
		call mult(ss,wk3,2,wk2,itypf,1.0_IDP,c2)

		if (myPE == numPEsm1) then
			do l=1,lmax
				ss(mj,l)=(ss(mjm1,l)*(r(mj)-r(mjm2))-ss(mjm2,l)*(r(mj)-r(mjm1)))/(r(mjm1)-r(mjm2))
			end do
		else if (myPE == 0) then
			do l=1,lmax
				ss(0,l)=0.0_IDP
				if (mm(l) == 0) ss(0,l)=(r(2)**2*ss(1,l)-r(1)**2*ss(2,l))/(r(2)**2-r(1)**2)
			end do
		end if

		if (myPE < numPEsm1) then
			do l=1,lmax
				ssend(l)=ss(mj_end-1,l)
			end do
			tag=myPE+1
			call MPI_SEND(ssend,lmax,MPI_DOUBLE_PRECISION,myPE+1,tag,MPI_COMM_WORLD,ierr)
		end if
		if (myPE > 0) then
			tag=myPE
			call MPI_RECV(srecv,lmax,MPI_DOUBLE_PRECISION,myPE-1,tag,MPI_COMM_WORLD,status,ierr)
			do l=1,lmax
				ss(mj_start,l)=srecv(l)
			end do
		end if
		if (myPE > 0) then
			do l=1,lmax
				ssend(l)=ss(mj_start+1,l)
			end do
			tag=numPEs+myPE
			call MPI_SEND(ssend,lmax,MPI_DOUBLE_PRECISION,myPE-1,tag,MPI_COMM_WORLD,ierr)
		end if
		if (myPE < numPEsm1) then
			tag=numPEs+myPE+1
			call MPI_RECV(srecv,lmax,MPI_DOUBLE_PRECISION,myPE+1,tag,MPI_COMM_WORLD,status,ierr)
			do l=1,lmax
				ss(mj_end,l)=srecv(l)
			end do
		end if

		wk1(mj_start:mj_end,0)=0. 
		wk2(mj_start:mj_end,0)=0. 
		wk3(mj_start:mj_end,0)=0. 
		ss(mj_start:mj_end,0)=0. 
		ff(mj_start:mj_end,0)=0. 

	end subroutine dlstar

	subroutine del2c(d,a,c1,c2,k)

		use param
		use domain
		implicit none

		integer :: k,l,m,j
		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: a,d
		real(IDP), dimension(0:mj) :: dold

		a(:,0)=0. 
!$OMP PARALLEL DO PRIVATE(dold,j,m)
		do l=1,lmax
			dold=d(0:mj,l)
			m=abs(mm(l))
			do j=1,mjm1
				d(j,l)=del2cm(j)*a(j-1,l)-(del2cm(j)+del2cp(j)+(m*rinv(j))**2)*a(j,l)+del2cp(j)*a(j+1,l)
			end do
			d(mj,l)=(d(mjm1,l)*(r(mj)-r(mj-2))-d(mj-2,l)*(r(mj)-r(mjm1)))/(r(mjm1)-r(mj-2))
			d(0,l)=0.
			if (k == 1 .or. k == 2) m=abs(m-1)
			if (m == 0) d(0,l)=(a(1,l)-a(0,l))*4.*rinv(1)**2
			d(0:mj,l)=c1*dold+c2*d(0:mj,l)
		end do
!$OMP END PARALLEL DO
		d(:,0)=0. 

	end subroutine del2c

	subroutine del2cl(d,a,c1,c2,k,l)

		use param
		use domain
		implicit none

		integer :: k,l,m,j
		real(IDP) :: c1,c2
		real(IDP), dimension(0:) :: a,d
		real(IDP), dimension(0:mj) :: dold

		dold=d 
		m=abs(mm(l))
		if (k == 1 .or. k == 2) m=abs(m-1)
		do j=1,mjm1
			d(j)=del2cm(j)*a(j-1)-(del2cm(j)+del2cp(j)+(m*rinv(j))**2)*a(j)+del2cp(j)*a(j+1)
		end do
		d(mj)=(d(mjm1)*(r(mj)-r(mj-2))-d(mj-2)*(r(mj)-r(mjm1)))/(r(mjm1)-r(mj-2))
		d(0)=0.
		if (k == 1 .or. k == 2) m=abs(m-1)
		if (m == 0) d(0)=(a(1)-a(0))*4.*rinv(1)**2
		d=c1*dold+c2*d

	end subroutine del2cl