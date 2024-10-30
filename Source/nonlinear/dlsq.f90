	subroutine dlsq(ss,ff,itypf,wk1,wk2,wk3,c1,c2)

		use mpi
		use param
		use processor
		use var_para
		use domain
		use equil
		implicit none

		integer :: itypf,l,j,le,tag,ierr
		real(IDP) :: c1,c2,xm,xn,coefs
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
			subroutine dbydzt_par(d,a,ltype,c1,c2)
				use param
				use var_para
				implicit none
				integer :: ltype
				real(IDP) :: c1,c2
				real(IDP), dimension(mj_start:,0:) :: a,d
			end subroutine dbydzt_par
		end interface

		ff(mj_start:mj_end,0)=0. 
		ss(mj_start:mj_end,0)=0. 
		wk1=0. 
		wk2=0. 
		wk3=0. 

		if (myPE == 0) then
			do l=1,lmax
				do j=1,mj_end-1
					wk2(j,l)=dc2p(j)*(ff(j+1,l)-ff(j,l))+dc2m(j)*(ff(j-1,l)-ff(j,l))
				end do
			end do
		else if (myPE == numPEsm1) then
			do l=1,lmax
				do j=mj_start+1,mjm1
					wk2(j,l)=dc2p(j)*(ff(j+1,l)-ff(j,l))+dc2m(j)*(ff(j-1,l)-ff(j,l))
				end do
			end do
		else
			do l=1,lmax
				do j=mj_start+1,mj_end-1
					wk2(j,l)=dc2p(j)*(ff(j+1,l)-ff(j,l))+dc2m(j)*(ff(j-1,l)-ff(j,l))
				end do
			end do
		end if
		call eqtodyn_par(wk3,lplrr,0.0_IDP,1.0_IDP)
		call mult(ss,wk2,itypf,wk3,2,c1,c2)
		call eqtodyn_par(wk1,lplr,0.0_IDP,1.0_IDP)
		call dbydr_par(wk2,ff,0.0_IDP,1.0_IDP,0)
		call mult(ss,wk2,itypf,wk1,2,1.0_IDP,c2)
		call eqtodyn_par(wk3,lplrt,0.0_IDP,1.0_IDP)
		call dbydth_par(wk1,wk2,itypf,0.0_IDP,1.0_IDP,0)
		call mult(ss,wk1,-itypf,wk3,-2,1.0_IDP,c2)
		call eqtodyn_par(wk3,lplrz,0.0_IDP,1.0_IDP)
		call dbydzt_par(wk1,wk2,itypf,0.0_IDP,1.0_IDP)
		call mult(ss,wk1,-itypf,wk3,-2,1.0_IDP,c2)
		call eqtodyn_par(wk1,lplt,0.0_IDP,1.0_IDP)
		call dbydth_par(wk2,ff,itypf,0.0_IDP,1.0_IDP,0)
		call mult(ss,wk2,-itypf,wk1,-2,1.0_IDP,c2)
		call eqtodyn_par(wk3,lpltz,0.0_IDP,1.0_IDP)
		call dbydzt_par(wk1,wk2,-itypf,0.0_IDP,1.0_IDP)
		call mult(ss,wk1,itypf,wk3,2,1.0_IDP,c2)
		call eqtodyn_par(wk3,lpltt,0.0_IDP,1.0_IDP)
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
		call mult(ss,wk2,itypf,wk3,2,1.0_IDP,c2)
		call eqtodyn_par(wk3,lplzz,0.0_IDP,1.0_IDP)
		do l=1,lmax
			xn=nn(l)
			wk2(mj_start:mj_end,l)=-xn*xn*ff(mj_start:mj_end,l)
		end do
		call mult(ss,wk2,itypf,wk3,2,1.0_IDP,c2)
		call eqtodyn_par(wk1,lplz,0.0_IDP,1.0_IDP)
		call dbydzt_par(wk2,ff,itypf,0.0_IDP,1.0_IDP)
		call mult(ss,wk2,-itypf,wk1,-2,1.0_IDP,c2)

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

	end subroutine dlsq