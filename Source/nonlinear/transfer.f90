	subroutine trnsfr0(f,ftype)

		use mpi
		use param
		use processor
		use var_para
		use domain
		implicit none

		integer :: ftype
		real(IDP), dimension(mj_start:,0:) :: f
		integer :: i,j,l,l1,l1t,lp,ln,iPE,iPE1,ierr,tag
		integer, dimension(MPI_STATUS_SIZE) :: status

		if (myPE > 0) then

			do iPE=0,numPEsm1
				l1=0
				do i=nr_st(iPE),nr_nd(iPE)
					do l1t=1,mnumn(i)
						l=l1t+lnumn(i-1)
						lp=lln(l)
						l1=l1+1
						do j=mj_br(myPE),mj_br(myPE)+mj_inc(myPE)-1
							scu((j-mj_br(myPE))*lmaxPE(iPE)+l1,iPE)=f(j,lp)
						end do
					end do
				end do
				do i=nc_st(iPE),nc_nd(iPE)
					do l1t=1,mnumn(i)
						l=l1t+lnumn(i-1)
						if (ftype > 0) then
							lp=lln(l)
							ln=lln(l+mnumn(i))
						else
							ln=lln(l)
							lp=lln(l+mnumn(i))
						end if
						l1=l1+1
						do j=mj_br(myPE),mj_br(myPE)+mj_inc(myPE)-1
							scu((j-mj_br(myPE))*lmaxPE(iPE)+l1,iPE)=f(j,lp)
						end do
						l1=l1+1
						do j=mj_br(myPE),mj_br(myPE)+mj_inc(myPE)-1
							scu((j-mj_br(myPE))*lmaxPE(iPE)+l1,iPE)=f(j,ln)
						end do
					end do
				end do
			end do

!	scu(:,iPE), iPE=0,numPEsm1, contains f(j1:j2,l1:l2), where j1:j2 are the j-indexes of local process myPE and l1:l2 are the 
!	l-indexes of process iPE. Each process (except 0) sends scu(:,iPE), iPE=0,numPEsm1, to process 0

			do iPE=0,numPEsm1
				tag=numPEs*iPE+myPE
				call MPI_SEND(scu(1,iPE),mj_inc(myPE)*lmaxPE(iPE),MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,ierr)
			end do
					
		else

!	Process 0 receives scu(:,iPE), iPE=0,numPEsm1, from the other processes, and get the whole f(0:mj,1:lmax)

			do iPE=1,numPEsm1
				do iPE1=0,numPEsm1
					tag=numPEs*iPE1+iPE
					call MPI_RECV(scu(1,iPE1),mj_inc(iPE)*lmaxPE(iPE1),MPI_DOUBLE_PRECISION,iPE,tag,MPI_COMM_WORLD,status,ierr)
					l1=0
					do i=nr_st(iPE1),nr_nd(iPE1)
						do l1t=1,mnumn(i)
							l=l1t+lnumn(i-1)
							lp=lln(l)
							l1=l1+1
							do j=mj_br(iPE),mj_br(iPE)+mj_inc(iPE)-1
								f(j,lp)=scu((j-mj_br(iPE))*lmaxPE(iPE1)+l1,iPE1)
							end do
						end do
					end do
					do i=nc_st(iPE1),nc_nd(iPE1)
						do l1t=1,mnumn(i)
							l=l1t+lnumn(i-1)
							if (ftype > 0) then
								lp=lln(l)
								ln=lln(l+mnumn(i))
							else
								ln=lln(l)
								lp=lln(l+mnumn(i))
							end if
							l1=l1+1
							do j=mj_br(iPE),mj_br(iPE)+mj_inc(iPE)-1
								f(j,lp)=scu((j-mj_br(iPE))*lmaxPE(iPE1)+l1,iPE1)
							end do
							l1=l1+1
							do j=mj_br(iPE),mj_br(iPE)+mj_inc(iPE)-1
								f(j,ln)=scu((j-mj_br(iPE))*lmaxPE(iPE1)+l1,iPE1)
							end do
						end do
					end do
				end do
			end do
						
		end if

	end subroutine trnsfr0

	subroutine trnsfr0l(f)

		use mpi
		use param
		use processor
		use var_para
		use domain
		implicit none

		real(IDP), dimension(mj_start:,0:) :: f
		integer :: i,j,l,l1,l1t,iPE,iPE1,ierr,tag
		integer, dimension(MPI_STATUS_SIZE) :: status

		if (myPE > 0) then

			do iPE=0,numPEsm1
				l1=0
				do i=nr_st(iPE),nr_nd(iPE)
					do l1t=1,mnumn(i)
						l=l1t+lnumn(i-1)
						l1=l1+1
						do j=mj_br(myPE),mj_br(myPE)+mj_inc(myPE)-1
							scu((j-mj_br(myPE))*lmaxPE(iPE)+l1,iPE)=f(j,l)
						end do
					end do
				end do
				do i=nc_st(iPE),nc_nd(iPE)
					do l1t=1,2*mnumn(i)
						l=l1t+lnumn(i-1)
						l1=l1+1
						do j=mj_br(myPE),mj_br(myPE)+mj_inc(myPE)-1
							scu((j-mj_br(myPE))*lmaxPE(iPE)+l1,iPE)=f(j,l)
						end do
					end do
				end do
			end do

!	scu(:,iPE), iPE=0,numPEsm1, contains f(j1:j2,l1:l2), where j1:j2 are the j-indexes of local process myPE and l1:l2 are the 
!	l-indexes of process iPE. Each process (except 0) sends scu(:,iPE), iPE=0,numPEsm1, to process 0

			do iPE=0,numPEsm1
				tag=numPEs*iPE+myPE
				call MPI_SEND(scu(1,iPE),mj_inc(myPE)*lmaxPE(iPE),MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,ierr)
			end do
					
		else

!	Process 0 receives scu(:,iPE), iPE=0,numPEsm1, from the other processes, and get the whole f(0:mj,1:lmax)

			do iPE=1,numPEsm1
				do iPE1=0,numPEsm1
					tag=numPEs*iPE1+iPE
					call MPI_RECV(scu(1,iPE1),mj_inc(iPE)*lmaxPE(iPE1),MPI_DOUBLE_PRECISION,iPE,tag,MPI_COMM_WORLD,status,ierr)
					l1=0
					do i=nr_st(iPE1),nr_nd(iPE1)
						do l1t=1,mnumn(i)
							l=l1t+lnumn(i-1)
							l1=l1+1
							do j=mj_br(iPE),mj_br(iPE)+mj_inc(iPE)-1
								f(j,l)=scu((j-mj_br(iPE))*lmaxPE(iPE1)+l1,iPE1)
							end do
						end do
					end do
					do i=nc_st(iPE1),nc_nd(iPE1)
						do l1t=1,2*mnumn(i)
							l=l1t+lnumn(i-1)
							l1=l1+1
							do j=mj_br(iPE),mj_br(iPE)+mj_inc(iPE)-1
								f(j,l)=scu((j-mj_br(iPE))*lmaxPE(iPE1)+l1,iPE1)
							end do
						end do
					end do
				end do
			end do
						
		end if

	end subroutine trnsfr0l

	subroutine trnsfr0e(f,ftype)

		use mpi
		use param
		use processor
		use var_para
		use domain
		use scratch
		implicit none

		integer :: ftype
		real(IDP), dimension(mj_start:,0:) :: f
		integer :: l,le,iPE,ierr,tag
		integer, dimension(MPI_STATUS_SIZE) :: status

		interface
			subroutine trnsfr0(f,ftype)
				use param
				use var_para
				implicit none
				integer :: ftype
				real(IDP), dimension(mj_start:,0:) :: f
			end subroutine trnsfr0
		end interface

		call trnsfr0(f,ftype)

		do le=m0dy+1,leqmax

			l=ll(mmeq(le),nneq(le))
			if (l == 0) cycle

			if (myPE > 0) then
				sd1(mj_br(myPE):mj_br(myPE)+mj_inc(myPE)-1)=f(mj_br(myPE):mj_br(myPE)+mj_inc(myPE)-1,l)
				tag=le*numPEs+myPE
				call MPI_SEND(sd1(mj_br(myPE)),mj_inc(myPE),MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,ierr)
			else
				do iPE=1,numPEsm1
					tag=le*numPEs+iPE
					call MPI_RECV(sd1(mj_br(iPE)),mj_inc(iPE),MPI_DOUBLE_PRECISION,iPE,tag,MPI_COMM_WORLD,status,ierr)
				end do
				f(mj_br(1):mj,l)=sd1(mj_br(1):mj)
			end if

		end do

	end subroutine trnsfr0e

	subroutine trnsfr(f,ftype,idir)

		use mpi
		use param
		use processor
		use var_para
		use domain
		implicit none

		integer :: ftype,idir
		real(IDP), dimension(mj_start:,0:) :: f
		integer :: i,j,l,l1,l1t,lp,ln,dPE,nPE1,nPE2,ierr,tag1,tag2
		integer, dimension(MPI_STATUS_SIZE) :: status

		interface
			subroutine pck(f,ftype)
				use param
				use var_para
				implicit none
				integer :: ftype
				real(IDP), dimension(mj_start:,0:) :: f
			end subroutine pck
			subroutine unpck(f,ftype)
				use param
				use var_para
				implicit none
				integer :: ftype
				real(IDP), dimension(mj_start:,0:) :: f
			end subroutine unpck
		end interface

		select case (idir)
			case (1)
				l1=0
				do i=nr_start,nr_end
					do l1t=1,mnumn(i)
						l=l1t+lnumn(i-1)
						lp=lln(l)
						l1=l1+1
						do j=mj_start,mj_end
							f(j,lp)=scp(l1,j)
						end do
					end do
				end do
				do i=nc_start,nc_end
					do l1t=1,mnumn(i)
						l=l1t+lnumn(i-1)
						if (ftype > 0) then
							lp=lln(l)
							ln=lln(l+mnumn(i))
						else
							ln=lln(l)
							lp=lln(l+mnumn(i))
						end if
						l1=l1+1
						do j=mj_start,mj_end
							f(j,lp)=scp(l1,j)
						end do
						do j=mj_start,mj_end
							f(j,ln)=scp(l1+mnumn(i),j)
						end do
					end do
				end do

!	Each process sends vector(l1:l2,j1:j2) to other processes, where l1:l2 are the l-indexes of the process sending the vector
!	and j1:j2 are the j-indexes of the process receiving the vector. After communication is completed, each process has all the 
!	variable(j1:j2,1:lmax), where j1:j2 are the j-indexes of the process.

				do dPE=1,numPEsm1
					nPE1=myPE+dPE
					if (nPE1 > numPEsm1) nPE1=nPE1-numPEs
					tag1=numPEs*nPE1+myPE
					nPE2=myPE-dPE
					if (nPE2 < 0) nPE2=nPE2+numPEs
					tag2=numPEs*myPE+nPE2
					call MPI_SENDRECV(scp(1,mj_st(nPE1)),mj_dl(nPE1)*lmaxPE(myPE),MPI_DOUBLE_PRECISION,nPE1,tag1, &
						scu(1,nPE2),mj_dl(myPE)*lmaxPE(nPE2),MPI_DOUBLE_PRECISION,nPE2,tag2,MPI_COMM_WORLD,status,ierr)
				end do
				call unpck(f,ftype)
			case (2)
				l1=0
				do i=nr_start,nr_end
					do l1t=1,mnumn(i)
						l=l1t+lnumn(i-1)
						lp=lln(l)
						l1=l1+1
						do j=mj_br(myPE),mj_br(myPE)+mj_inc(myPE)-1
							scp(l1,j)=f(j,lp)
						end do
					end do
				end do 
                              
				do i=nc_start,nc_end
					do l1t=1,mnumn(i)
						l=l1t+lnumn(i-1)
						if (ftype > 0) then
							lp=lln(l)
							ln=lln(l+mnumn(i))
						else
							ln=lln(l)
							lp=lln(l+mnumn(i))
						end if
						l1=l1+1
						do j=mj_br(myPE),mj_br(myPE)+mj_inc(myPE)-1
							scp(l1,j)=f(j,lp)
						end do
						do j=mj_br(myPE),mj_br(myPE)+mj_inc(myPE)-1
							scp(l1+mnumn(i),j)=f(j,ln)
						end do

					end do
				end do

				call pck(f,ftype)				

!	Each process sends variable(l1:l2,j1:j2) to other processes, where j1:j2 are the j-indexes of the process sending the vector
!	and l1:l2 are the l-indexes of the process receiving the vector. After communication is completed, each process has all the 
!	vector(l1:l2,1:mj), where l1:l2 are the l-indexes of the process.

				do dPE=1,numPEsm1
					nPE1=myPE+dPE
					if (nPE1 > numPEsm1) nPE1=nPE1-numPEs
					tag1=numPEs*nPE1+myPE
					nPE2=myPE-dPE
					if (nPE2 < 0) nPE2=nPE2+numPEs
					tag2=numPEs*myPE+nPE2
					call MPI_SENDRECV(scu(1,nPE1),mj_inc(myPE)*lmaxPE(nPE1),MPI_DOUBLE_PRECISION,nPE1,tag1, &
						scp(1,mj_br(nPE2)),mj_inc(nPE2)*lmaxPE(myPE),MPI_DOUBLE_PRECISION,nPE2,tag2,MPI_COMM_WORLD,status,ierr)
				end do
		end select

	end subroutine trnsfr

	subroutine pck(f,ftype)

		use param
		use processor
		use var_para
		use domain
		implicit none

		integer :: ftype,iPE,i,j,l,l1,l1t,lp,ln
		real(IDP), dimension(mj_start:,0:) :: f

		do iPE=0,numPEsm1
			if (iPE /= myPE) then
				l1=0
				do i=nr_st(iPE),nr_nd(iPE)
					do l1t=1,mnumn(i)
						l=l1t+lnumn(i-1)
						lp=lln(l)
						l1=l1+1
						do j=mj_br(myPE),mj_br(myPE)+mj_inc(myPE)-1
							scu((j-mj_br(myPE))*lmaxPE(iPE)+l1,iPE)=f(j,lp)
						end do
					end do
				end do
				do i=nc_st(iPE),nc_nd(iPE)
					do l1t=1,mnumn(i)
						l=l1t+lnumn(i-1)
						if (ftype > 0) then
							lp=lln(l)
							ln=lln(l+mnumn(i))
						else
							ln=lln(l)
							lp=lln(l+mnumn(i))
						end if
						l1=l1+1
						do j=mj_br(myPE),mj_br(myPE)+mj_inc(myPE)-1
							scu((j-mj_br(myPE))*lmaxPE(iPE)+l1,iPE)=f(j,lp)
						end do
						do j=mj_br(myPE),mj_br(myPE)+mj_inc(myPE)-1
							scu((j-mj_br(myPE))*lmaxPE(iPE)+l1+mnumn(i),iPE)=f(j,ln)
						end do
					end do
				end do
			end if
		end do

	end subroutine pck

	subroutine unpck(f,ftype)

		use param
		use processor
		use var_para
		use domain
		implicit none

		integer :: ftype,iPE,i,j,l,l1,l1t,lp,ln
		real(IDP), dimension(mj_start:,0:) :: f

		do iPE=0,numPEsm1
			if (iPE /= myPE) then
				l1=0
				do i=nr_st(iPE),nr_nd(iPE)
					do l1t=1,mnumn(i)
						l=l1t+lnumn(i-1)
						lp=lln(l)
						l1=l1+1
						do j=mj_start,mj_end
							f(j,lp)=scu((j-mj_start)*lmaxPE(iPE)+l1,iPE)
						end do
					end do
				end do
				do i=nc_st(iPE),nc_nd(iPE)
					do l1t=1,mnumn(i)
						l=l1t+lnumn(i-1)
						if (ftype > 0) then
							lp=lln(l)
							ln=lln(l+mnumn(i))
						else
							ln=lln(l)
							lp=lln(l+mnumn(i))
						end if
						l1=l1+1
						do j=mj_start,mj_end
							f(j,lp)=scu((j-mj_start)*lmaxPE(iPE)+l1,iPE)
						end do
						do j=mj_start,mj_end
							f(j,ln)=scu((j-mj_start)*lmaxPE(iPE)+l1+mnumn(i),iPE)
						end do
					end do
				end do
			end if
		end do

	end subroutine unpck