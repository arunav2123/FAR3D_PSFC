	subroutine om(tx,itx,ieqn,ivar,ith,ir,izt,coef)

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		integer :: itx,ieqn,ivar,ith,ir,izt,l
		real(IDP) :: coef
		real(IDP), dimension(0:,0:) :: tx
		real(IDP), dimension(0:jdim,0:leqdim) :: wrk
		
		interface
			subroutine blockj(tx,itx,ieqn,ivar,ith,ir,izt,coef)
				use param
				implicit none
				integer :: itx,ieqn,ivar,ith,ir,izt
				real(IDP) :: coef
				real(IDP), dimension(0:,0:) :: tx
			end subroutine blockj
		end interface

		do l=1,leqmax
			wrk(1:mjm1,l)=tx(1:mjm1,l)
		end do
		call blockj(wrk,itx,ieqn,ivar,ith,ir,izt+1,coef)
		do l=1,leqmax
			wrk(1:mjm1,l)=r(1:mjm1)*wrk(1:mjm1,l)/qq(1:mjm1)
		end do
		call blockj(wrk,itx,ieqn,ivar,ith+1,ir,izt,-coef)

	end subroutine om

	subroutine omc(tx,itx,ieqn,ivar,ith,ir,izt,coef)

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		integer :: itx,ieqn,ivar,ith,ir,izt,l
		real(IDP) :: coef
		real(IDP), dimension(0:,0:) :: tx
		real(IDP), dimension(0:jdim,0:leqdim) :: wrk

		interface
			subroutine b2lx(tx,itx,ieqn,ivar,ith,ir,izt,coef)
				use param
				implicit none
				integer :: itx,ieqn,ivar,ith,ir,izt
				real(IDP) :: coef
				real(IDP), dimension(0:,0:) :: tx
			end subroutine b2lx
		end interface

		do l=1,leqmax
			wrk(1:mjm1,l)=tx(1:mjm1,l)
		end do
		call b2lx(wrk,itx,ieqn,ivar,ith,ir,izt+1,coef)
		do l=1,leqmax
			wrk(1:mjm1,l)=r(1:mjm1)*wrk(1:mjm1,l)/qq(1:mjm1)
		end do
		call b2lx(wrk,itx,ieqn,ivar,ith+1,ir,izt,-coef)

	end subroutine omc

	subroutine om2(tx,itx,ieqn,ivar,ith,ir,izt,coef)

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		interface
			subroutine blockj(tx,itx,ieqn,ivar,ith,ir,izt,coef)
				use param
				implicit none
				integer :: itx,ieqn,ivar,ith,ir,izt
				real(IDP) :: coef
				real(IDP), dimension(0:,0:) :: tx
			end subroutine blockj
		end interface

		integer :: itx,ieqn,ivar,ith,ir,izt,l
		real(IDP) :: coef
		real(IDP), dimension(0:,0:) :: tx
		real(IDP), dimension(0:jdim,0:leqdim) :: wrk

		do l=1,leqmax
			wrk(1:mjm1,l)=tx(1:mjm1,l)
		end do
		call blockj(wrk,itx,ieqn,ivar,ith,ir,izt+2,coef)
		do l=1,leqmax
			wrk(1:mjm1,l)=r(1:mjm1)*wrk(1:mjm1,l)/qq(1:mjm1)
		end do
		call blockj(wrk,itx,ieqn,ivar,ith+1,ir,izt+1,-2.*coef)
		do l=1,leqmax
			wrk(1:mjm1,l)=r(1:mjm1)*wrk(1:mjm1,l)/qq(1:mjm1)
		end do
		call blockj(wrk,itx,ieqn,ivar,ith+2,ir,izt,coef)

	end subroutine om2

	subroutine om2c(tx,itx,ieqn,ivar,ith,ir,izt,coef)

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		integer :: itx,ieqn,ivar,ith,ir,izt,l
		real(IDP) :: coef
		real(IDP), dimension(0:,0:) :: tx
		real(IDP), dimension(0:jdim,0:leqdim) :: wrk

		interface
			subroutine b2lx(tx,itx,ieqn,ivar,ith,ir,izt,coef)
				use param
				implicit none
				integer :: itx,ieqn,ivar,ith,ir,izt
				real(IDP) :: coef
				real(IDP), dimension(0:,0:) :: tx
			end subroutine b2lx
		end interface

		do l=1,leqmax
			wrk(1:mjm1,l)=tx(1:mjm1,l)
		end do
		call b2lx(wrk,itx,ieqn,ivar,ith,ir,izt+2,coef)
		do l=1,leqmax
			wrk(1:mjm1,l)=r(1:mjm1)*wrk(1:mjm1,l)/qq(1:mjm1)
		end do
		call b2lx(wrk,itx,ieqn,ivar,ith+1,ir,izt+1,-2.*coef)
		do l=1,leqmax
			wrk(1:mjm1,l)=r(1:mjm1)*wrk(1:mjm1,l)/qq(1:mjm1)
		end do
		call b2lx(wrk,itx,ieqn,ivar,ith+2,ir,izt,coef)

	end subroutine om2c
	
	subroutine om0(tx,ieqn,ivar,ith,ir,izt,coef)

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		integer :: ieqn,ivar,ith,ir,izt
		real(IDP) :: coef
		real(IDP), dimension(0:) :: tx
		real(IDP), dimension(0:jdim) :: wrk

		interface
			subroutine block0(tx,ieqn,ivar,ith,ir,izt,coef)
				use param
				implicit none
				integer :: ieqn,ivar,ith,ir,izt
				real(IDP) :: coef
				real(IDP), dimension(0:) :: tx
			end subroutine block0
		end interface

		wrk(1:mjm1)=tx(1:mjm1)
		call block0(wrk,ieqn,ivar,ith,ir,izt+1,coef)
		wrk(1:mjm1)=r(1:mjm1)*wrk(1:mjm1)/qq(1:mjm1)
		call block0(wrk,ieqn,ivar,ith+1,ir,izt,-coef)

	end subroutine om0

	subroutine omc0(tx,ieqn,ivar,ith,ir,izt,coef)

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		integer :: ieqn,ivar,ith,ir,izt
		real(IDP) :: coef
		real(IDP), dimension(0:) :: tx
		real(IDP), dimension(0:jdim) :: wrk

		interface
			subroutine b2lx0(tx,ieqn,ivar,ith,ir,izt,coef)
				use param
				implicit none
				integer :: ieqn,ivar,ith,ir,izt
				real(IDP) :: coef
				real(IDP), dimension(0:) :: tx
			end subroutine b2lx0
		end interface

		wrk(1:mjm1)=tx(1:mjm1)
		call b2lx0(wrk,ieqn,ivar,ith,ir,izt+1,coef)
		wrk(1:mjm1)=r(1:mjm1)*wrk(1:mjm1)/qq(1:mjm1)
		call b2lx0(wrk,ieqn,ivar,ith+1,ir,izt,-coef)

	end subroutine omc0