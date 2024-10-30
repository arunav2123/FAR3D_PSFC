	subroutine eqtodyn(adyn,aeq,c1,c2)

		use param
		use domain
		implicit none

		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: adyn,aeq
		integer :: le,l

		adyn=c1*adyn
		do le=1,leqmax
			l=ll(mmeq(le),nneq(le))
			if (l == 0) cycle
			adyn(:,l)=adyn(:,l)+c2*aeq(:,le) 
		end do
		adyn(:,0)=0. 
		aeq(:,0)=0. 

	end subroutine eqtodyn

	subroutine eqtodyn_par(adyn,aeq,c1,c2)

		use param
		use var_para
		use domain
		implicit none

		real(IDP) :: c1,c2
		real(IDP), dimension(mj_start:,0:) :: adyn
		real(IDP), dimension(0:,0:) :: aeq
		integer :: le,l

		do l=1,lmax
			adyn(mj_start:mj_end,l)=c1*adyn(mj_start:mj_end,l)
		end do
		do le=1,leqmax
			l=ll(mmeq(le),nneq(le))
			if (l == 0) cycle
			adyn(mj_start:mj_end,l)=adyn(mj_start:mj_end,l)+c2*aeq(mj_start:mj_end,le) 
		end do
		adyn(mj_start:mj_end,0)=0. 
		aeq(mj_start:mj_end,0)=0. 

	end subroutine eqtodyn_par

	subroutine dyntoeq(aeq,adyn,c1,c2)

		use param
		use domain
		implicit none

		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: aeq,adyn
		integer :: le,l

		aeq=c1*aeq
		do le=1,leqmax
			l=ll(mmeq(le),nneq(le))
			if (l == 0) cycle
			aeq(:,le)=aeq(:,le)+c2*adyn(:,l) 
		end do
		aeq(:,0)=0. 
		adyn(:,0)=0. 

	end subroutine dyntoeq

	subroutine rhotodyn(adyn,arho,c1,c2)

		use param
		use domain
		implicit none

		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: adyn
		real(IDP), dimension(0:) :: arho

		adyn=c1*adyn
		adyn(:,l0)=adyn(:,l0)+c2*arho
		adyn(:,0)=0. 

	end subroutine rhotodyn