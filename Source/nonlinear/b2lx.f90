	subroutine b2lx(tx,itx,ieqn,ivar,ith,ir,izt,coef)

		use param
		use var_para
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		integer :: itx,ieqn,ivar,ith,ir,izt
		real(IDP) :: coef,x
		real(IDP), dimension(0:,0:) :: tx
		real(IDP), dimension(jdim) :: xa
		integer :: ity,isv,ips,iph,itz,icn,icg,ilp,ich,i,mnum2,mnum3,l,m,n,l1,l2,l3,l1t,l2t,l3n,l1p,l2p,lp,lpp,lpn,ibnd,j,imt1,imt2,il, &
			   lpq,lq,lq2,imtq

		ity= 1
		if (ivar == 2 .or. ivar == 4 .or. ivar == 6 .or. ivar == 7 .or. ivar == ivalp .or. ivar > iq) ity=-1
		isv=1
		ips=(-1)**((ith+izt+1)/2)
		iph=ity**(ith+izt)*ips
		itz=(-1)**(ith+izt)
		ity=ity*itz
		icn=itx*itz
		ilp=(1-icn)/2
		ich=itx*ity
		if (ieqn == 2 .or. ieqn == 4 .or. ieqn == 6 .or. ieqn == 7 .or. ieqn == ivalp .or. ieqn > iq) ich=-ich
		do i=nr_start,nr_end
			mnum3=noeqn*mnumn(i)
			do l2t=1,mnumn(i)
				l=l2t+lnumn(i-1)
				l2=l2t+(ivar-1)*mnumn(i)
				m=mm(lln(l))
				n=nn(lln(l))
				ibnd=1
				if (isv == 1 .and. m == 0) ibnd=2
				do l1t=1,mnumn(i)
					lp=l1t+lnumn(i-1)
					lpp=lp
					if (ich == -1) lpp=lo(lp)
					if (lpp == 0) cycle
					l1=lpp+(ieqn-1)*mnumn(i)-lnumn(i-1)
					lpn=nskp2n(i)+mnumn(i)*(l2t-1)+l1t
					do l3=1,leqmax
						l3n=l3
						if (ity ==  1 .and. itx ==  1) x=coef*cmapp(lpn,l3n)
						if (ity ==  1 .and. itx == -1) x=coef*cmamp(lpn,l3n)
						if (ity == -1 .and. itx == -1) x=coef*cmamm(lpn,l3n)
						if (ity == -1 .and. itx ==  1) x=coef*cmapm(lpn,l3n)
						do il=1,ith
							x=x*m
						end do
						do il=1,izt
							x=x*n
						end do
						if (x == 0.) cycle
						x=x*iph
						do j=1,mjm1
							xa(j)=x*(rinv(j))**ith
						end do
						if (ir == 0) then
							do j=1,mjm1
								imt1=l1+mnum3*(j-1)+nskpxn(i)
								imt2=l2+mnum3*(j-1)+nskpxn(i)
								yt(imt1)=yt(imt1)+xa(j)*tx(j,l3n)*xt(imt2)
							end do
						else if (ir == 1) then
							imt1=l1+nskpxn(i)
							imt2=l2+nskpxn(i)
							yt(imt1)=yt(imt1)+wt10(1,ibnd)*xa(1)*tx(1,l3n)*xt(imt2)
							yt(imt1)=yt(imt1)+wt1p(1,ibnd)*xa(1)*tx(1,l3n)*xt(imt2+mnum3)
							do j=2,mjm2
								imt1=l1+mnum3*(j-1)+nskpxn(i)
								imt2=l2+mnum3*(j-1)+nskpxn(i)
								yt(imt1)=yt(imt1)+wt10(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2)
								yt(imt1)=yt(imt1)+wt1m(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2-mnum3)
								yt(imt1)=yt(imt1)+wt1p(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2+mnum3)
							end do
							imt1=l1+mnum3*mjm2+nskpxn(i)
							imt2=l2+mnum3*mjm2+nskpxn(i)
							if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
								yt(imt1)=yt(imt1)+wt10(mjm1,2)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2)
								yt(imt1)=yt(imt1)+wt1m(mjm1,2)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2-mnum3)
!							else if (ivar == 4) then
!								yt(imt1)=yt(imt1)+wt10(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2)
!								yt(imt1)=yt(imt1)+wt1m(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2-mnum3)
!								yt(imt1)=yt(imt1)+dc1p(mjm1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2)*(r(mj)-r(mjm2))/ &
!										  (r(mjm1)-r(mjm2))
!								yt(imt1)=yt(imt1)-dc1p(mjm1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2-mnum3)*(r(mj)-r(mjm1))/ &
!										  (r(mjm1)-r(mjm2))
							else
								yt(imt1)=yt(imt1)+wt10(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2)
								yt(imt1)=yt(imt1)+wt1m(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2-mnum3)
							end if
						else
							imt1=l1+nskpxn(i)
							imt2=l2+nskpxn(i)
							yt(imt1)=yt(imt1)+wt20(1,ibnd)*xa(1)*tx(1,l3n)*xt(imt2)
							yt(imt1)=yt(imt1)+wt2p(1,ibnd)*xa(1)*tx(1,l3n)*xt(imt2+mnum3)
							do j=2,mjm2
								imt1=l1+mnum3*(j-1)+nskpxn(i)
								imt2=l2+mnum3*(j-1)+nskpxn(i)
								yt(imt1)=yt(imt1)+wt20(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2)
								yt(imt1)=yt(imt1)+wt2m(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2-mnum3)
								yt(imt1)=yt(imt1)+wt2p(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2+mnum3)
							end do
							imt1=l1+mnum3*mjm2+nskpxn(i)
							imt2=l2+mnum3*mjm2+nskpxn(i)
							if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
								yt(imt1)=yt(imt1)+wt20(mjm1,2)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2)
								yt(imt1)=yt(imt1)+wt2m(mjm1,2)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2-mnum3)
!							else if (ivar == 4) then
!								yt(imt1)=yt(imt1)+wt20(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2)
!								yt(imt1)=yt(imt1)+wt2m(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2-mnum3)
!								yt(imt1)=yt(imt1)+dc2p(mjm1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2)*(r(mj)-r(mjm2))/ &
!										  (r(mjm1)-r(mjm2))
!								yt(imt1)=yt(imt1)-dc2p(mjm1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2-mnum3)*(r(mj)-r(mjm1))/ &
!										  (r(mjm1)-r(mjm2))
							else
								yt(imt1)=yt(imt1)+wt20(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2)
								yt(imt1)=yt(imt1)+wt2m(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2-mnum3)
							end if
						end if
					end do
				end do
			end do
		end do
		do i=nc_start,nc_end
			mnum2=noeqn*mnumn(i)
			mnum3=2*noeqn*mnumn(i)
			do l2t=1,mnumn(i)
				l=l2t+lnumn(i-1)
				l2=l2t+(ivar-1)*mnumn(i)
				l2p=l2+mnum2
				m=mm(lln(l))
				n=nn(lln(l))
				ibnd=1
				if (isv == 1 .and. m == 0) ibnd=2
				do l1t=1,mnumn(i)
					l1=l1t+(ieqn-1)*mnumn(i)+ilp*mnum2
					l1p=l1t+(ieqn-1)*mnumn(i)+(1-ilp)*mnum2
					lpn=nskp2n(i)+mnumn(i)*(l2t-1)+l1t
					do l3=1,leqmax
						l3n=l3
						icg=icn
						if (itz ==  1 .and. itx ==  1) x=coef*cmapp(lpn,l3n)
						if (itz ==  1 .and. itx == -1) x=coef*cmamp(lpn,l3n)
						if (itz == -1 .and. itx == -1) x=coef*cmamm(lpn,l3n)
						if (itz == -1 .and. itx ==  1) x=coef*cmapm(lpn,l3n)
						if (l3 == leq0 .and. itx == -1 .and. l1t == l2t) then
							icg=-icn
							l1=l1t+(ieqn-1)*mnumn(i)+(1-ilp)*mnum2
							l1p=l1t+(ieqn-1)*mnumn(i)+ilp*mnum2
							x=coef
						end if
						do il=1,ith
							x=x*m
						end do
						do il=1,izt
							x=x*n
						end do
						if (x == 0.) cycle
						x=x*ips
						do j=1,mjm1
							xa(j)=x*(rinv(j))**ith
						end do
						if (ir == 0) then
							do j=1,mjm1
								imt1=l1+mnum3*(j-1)+nskpxn(i)
								imt2=l2+mnum3*(j-1)+nskpxn(i)
								yt(imt1)=yt(imt1)+xa(j)*tx(j,l3n)*xt(imt2)
							end do
							do j=1,mjm1
								imt1=l1p+mnum3*(j-1)+nskpxn(i)
								imt2=l2p+mnum3*(j-1)+nskpxn(i)
								yt(imt1)=yt(imt1)+icg*xa(j)*tx(j,l3n)*xt(imt2)
							end do
						else if (ir == 1) then
							imt1=l1+nskpxn(i)
							imt2=l2+nskpxn(i)
							yt(imt1)=yt(imt1)+wt10(1,ibnd)*xa(1)*tx(1,l3n)*xt(imt2)
							yt(imt1)=yt(imt1)+wt1p(1,ibnd)*xa(1)*tx(1,l3n)*xt(imt2+mnum3)
							do j=2,mjm2
								imt1=l1+mnum3*(j-1)+nskpxn(i)
								imt2=l2+mnum3*(j-1)+nskpxn(i)
								yt(imt1)=yt(imt1)+wt10(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2)
								yt(imt1)=yt(imt1)+wt1m(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2-mnum3)
								yt(imt1)=yt(imt1)+wt1p(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2+mnum3)
							end do
							imt1=l1+mnum3*mjm2+nskpxn(i)
							imt2=l2+mnum3*mjm2+nskpxn(i)
							if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
								yt(imt1)=yt(imt1)+wt10(mjm1,2)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2)
								yt(imt1)=yt(imt1)+wt1m(mjm1,2)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2-mnum3)
!							else if (ivar == 4) then
!								yt(imt1)=yt(imt1)+wt10(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2)
!								yt(imt1)=yt(imt1)+wt1m(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2-mnum3)
!								yt(imt1)=yt(imt1)+dc1p(mjm1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2)*(r(mj)-r(mjm2))/ &
!										  (r(mjm1)-r(mjm2))
!								yt(imt1)=yt(imt1)-dc1p(mjm1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2-mnum3)*(r(mj)-r(mjm1))/ &
!										  (r(mjm1)-r(mjm2))
							else
								yt(imt1)=yt(imt1)+wt10(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2)
								yt(imt1)=yt(imt1)+wt1m(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2-mnum3)
							end if
							imt1=l1p+nskpxn(i)
							imt2=l2p+nskpxn(i)
							yt(imt1)=yt(imt1)+icg*wt10(1,ibnd)*xa(1)*tx(1,l3n)*xt(imt2)
							yt(imt1)=yt(imt1)+icg*wt1p(1,ibnd)*xa(1)*tx(1,l3n)*xt(imt2+mnum3)
							do j=2,mjm2
								imt1=l1p+mnum3*(j-1)+nskpxn(i)
								imt2=l2p+mnum3*(j-1)+nskpxn(i)
								yt(imt1)=yt(imt1)+icg*wt10(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2)
								yt(imt1)=yt(imt1)+icg*wt1m(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2-mnum3)
								yt(imt1)=yt(imt1)+icg*wt1p(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2+mnum3)
							end do
							imt1=l1p+mnum3*mjm2+nskpxn(i)
							imt2=l2p+mnum3*mjm2+nskpxn(i)
							if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
								yt(imt1)=yt(imt1)+icg*wt10(mjm1,2)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2)
								yt(imt1)=yt(imt1)+icg*wt1m(mjm1,2)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2-mnum3)
!							else if (ivar == 4) then
!								yt(imt1)=yt(imt1)+icg*wt10(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2)
!								yt(imt1)=yt(imt1)+icg*wt1m(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2-mnum3)
!								yt(imt1)=yt(imt1)+icg*dc1p(mjm1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2)*(r(mj)-r(mjm2))/ &
!										  (r(mjm1)-r(mjm2))
!								yt(imt1)=yt(imt1)-icg*dc1p(mjm1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2-mnum3)*(r(mj)-r(mjm1))/ &
!										  (r(mjm1)-r(mjm2))
							else
								yt(imt1)=yt(imt1)+icg*wt10(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2)
								yt(imt1)=yt(imt1)+icg*wt1m(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2-mnum3)
							end if
						else
							imt1=l1+nskpxn(i)
							imt2=l2+nskpxn(i)
							yt(imt1)=yt(imt1)+wt20(1,ibnd)*xa(1)*tx(1,l3n)*xt(imt2)
							yt(imt1)=yt(imt1)+wt2p(1,ibnd)*xa(1)*tx(1,l3n)*xt(imt2+mnum3)
							do j=2,mjm2
								imt1=l1+mnum3*(j-1)+nskpxn(i)
								imt2=l2+mnum3*(j-1)+nskpxn(i)
								yt(imt1)=yt(imt1)+wt20(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2)
								yt(imt1)=yt(imt1)+wt2m(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2-mnum3)
								yt(imt1)=yt(imt1)+wt2p(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2+mnum3)
							end do
							imt1=l1+mnum3*mjm2+nskpxn(i)
							imt2=l2+mnum3*mjm2+nskpxn(i)
							if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
								yt(imt1)=yt(imt1)+wt20(mjm1,2)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2)
								yt(imt1)=yt(imt1)+wt2m(mjm1,2)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2-mnum3)
!							else if (ivar == 4) then
!								yt(imt1)=yt(imt1)+wt20(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2)
!								yt(imt1)=yt(imt1)+wt2m(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2-mnum3)
!								yt(imt1)=yt(imt1)+dc2p(mjm1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2)*(r(mj)-r(mjm2))/ &
!										  (r(mjm1)-r(mjm2))
!								yt(imt1)=yt(imt1)-dc2p(mjm1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2-mnum3)*(r(mj)-r(mjm1))/ &
!										  (r(mjm1)-r(mjm2))
							else
								yt(imt1)=yt(imt1)+wt20(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2)
								yt(imt1)=yt(imt1)+wt2m(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2-mnum3)
							end if
							imt1=l1p+nskpxn(i)
							imt2=l2p+nskpxn(i)
							yt(imt1)=yt(imt1)+icg*wt20(1,ibnd)*xa(1)*tx(1,l3n)*xt(imt2)
							yt(imt1)=yt(imt1)+icg*wt2p(1,ibnd)*xa(1)*tx(1,l3n)*xt(imt2+mnum3)
							do j=2,mjm2
								imt1=l1p+mnum3*(j-1)+nskpxn(i)
								imt2=l2p+mnum3*(j-1)+nskpxn(i)
								yt(imt1)=yt(imt1)+icg*wt20(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2)
								yt(imt1)=yt(imt1)+icg*wt2m(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2-mnum3)
								yt(imt1)=yt(imt1)+icg*wt2p(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2+mnum3)
							end do
							imt1=l1p+mnum3*mjm2+nskpxn(i)
							imt2=l2p+mnum3*mjm2+nskpxn(i)
							if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
								yt(imt1)=yt(imt1)+icg*wt20(mjm1,2)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2)
								yt(imt1)=yt(imt1)+icg*wt2m(mjm1,2)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2-mnum3)
!							else if (ivar == 4) then
!								yt(imt1)=yt(imt1)+icg*wt20(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2)
!								yt(imt1)=yt(imt1)+icg*wt2m(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2-mnum3)
!								yt(imt1)=yt(imt1)+icg*dc2p(mjm1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2)*(r(mj)-r(mjm2))/ &
!										  (r(mjm1)-r(mjm2))
!								yt(imt1)=yt(imt1)-icg*dc2p(mjm1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2-mnum3)*(r(mj)-r(mjm1))/ &
!										  (r(mjm1)-r(mjm2))
							else
								yt(imt1)=yt(imt1)+icg*wt20(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2)
								yt(imt1)=yt(imt1)+icg*wt2m(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(imt2-mnum3)
							end if
						end if
					end do
				end do
			end do
		end do

	end subroutine b2lx
	
	subroutine b2lx0(tx,ieqn,ivar,ith,ir,izt,coef)

		use param
		use var_para
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		integer :: ieqn,ivar,ith,ir,izt
		real(IDP) :: coef,x
		real(IDP), dimension(0:) :: tx
		real(IDP), dimension(jdim) :: xa
		integer :: ity,isv,ips,iph,itz,ilp,ich,i,mnum2,mnum3,l,m,n,l1,l2,l1p,l2p,l2t,lp,ibnd,j,imt1,imt2,il,lpq,lq,lq2,imtq

		ity= 1
		if (ivar == 2 .or. ivar == 4 .or. ivar == 6 .or. ivar == 7 .or. ivar == ivalp .or. ivar > iq) ity=-1
		isv=1
		ips=(-1)**((ith+izt+1)/2)
		iph=ity**(ith+izt)*ips
		itz=(-1)**(ith+izt)
		ity=ity*itz
		ilp=(1-itz)/2
		ich=ity
		if (ieqn == 2 .or. ieqn == 4 .or. ieqn == 6 .or. ieqn == 7 .or. ieqn == ivalp .or. ieqn > iq) ich=-ich
!$OMP PARALLEL DO PRIVATE(mnum3,l2t,l,l2,m,n,ibnd,lp,l1,x,il,j,xa,imt1,imt2)
		do i=nr_start,nr_end
			mnum3=noeqn*mnumn(i)
			do l2t=1,mnumn(i)
				l=l2t+lnumn(i-1)
				l2=l2t+(ivar-1)*mnumn(i)
				m=mm(lln(l))
				n=nn(lln(l))
				ibnd=1
				if (isv == 1 .and. m == 0) ibnd=2
				lp=l
				if (ich == -1) lp=lo(l)
				if (lp == 0) cycle
				l1=lp+(ieqn-1)*mnumn(i)-lnumn(i-1)
				x=coef
				do il=1,ith
					x=x*m
				end do
				do il=1,izt
					x=x*n
				end do
				if (x == 0.) cycle
				x=x*iph
				do j=1,mjm1
					xa(j)=x*(rinv(j))**ith
				end do
				if (ir == 0) then
					do j=1,mjm1
						imt1=l1+mnum3*(j-1)+nskpxn(i)
						imt2=l2+mnum3*(j-1)+nskpxn(i)
						yt(imt1)=yt(imt1)+xa(j)*tx(j)*xt(imt2)
					end do
				else if (ir == 1) then
					imt1=l1+nskpxn(i)
					imt2=l2+nskpxn(i)
					yt(imt1)=yt(imt1)+wt10(1,ibnd)*xa(1)*tx(1)*xt(imt2)
					yt(imt1)=yt(imt1)+wt1p(1,ibnd)*xa(1)*tx(1)*xt(imt2+mnum3)
					do j=2,mjm2
						imt1=l1+mnum3*(j-1)+nskpxn(i)
						imt2=l2+mnum3*(j-1)+nskpxn(i)
						yt(imt1)=yt(imt1)+wt10(j,ibnd)*xa(j)*tx(j)*xt(imt2)
						yt(imt1)=yt(imt1)+wt1m(j,ibnd)*xa(j)*tx(j)*xt(imt2-mnum3)
						yt(imt1)=yt(imt1)+wt1p(j,ibnd)*xa(j)*tx(j)*xt(imt2+mnum3)
					end do
					imt1=l1+mnum3*mjm2+nskpxn(i)
					imt2=l2+mnum3*mjm2+nskpxn(i)
					if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
						yt(imt1)=yt(imt1)+wt10(mjm1,2)*xa(mjm1)*tx(mjm1)*xt(imt2)
						yt(imt1)=yt(imt1)+wt1m(mjm1,2)*xa(mjm1)*tx(mjm1)*xt(imt2-mnum3)
!					else if (ivar == 4) then
!						yt(imt1)=yt(imt1)+wt10(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(imt2)
!						yt(imt1)=yt(imt1)+wt1m(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(imt2-mnum3)
!						yt(imt1)=yt(imt1)+dc1p(mjm1)*xa(mjm1)*tx(mjm1)*xt(imt2)*(r(mj)-r(mjm2))/ &
!								  (r(mjm1)-r(mjm2))
!						yt(imt1)=yt(imt1)-dc1p(mjm1)*xa(mjm1)*tx(mjm1)*xt(imt2-mnum3)*(r(mj)-r(mjm1))/ &
!								  (r(mjm1)-r(mjm2))
					else
						yt(imt1)=yt(imt1)+wt10(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(imt2)
						yt(imt1)=yt(imt1)+wt1m(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(imt2-mnum3)
					end if
				else
					imt1=l1+nskpxn(i)
					imt2=l2+nskpxn(i)
					yt(imt1)=yt(imt1)+wt20(1,ibnd)*xa(1)*tx(1)*xt(imt2)
					yt(imt1)=yt(imt1)+wt2p(1,ibnd)*xa(1)*tx(1)*xt(imt2+mnum3)
					do j=2,mjm2
						imt1=l1+mnum3*(j-1)+nskpxn(i)
						imt2=l2+mnum3*(j-1)+nskpxn(i)
						yt(imt1)=yt(imt1)+wt20(j,ibnd)*xa(j)*tx(j)*xt(imt2)
						yt(imt1)=yt(imt1)+wt2m(j,ibnd)*xa(j)*tx(j)*xt(imt2-mnum3)
						yt(imt1)=yt(imt1)+wt2p(j,ibnd)*xa(j)*tx(j)*xt(imt2+mnum3)
					end do
					imt1=l1+mnum3*mjm2+nskpxn(i)
					imt2=l2+mnum3*mjm2+nskpxn(i)
					if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
						yt(imt1)=yt(imt1)+wt20(mjm1,2)*xa(mjm1)*tx(mjm1)*xt(imt2)
						yt(imt1)=yt(imt1)+wt2m(mjm1,2)*xa(mjm1)*tx(mjm1)*xt(imt2-mnum3)
!					else if (ivar == 4) then
!						yt(imt1)=yt(imt1)+wt20(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(imt2)
!						yt(imt1)=yt(imt1)+wt2m(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(imt2-mnum3)
!						yt(imt1)=yt(imt1)+dc2p(mjm1)*xa(mjm1)*tx(mjm1)*xt(imt2)*(r(mj)-r(mjm2))/ &
!								  (r(mjm1)-r(mjm2))
!						yt(imt1)=yt(imt1)-dc2p(mjm1)*xa(mjm1)*tx(mjm1)*xt(imt2-mnum3)*(r(mj)-r(mjm1))/ &
!								  (r(mjm1)-r(mjm2))
					else
						yt(imt1)=yt(imt1)+wt20(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(imt2)
						yt(imt1)=yt(imt1)+wt2m(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(imt2-mnum3)
					end if
				end if
			end do
		end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(mnum2,mnum3,l2t,l,l2,l2p,m,n,ibnd,lp,l1,l1p,x,il,j,xa,imt1,imt2)
		do i=nc_start,nc_end
			mnum2=noeqn*mnumn(i)
			mnum3=2*noeqn*mnumn(i)
			do l2t=1,mnumn(i)
				l=l2t+lnumn(i-1)
				l2=l2t+(ivar-1)*mnumn(i)
				l2p=l2+mnum2
				m=mm(lln(l))
				n=nn(lln(l))
				ibnd=1
				if (isv == 1 .and. m == 0) ibnd=2
				l1=l2t+(ieqn-1)*mnumn(i)+ilp*mnum2
				l1p=l2t+(ieqn-1)*mnumn(i)+(1-ilp)*mnum2
				x=coef
				do il=1,ith
					x=x*m
				end do
				do il=1,izt
					x=x*n
				end do
				if (x == 0.) cycle
				x=x*ips
				do j=1,mjm1
					xa(j)=x*(rinv(j))**ith
				end do
				if (ir == 0) then
					do j=1,mjm1
						imt1=l1+mnum3*(j-1)+nskpxn(i)
						imt2=l2+mnum3*(j-1)+nskpxn(i)
						yt(imt1)=yt(imt1)+xa(j)*tx(j)*xt(imt2)
					end do
					do j=1,mjm1
						imt1=l1p+mnum3*(j-1)+nskpxn(i)
						imt2=l2p+mnum3*(j-1)+nskpxn(i)
						yt(imt1)=yt(imt1)+itz*xa(j)*tx(j)*xt(imt2)
					end do
				else if (ir == 1) then
					imt1=l1+nskpxn(i)
					imt2=l2+nskpxn(i)
					yt(imt1)=yt(imt1)+wt10(1,ibnd)*xa(1)*tx(1)*xt(imt2)
					yt(imt1)=yt(imt1)+wt1p(1,ibnd)*xa(1)*tx(1)*xt(imt2+mnum3)
					do j=2,mjm2
						imt1=l1+mnum3*(j-1)+nskpxn(i)
						imt2=l2+mnum3*(j-1)+nskpxn(i)
						yt(imt1)=yt(imt1)+wt10(j,ibnd)*xa(j)*tx(j)*xt(imt2)
						yt(imt1)=yt(imt1)+wt1m(j,ibnd)*xa(j)*tx(j)*xt(imt2-mnum3)
						yt(imt1)=yt(imt1)+wt1p(j,ibnd)*xa(j)*tx(j)*xt(imt2+mnum3)
					end do
					imt1=l1+mnum3*mjm2+nskpxn(i)
					imt2=l2+mnum3*mjm2+nskpxn(i)
					if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
						yt(imt1)=yt(imt1)+wt10(mjm1,2)*xa(mjm1)*tx(mjm1)*xt(imt2)
						yt(imt1)=yt(imt1)+wt1m(mjm1,2)*xa(mjm1)*tx(mjm1)*xt(imt2-mnum3)
!					else if (ivar == 4) then
!						yt(imt1)=yt(imt1)+wt10(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(imt2)
!						yt(imt1)=yt(imt1)+wt1m(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(imt2-mnum3)
!						yt(imt1)=yt(imt1)+dc1p(mjm1)*xa(mjm1)*tx(mjm1)*xt(imt2)*(r(mj)-r(mjm2))/ &
!								  (r(mjm1)-r(mjm2))
!						yt(imt1)=yt(imt1)-dc1p(mjm1)*xa(mjm1)*tx(mjm1)*xt(imt2-mnum3)*(r(mj)-r(mjm1))/ &
!								  (r(mjm1)-r(mjm2))
					else
						yt(imt1)=yt(imt1)+wt10(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(imt2)
						yt(imt1)=yt(imt1)+wt1m(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(imt2-mnum3)
					end if
					imt1=l1p+nskpxn(i)
					imt2=l2p+nskpxn(i)
					yt(imt1)=yt(imt1)+itz*wt10(1,ibnd)*xa(1)*tx(1)*xt(imt2)
					yt(imt1)=yt(imt1)+itz*wt1p(1,ibnd)*xa(1)*tx(1)*xt(imt2+mnum3)
					do j=2,mjm2
						imt1=l1p+mnum3*(j-1)+nskpxn(i)
						imt2=l2p+mnum3*(j-1)+nskpxn(i)
						yt(imt1)=yt(imt1)+itz*wt10(j,ibnd)*xa(j)*tx(j)*xt(imt2)
						yt(imt1)=yt(imt1)+itz*wt1m(j,ibnd)*xa(j)*tx(j)*xt(imt2-mnum3)
						yt(imt1)=yt(imt1)+itz*wt1p(j,ibnd)*xa(j)*tx(j)*xt(imt2+mnum3)
					end do
					imt1=l1p+mnum3*mjm2+nskpxn(i)
					imt2=l2p+mnum3*mjm2+nskpxn(i)
					if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
						yt(imt1)=yt(imt1)+itz*wt10(mjm1,2)*xa(mjm1)*tx(mjm1)*xt(imt2)
						yt(imt1)=yt(imt1)+itz*wt1m(mjm1,2)*xa(mjm1)*tx(mjm1)*xt(imt2-mnum3)
!					else if (ivar == 4) then
!						yt(imt1)=yt(imt1)+itz*wt10(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(imt2)
!						yt(imt1)=yt(imt1)+itz*wt1m(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(imt2-mnum3)
!						yt(imt1)=yt(imt1)+itz*dc1p(mjm1)*xa(mjm1)*tx(mjm1)*xt(imt2)*(r(mj)-r(mjm2))/ &
!								  (r(mjm1)-r(mjm2))
!						yt(imt1)=yt(imt1)-itz*dc1p(mjm1)*xa(mjm1)*tx(mjm1)*xt(imt2-mnum3)*(r(mj)-r(mjm1))/ &
!								  (r(mjm1)-r(mjm2))
					else
						yt(imt1)=yt(imt1)+itz*wt10(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(imt2)
						yt(imt1)=yt(imt1)+itz*wt1m(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(imt2-mnum3)
					end if
				else
					imt1=l1+nskpxn(i)
					imt2=l2+nskpxn(i)
					yt(imt1)=yt(imt1)+wt20(1,ibnd)*xa(1)*tx(1)*xt(imt2)
					yt(imt1)=yt(imt1)+wt2p(1,ibnd)*xa(1)*tx(1)*xt(imt2+mnum3)
					do j=2,mjm2
						imt1=l1+mnum3*(j-1)+nskpxn(i)
						imt2=l2+mnum3*(j-1)+nskpxn(i)
						yt(imt1)=yt(imt1)+wt20(j,ibnd)*xa(j)*tx(j)*xt(imt2)
						yt(imt1)=yt(imt1)+wt2m(j,ibnd)*xa(j)*tx(j)*xt(imt2-mnum3)
						yt(imt1)=yt(imt1)+wt2p(j,ibnd)*xa(j)*tx(j)*xt(imt2+mnum3)
					end do
					imt1=l1+mnum3*mjm2+nskpxn(i)
					imt2=l2+mnum3*mjm2+nskpxn(i)
					if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
						yt(imt1)=yt(imt1)+wt20(mjm1,2)*xa(mjm1)*tx(mjm1)*xt(imt2)
						yt(imt1)=yt(imt1)+wt2m(mjm1,2)*xa(mjm1)*tx(mjm1)*xt(imt2-mnum3)
!					else if (ivar == 4) then
!						yt(imt1)=yt(imt1)+wt20(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(imt2)
!						yt(imt1)=yt(imt1)+wt2m(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(imt2-mnum3)
!						yt(imt1)=yt(imt1)+dc2p(mjm1)*xa(mjm1)*tx(mjm1)*xt(imt2)*(r(mj)-r(mjm2))/ &
!								  (r(mjm1)-r(mjm2))
!						yt(imt1)=yt(imt1)-dc2p(mjm1)*xa(mjm1)*tx(mjm1)*xt(imt2-mnum3)*(r(mj)-r(mjm1))/ &
!								  (r(mjm1)-r(mjm2))
					else
						yt(imt1)=yt(imt1)+wt20(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(imt2)
						yt(imt1)=yt(imt1)+wt2m(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(imt2-mnum3)
					end if
					imt1=l1p+nskpxn(i)
					imt2=l2p+nskpxn(i)
					yt(imt1)=yt(imt1)+itz*wt20(1,ibnd)*xa(1)*tx(1)*xt(imt2)
					yt(imt1)=yt(imt1)+itz*wt2p(1,ibnd)*xa(1)*tx(1)*xt(imt2+mnum3)
					do j=2,mjm2
						imt1=l1p+mnum3*(j-1)+nskpxn(i)
						imt2=l2p+mnum3*(j-1)+nskpxn(i)
						yt(imt1)=yt(imt1)+itz*wt20(j,ibnd)*xa(j)*tx(j)*xt(imt2)
						yt(imt1)=yt(imt1)+itz*wt2m(j,ibnd)*xa(j)*tx(j)*xt(imt2-mnum3)
						yt(imt1)=yt(imt1)+itz*wt2p(j,ibnd)*xa(j)*tx(j)*xt(imt2+mnum3)
					end do
					imt1=l1p+mnum3*mjm2+nskpxn(i)
					imt2=l2p+mnum3*mjm2+nskpxn(i)
					if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
						yt(imt1)=yt(imt1)+itz*wt20(mjm1,2)*xa(mjm1)*tx(mjm1)*xt(imt2)
						yt(imt1)=yt(imt1)+itz*wt2m(mjm1,2)*xa(mjm1)*tx(mjm1)*xt(imt2-mnum3)
!					else if (ivar == 4) then
!						yt(imt1)=yt(imt1)+itz*wt20(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(imt2)
!						yt(imt1)=yt(imt1)+itz*wt2m(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(imt2-mnum3)
!						yt(imt1)=yt(imt1)+itz*dc2p(mjm1)*xa(mjm1)*tx(mjm1)*xt(imt2)*(r(mj)-r(mjm2))/ &
!								  (r(mjm1)-r(mjm2))
!						yt(imt1)=yt(imt1)-itz*dc2p(mjm1)*xa(mjm1)*tx(mjm1)*xt(imt2-mnum3)*(r(mj)-r(mjm1))/ &
!								  (r(mjm1)-r(mjm2))
					else
						yt(imt1)=yt(imt1)+itz*wt20(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(imt2)
						yt(imt1)=yt(imt1)+itz*wt2m(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(imt2-mnum3)
					end if
				end if
			end do
		end do
!$OMP END PARALLEL DO

	end subroutine b2lx0

	subroutine b2lx0_dlsq(tx,ieqn,ivar,coef)

		use param
		use var_para
		use cotrol
		use domain
		use dynamo
		implicit none

		integer :: ieqn,ivar
		real(IDP) :: coef
		real(IDP), dimension(0:) :: tx
		integer :: i,l,m,l1,l2,l2t,j,mnum2,mnum3,imt1,imt2

		do i=nr_start,nr_end
			mnum3=noeqn*mnumn(i)
			do l2t=1,mnumn(i)
				l=l2t+lnumn(i-1)
				l2=l2t+(ivar-1)*mnumn(i)
				m=mm(lln(l))
				l1=l2t+(ieqn-1)*mnumn(i)
				do j=1,mjm1
					imt1=l1+mnum3*(j-1)+nskpxn(i)
					imt2=l2+mnum3*(j-1)+nskpxn(i)
					yt(imt1)=yt(imt1)-coef*m*m*rinv(j)*rinv(j)*tx(j)*xt(imt2)
				end do
				imt1=l1+nskpxn(i)
				imt2=l2+nskpxn(i)
				if (m == 0) then
					yt(imt1)=yt(imt1)-coef*tx(1)*(del2cm(1)+del2cp(1)-r(2)**2*del2cm(1)/(r(2)**2-r(1)**2))*xt(imt2)
					yt(imt1)=yt(imt1)+coef*tx(1)*(del2cp(1)-r(1)**2*del2cm(1)/(r(2)**2-r(1)**2))*xt(imt2+mnum3)
				else
					yt(imt1)=yt(imt1)-coef*tx(1)*(del2cm(1)+del2cp(1))*xt(imt2)
					yt(imt1)=yt(imt1)+coef*tx(1)*del2cp(1)*xt(imt2+mnum3)
				end if
				do j=2,mjm2
					imt1=l1+mnum3*(j-1)+nskpxn(i)
					imt2=l2+mnum3*(j-1)+nskpxn(i)
					yt(imt1)=yt(imt1)-coef*tx(j)*(del2cm(j)+del2cp(j))*xt(imt2)
					yt(imt1)=yt(imt1)+coef*tx(j)*del2cm(j)*xt(imt2-mnum3)
					yt(imt1)=yt(imt1)+coef*tx(j)*del2cp(j)*xt(imt2+mnum3)
				end do
				imt1=l1+mnum3*mjm2+nskpxn(i)
				imt2=l2+mnum3*mjm2+nskpxn(i)
				if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
					yt(imt1)=yt(imt1)-coef*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1)-del2cp(mjm1)* &
							  (r(mj)-r(mjm2))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2))*xt(imt2)
					yt(imt1)=yt(imt1)+coef*tx(mjm1)*(del2cm(mjm1)-del2cp(mjm1)* &
							  (r(mj)-r(mjm1))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2))*xt(imt2-mnum3)
!				else if (ivar == 4) then
!					yt(imt1)=yt(imt1)-coef*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1)-del2cp(mjm1)* &
!							  (r(mj)-r(mjm2))/(r(mjm1)-r(mjm2)))*xt(imt2)
!					yt(imt1)=yt(imt1)+coef*tx(mjm1)*(del2cm(mjm1)-del2cp(mjm1)* &
!							  (r(mj)-r(mjm1))/(r(mjm1)-r(mjm2)))*xt(imt2-mnum3)
				else
					yt(imt1)=yt(imt1)-coef*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1))*xt(imt2)
					yt(imt1)=yt(imt1)+coef*tx(mjm1)*del2cm(mjm1)*xt(imt2-mnum3)
				end if
			end do
		end do
		do i=nc_start,nc_end
			mnum2=noeqn*mnumn(i)
			mnum3=2*noeqn*mnumn(i)
			do l2t=1,mnumn(i)
				l=l2t+lnumn(i-1)
				l2=l2t+(ivar-1)*mnumn(i)
				m=mm(lln(l))
				l1=l2t+(ieqn-1)*mnumn(i)
				do j=1,mjm1
					imt1=l1+mnum3*(j-1)+nskpxn(i)
					imt2=l2+mnum3*(j-1)+nskpxn(i)
					yt(imt1)=yt(imt1)-coef*m*m*rinv(j)*rinv(j)*tx(j)*xt(imt2)
					yt(imt1+mnum2)=yt(imt1+mnum2)-coef*m*m*rinv(j)*rinv(j)*tx(j)*xt(imt2+mnum2)
				end do
				imt1=l1+nskpxn(i)
				imt2=l2+nskpxn(i)
				if (m == 0) then
					yt(imt1)=yt(imt1)-coef*tx(1)*(del2cm(1)+del2cp(1)-r(2)**2*del2cm(1)/(r(2)**2-r(1)**2))*xt(imt2)
					yt(imt1)=yt(imt1)+coef*tx(1)*(del2cp(1)-r(1)**2*del2cm(1)/(r(2)**2-r(1)**2))*xt(imt2+mnum3)
				else
					yt(imt1)=yt(imt1)-coef*tx(1)*(del2cm(1)+del2cp(1))*xt(imt2)
					yt(imt1)=yt(imt1)+coef*tx(1)*del2cp(1)*xt(imt2+mnum3)
				end if
				do j=2,mjm2
					imt1=l1+mnum3*(j-1)+nskpxn(i)
					imt2=l2+mnum3*(j-1)+nskpxn(i)
					yt(imt1)=yt(imt1)-coef*tx(j)*(del2cm(j)+del2cp(j))*xt(imt2)
					yt(imt1)=yt(imt1)+coef*tx(j)*del2cm(j)*xt(imt2-mnum3)
					yt(imt1)=yt(imt1)+coef*tx(j)*del2cp(j)*xt(imt2+mnum3)
				end do
				imt1=l1+mnum3*mjm2+nskpxn(i)
				imt2=l2+mnum3*mjm2+nskpxn(i)
				if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
					yt(imt1)=yt(imt1)-coef*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1)-del2cp(mjm1)* &
							  (r(mj)-r(mjm2))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2))*xt(imt2)
					yt(imt1)=yt(imt1)+coef*tx(mjm1)*(del2cm(mjm1)-del2cp(mjm1)* &
							  (r(mj)-r(mjm1))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2))*xt(imt2-mnum3)
!				else if (ivar == 4) then
!					yt(imt1)=yt(imt1)-coef*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1)-del2cp(mjm1)* &
!							  (r(mj)-r(mjm2))/(r(mjm1)-r(mjm2)))*xt(imt2)
!					yt(imt1)=yt(imt1)+coef*tx(mjm1)*(del2cm(mjm1)-del2cp(mjm1)* &
!							  (r(mj)-r(mjm1))/(r(mjm1)-r(mjm2)))*xt(imt2-mnum3)
				else
					yt(imt1)=yt(imt1)-coef*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1))*xt(imt2)
					yt(imt1)=yt(imt1)+coef*tx(mjm1)*del2cm(mjm1)*xt(imt2-mnum3)
				end if
				imt1=l1+mnum2+nskpxn(i)
				imt2=l2+mnum2+nskpxn(i)
				if (m == 0) then
					yt(imt1)=yt(imt1)-coef*tx(1)*(del2cm(1)+del2cp(1)-r(2)**2*del2cm(1)/(r(2)**2-r(1)**2))*xt(imt2)
					yt(imt1)=yt(imt1)+coef*tx(1)*(del2cp(1)-r(1)**2*del2cm(1)/(r(2)**2-r(1)**2))*xt(imt2+mnum3)
				else
					yt(imt1)=yt(imt1)-coef*tx(1)*(del2cm(1)+del2cp(1))*xt(imt2)
					yt(imt1)=yt(imt1)+coef*tx(1)*del2cp(1)*xt(imt2+mnum3)
				end if
				do j=2,mjm2
					imt1=l1+mnum2+mnum3*(j-1)+nskpxn(i)
					imt2=l2+mnum2+mnum3*(j-1)+nskpxn(i)
					yt(imt1)=yt(imt1)-coef*tx(j)*(del2cm(j)+del2cp(j))*xt(imt2)
					yt(imt1)=yt(imt1)+coef*tx(j)*del2cm(j)*xt(imt2-mnum3)
					yt(imt1)=yt(imt1)+coef*tx(j)*del2cp(j)*xt(imt2+mnum3)
				end do
				imt1=l1+mnum2+mnum3*mjm2+nskpxn(i)
				imt2=l2+mnum2+mnum3*mjm2+nskpxn(i)
				if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
					yt(imt1)=yt(imt1)-coef*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1)-del2cp(mjm1)* &
							  (r(mj)-r(mjm2))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2))*xt(imt2)
					yt(imt1)=yt(imt1)+coef*tx(mjm1)*(del2cm(mjm1)-del2cp(mjm1)* &
							  (r(mj)-r(mjm1))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2))*xt(imt2-mnum3)
!				else if (ivar == 4) then
!					yt(imt1)=yt(imt1)-coef*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1)-del2cp(mjm1)* &
!							  (r(mj)-r(mjm2))/(r(mjm1)-r(mjm2)))*xt(imt2)
!					yt(imt1)=yt(imt1)+coef*tx(mjm1)*(del2cm(mjm1)-del2cp(mjm1)* &
!							  (r(mj)-r(mjm1))/(r(mjm1)-r(mjm2)))*xt(imt2-mnum3)
				else
					yt(imt1)=yt(imt1)-coef*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1))*xt(imt2)
					yt(imt1)=yt(imt1)+coef*tx(mjm1)*del2cm(mjm1)*xt(imt2-mnum3)
				end if
			end do
		end do

	end subroutine b2lx0_dlsq

	subroutine b2lx_dlsq(ieqn,ivar,coef)

		use param
		use equil
		implicit none

		integer :: ieqn,ivar
		real(IDP) :: coef

		interface
			subroutine b2lx(tx,itx,ieqn,ivar,ith,ir,izt,coef)
				use param
				implicit none
				integer :: itx,ieqn,ivar,ith,ir,izt
				real(IDP) :: coef
				real(IDP), dimension(0:,0:) :: tx
			end subroutine b2lx
		end interface

		call b2lx(lplrr,1,ieqn,ivar,0,2,0,coef)
		call b2lx(lplrt,-1,ieqn,ivar,1,1,0,coef)
		call b2lx(lplrz,-1,ieqn,ivar,0,1,1,coef)
		call b2lx(lpltt,1,ieqn,ivar,2,0,0,coef)
		call b2lx(lpltz,1,ieqn,ivar,1,0,1,coef)
		call b2lx(lplzz,1,ieqn,ivar,0,0,2,coef)
		call b2lx(lplr,1,ieqn,ivar,0,1,0,coef)
		call b2lx(lplt,-1,ieqn,ivar,1,0,0,coef)
		call b2lx(lplz,-1,ieqn,ivar,0,0,1,coef)

	end subroutine b2lx_dlsq

	subroutine b2lx_landau_grad_parallel(tx,itx,ieqn,ivar,ith,ir,izt,coef)

		use param
		use var_para
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		integer :: itx,ieqn,ivar,ith,ir,izt
		real(IDP) :: coef,x
		real(IDP), dimension(0:,0:) :: tx
		real(IDP), dimension(jdim) :: xa
		integer :: ity,isv,ips,iph,itz,icn,icg,ilp,ich,i,mnum2,mnum3,l,m,n,l1,l2,l3,l1t,l2t,l3n,l1p,l2p,lp,lpp,lpn,ibnd,j,imt1,imt2,il

		ity= 1
		if (ivar == 2 .or. ivar == 4 .or. ivar == 6 .or. ivar == 7 .or. ivar == ivalp .or. ivar > iq) ity=-1
		isv=1
		ips=(-1)**((ith+izt+1)/2)
		iph=ity**(ith+izt)*ips
		itz=(-1)**(ith+izt)
		ity=ity*itz
		icn=itx*itz
		ilp=(1-icn)/2
		ich=itx*ity
		if (ieqn == 2 .or. ieqn == 4 .or. ieqn == 6 .or. ieqn == 7 .or. ieqn == ivalp .or. ieqn > iq) ich=-ich
		do i=nr_start,nr_end
			mnum3=noeqn*mnumn(i)
			do l2t=1,mnumn(i)
				l=l2t+lnumn(i-1)
				l2=l2t+(ivar-1)*mnumn(i)
				m=mm(lln(l))
				n=nn(lln(l))
				ibnd=1
				if (isv == 1 .and. m == 0) ibnd=2
				do l1t=1,mnumn(i)
					lp=l1t+lnumn(i-1)
					lpp=lp
					if (ich == -1) lpp=lo(lp)
					if (lpp == 0) cycle
					l1=lpp+(ieqn-1)*mnumn(i)-lnumn(i-1)
					lpn=nskp2n(i)+mnumn(i)*(l2t-1)+l1t
					do l3=1,leqmax
						l3n=l3
						if (ity ==  1 .and. itx ==  1) x=coef*cmapp(lpn,l3n)
						if (ity ==  1 .and. itx == -1) x=coef*cmamp(lpn,l3n)
						if (ity == -1 .and. itx == -1) x=coef*cmamm(lpn,l3n)
						if (ity == -1 .and. itx ==  1) x=coef*cmapm(lpn,l3n)
						do il=1,ith
							x=x*m
						end do
						do il=1,izt
							x=x*n
						end do
						if (x == 0.) cycle
						x=x*iph
						do j=1,mjm1
							xa(j)=x*(rinv(j))**ith*abs(n-m*qqinv(j))
						end do
						if (ir == 0) then
							do j=1,mjm1
								imt1=l1+mnum3*(j-1)+nskpxn(i)
								imt2=l2+mnum3*(j-1)+nskpxn(i)
								yt(imt1)=yt(imt1)+xa(j)*tx(j,l3n)*xt(imt2)
							end do
						else if (ir == 1) then
							imt1=l1+nskpxn(i)
							imt2=l2+nskpxn(i)
							yt(imt1)=yt(imt1)+wt10(1,ibnd)*xa(1)*tx(1,l3n)*xt(imt2)
							yt(imt1)=yt(imt1)+wt1p(1,ibnd)*xa(1)*tx(1,l3n)*xt(imt2+mnum3)
							do j=2,mjm1
								imt1=l1+mnum3*(j-1)+nskpxn(i)
								imt2=l2+mnum3*(j-1)+nskpxn(i)
								yt(imt1)=yt(imt1)+wt10(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2)
								yt(imt1)=yt(imt1)+wt1m(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2-mnum3)
								yt(imt1)=yt(imt1)+wt1p(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2+mnum3)
							end do
						else
							imt1=l1+nskpxn(i)
							imt2=l2+nskpxn(i)
							yt(imt1)=yt(imt1)+wt20(1,ibnd)*xa(1)*tx(1,l3n)*xt(imt2)
							yt(imt1)=yt(imt1)+wt2p(1,ibnd)*xa(1)*tx(1,l3n)*xt(imt2+mnum3)
							do j=2,mjm1
								imt1=l1+mnum3*(j-1)+nskpxn(i)
								imt2=l2+mnum3*(j-1)+nskpxn(i)
								yt(imt1)=yt(imt1)+wt20(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2)
								yt(imt1)=yt(imt1)+wt2m(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2-mnum3)
								yt(imt1)=yt(imt1)+wt2p(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2+mnum3)
							end do
						end if
					end do
				end do
			end do
		end do
		do i=nc_start,nc_end
			mnum2=noeqn*mnumn(i)
			mnum3=2*noeqn*mnumn(i)
			do l2t=1,mnumn(i)
				l=l2t+lnumn(i-1)
				l2=l2t+(ivar-1)*mnumn(i)
				l2p=l2+mnum2
				m=mm(lln(l))
				n=nn(lln(l))
				ibnd=1
				if (isv == 1 .and. m == 0) ibnd=2
				do l1t=1,mnumn(i)
					l1=l1t+(ieqn-1)*mnumn(i)+ilp*mnum2
					l1p=l1t+(ieqn-1)*mnumn(i)+(1-ilp)*mnum2
					lpn=nskp2n(i)+mnumn(i)*(l2t-1)+l1t
					do l3=1,leqmax
						l3n=l3
						icg=icn
						if (itz ==  1 .and. itx ==  1) x=coef*cmapp(lpn,l3n)
						if (itz ==  1 .and. itx == -1) x=coef*cmamp(lpn,l3n)
						if (itz == -1 .and. itx == -1) x=coef*cmamm(lpn,l3n)
						if (itz == -1 .and. itx ==  1) x=coef*cmapm(lpn,l3n)
						if (l3 == leq0 .and. itx == -1 .and. l1t == l2t) then
							icg=-icn
							l1=l1t+(ieqn-1)*mnumn(i)+(1-ilp)*mnum2
							l1p=l1t+(ieqn-1)*mnumn(i)+ilp*mnum2
							x=coef
						end if
						do il=1,ith
							x=x*m
						end do
						do il=1,izt
							x=x*n
						end do
						if (x == 0.) cycle
						x=x*ips
						do j=1,mjm1
							xa(j)=x*(rinv(j))**ith*abs(n-m*qqinv(j))
						end do
						if (ir == 0) then
							do j=1,mjm1
								imt1=l1+mnum3*(j-1)+nskpxn(i)
								imt2=l2+mnum3*(j-1)+nskpxn(i)
								yt(imt1)=yt(imt1)+xa(j)*tx(j,l3n)*xt(imt2)
							end do
							do j=1,mjm1
								imt1=l1p+mnum3*(j-1)+nskpxn(i)
								imt2=l2p+mnum3*(j-1)+nskpxn(i)
								yt(imt1)=yt(imt1)+icg*xa(j)*tx(j,l3n)*xt(imt2)
							end do
						else if (ir == 1) then
							imt1=l1+nskpxn(i)
							imt2=l2+nskpxn(i)
							yt(imt1)=yt(imt1)+wt10(1,ibnd)*xa(1)*tx(1,l3n)*xt(imt2)
							yt(imt1)=yt(imt1)+wt1p(1,ibnd)*xa(1)*tx(1,l3n)*xt(imt2+mnum3)
							do j=2,mjm1
								imt1=l1+mnum3*(j-1)+nskpxn(i)
								imt2=l2+mnum3*(j-1)+nskpxn(i)
								yt(imt1)=yt(imt1)+wt10(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2)
								yt(imt1)=yt(imt1)+wt1m(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2-mnum3)
								yt(imt1)=yt(imt1)+wt1p(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2+mnum3)
							end do
							imt1=l1p+nskpxn(i)
							imt2=l2p+nskpxn(i)
							yt(imt1)=yt(imt1)+icg*wt10(1,ibnd)*xa(1)*tx(1,l3n)*xt(imt2)
							yt(imt1)=yt(imt1)+icg*wt1p(1,ibnd)*xa(1)*tx(1,l3n)*xt(imt2+mnum3)
							do j=2,mjm1
								imt1=l1p+mnum3*(j-1)+nskpxn(i)
								imt2=l2p+mnum3*(j-1)+nskpxn(i)
								yt(imt1)=yt(imt1)+icg*wt10(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2)
								yt(imt1)=yt(imt1)+icg*wt1m(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2-mnum3)
								yt(imt1)=yt(imt1)+icg*wt1p(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2+mnum3)
							end do
						else
							imt1=l1+nskpxn(i)
							imt2=l2+nskpxn(i)
							yt(imt1)=yt(imt1)+wt20(1,ibnd)*xa(1)*tx(1,l3n)*xt(imt2)
							yt(imt1)=yt(imt1)+wt2p(1,ibnd)*xa(1)*tx(1,l3n)*xt(imt2+mnum3)
							do j=2,mjm1
								imt1=l1+mnum3*(j-1)+nskpxn(i)
								imt2=l2+mnum3*(j-1)+nskpxn(i)
								yt(imt1)=yt(imt1)+wt20(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2)
								yt(imt1)=yt(imt1)+wt2m(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2-mnum3)
								yt(imt1)=yt(imt1)+wt2p(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2+mnum3)
							end do
							imt1=l1p+nskpxn(i)
							imt2=l2p+nskpxn(i)
							yt(imt1)=yt(imt1)+icg*wt20(1,ibnd)*xa(1)*tx(1,l3n)*xt(imt2)
							yt(imt1)=yt(imt1)+icg*wt2p(1,ibnd)*xa(1)*tx(1,l3n)*xt(imt2+mnum3)
							do j=2,mjm1
								imt1=l1p+mnum3*(j-1)+nskpxn(i)
								imt2=l2p+mnum3*(j-1)+nskpxn(i)
								yt(imt1)=yt(imt1)+icg*wt20(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2)
								yt(imt1)=yt(imt1)+icg*wt2m(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2-mnum3)
								yt(imt1)=yt(imt1)+icg*wt2p(j,ibnd)*xa(j)*tx(j,l3n)*xt(imt2+mnum3)
							end do
						end if
					end do
				end do
			end do
		end do

	end subroutine b2lx_landau_grad_parallel	

	subroutine b2lx0_dlsqr(tx,ieqn,ivar,coef,prf)

		use param
		use var_para
		use cotrol
		use domain
		use dynamo
		implicit none

		integer :: ieqn,ivar
		real(IDP), dimension(:) :: coef
		real(IDP), dimension(0:) :: prf
		real(IDP), dimension(0:) :: tx
		integer :: i,l,m,l1,l2,l2t,j,mnum2,mnum3,imt1,imt2

		do i=nr_start,nr_end
			mnum3=noeqn*mnumn(i)
			do l2t=1,mnumn(i)
				l=l2t+lnumn(i-1)
				l2=l2t+(ivar-1)*mnumn(i)
				m=mm(lln(l))
				l1=l2t+(ieqn-1)*mnumn(i)
				do j=1,mjm1
					imt1=l1+mnum3*(j-1)+nskpxn(i)
					imt2=l2+mnum3*(j-1)+nskpxn(i)
					yt(imt1)=yt(imt1)-coef(i)*m*m*rinv(j)*rinv(j)*tx(j)*xt(imt2)*prf(j)
				end do
				imt1=l1+nskpxn(i)
				imt2=l2+nskpxn(i)
				if (m == 0) then
					yt(imt1)=yt(imt1)-coef(i)*tx(1)*(del2cm(1)+del2cp(1)-r(2)**2*del2cm(1)/(r(2)**2-r(1)**2))*xt(imt2)*prf(j)
					yt(imt1)=yt(imt1)+coef(i)*tx(1)*(del2cp(1)-r(1)**2*del2cm(1)/(r(2)**2-r(1)**2))*xt(imt2+mnum3)*prf(j)
				else
					yt(imt1)=yt(imt1)-coef(i)*tx(1)*(del2cm(1)+del2cp(1))*xt(imt2)*prf(j)
					yt(imt1)=yt(imt1)+coef(i)*tx(1)*del2cp(1)*xt(imt2+mnum3)*prf(j)
				end if
				do j=2,mjm2
					imt1=l1+mnum3*(j-1)+nskpxn(i)
					imt2=l2+mnum3*(j-1)+nskpxn(i)
					yt(imt1)=yt(imt1)-coef(i)*tx(j)*(del2cm(j)+del2cp(j))*xt(imt2)*prf(j)
					yt(imt1)=yt(imt1)+coef(i)*tx(j)*del2cm(j)*xt(imt2-mnum3)*prf(j)
					yt(imt1)=yt(imt1)+coef(i)*tx(j)*del2cp(j)*xt(imt2+mnum3)*prf(j)
				end do
				imt1=l1+mnum3*mjm2+nskpxn(i)
				imt2=l2+mnum3*mjm2+nskpxn(i)
				if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
					yt(imt1)=yt(imt1)-coef(i)*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1)-del2cp(mjm1)* &
							  (r(mj)-r(mjm2))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2))*xt(imt2)*prf(j)
					yt(imt1)=yt(imt1)+coef(i)*tx(mjm1)*(del2cm(mjm1)-del2cp(mjm1)* &
							  (r(mj)-r(mjm1))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2))*xt(imt2-mnum3)*prf(j)
!				else if (ivar == 4) then
!					yt(imt1)=yt(imt1)-coef(i)*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1)-del2cp(mjm1)* &
!							  (r(mj)-r(mjm2))/(r(mjm1)-r(mjm2)))*xt(imt2)*prf(j)
!					yt(imt1)=yt(imt1)+coef(i)*tx(mjm1)*(del2cm(mjm1)-del2cp(mjm1)* &
!							  (r(mj)-r(mjm1))/(r(mjm1)-r(mjm2)))*xt(imt2-mnum3)*prf(j)
				else
					yt(imt1)=yt(imt1)-coef(i)*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1))*xt(imt2)*prf(j)
					yt(imt1)=yt(imt1)+coef(i)*tx(mjm1)*del2cm(mjm1)*xt(imt2-mnum3)*prf(j)
				end if
			end do
		end do
		do i=nc_start,nc_end
			mnum2=noeqn*mnumn(i)
			mnum3=2*noeqn*mnumn(i)
			do l2t=1,mnumn(i)
				l=l2t+lnumn(i-1)
				l2=l2t+(ivar-1)*mnumn(i)
				m=mm(lln(l))
				l1=l2t+(ieqn-1)*mnumn(i)
				do j=1,mjm1
					imt1=l1+mnum3*(j-1)+nskpxn(i)
					imt2=l2+mnum3*(j-1)+nskpxn(i)
					yt(imt1)=yt(imt1)-coef(i)*m*m*rinv(j)*rinv(j)*tx(j)*xt(imt2)*prf(j)
					yt(imt1+mnum2)=yt(imt1+mnum2)-coef(i)*m*m*rinv(j)*rinv(j)*tx(j)*xt(imt2+mnum2)*prf(j)
				end do
				imt1=l1+nskpxn(i)
				imt2=l2+nskpxn(i)
				if (m == 0) then
					yt(imt1)=yt(imt1)-coef(i)*tx(1)*(del2cm(1)+del2cp(1)-r(2)**2*del2cm(1)/(r(2)**2-r(1)**2))*xt(imt2)*prf(j)
					yt(imt1)=yt(imt1)+coef(i)*tx(1)*(del2cp(1)-r(1)**2*del2cm(1)/(r(2)**2-r(1)**2))*xt(imt2+mnum3)*prf(j)
				else
					yt(imt1)=yt(imt1)-coef(i)*tx(1)*(del2cm(1)+del2cp(1))*xt(imt2)*prf(j)
					yt(imt1)=yt(imt1)+coef(i)*tx(1)*del2cp(1)*xt(imt2+mnum3)*prf(j)
				end if
				do j=2,mjm2
					imt1=l1+mnum3*(j-1)+nskpxn(i)
					imt2=l2+mnum3*(j-1)+nskpxn(i)
					yt(imt1)=yt(imt1)-coef(i)*tx(j)*(del2cm(j)+del2cp(j))*xt(imt2)*prf(j)
					yt(imt1)=yt(imt1)+coef(i)*tx(j)*del2cm(j)*xt(imt2-mnum3)*prf(j)
					yt(imt1)=yt(imt1)+coef(i)*tx(j)*del2cp(j)*xt(imt2+mnum3)*prf(j)
				end do
				imt1=l1+mnum3*mjm2+nskpxn(i)
				imt2=l2+mnum3*mjm2+nskpxn(i)
				if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
					yt(imt1)=yt(imt1)-coef(i)*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1)-del2cp(mjm1)* &
							  (r(mj)-r(mjm2))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2))*xt(imt2)*prf(j)
					yt(imt1)=yt(imt1)+coef(i)*tx(mjm1)*(del2cm(mjm1)-del2cp(mjm1)* &
							  (r(mj)-r(mjm1))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2))*xt(imt2-mnum3)*prf(j)
!				else if (ivar == 4) then
!					yt(imt1)=yt(imt1)-coef(i)*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1)-del2cp(mjm1)* &
!							  (r(mj)-r(mjm2))/(r(mjm1)-r(mjm2)))*xt(imt2)*prf(j)
!					yt(imt1)=yt(imt1)+coef(i)*tx(mjm1)*(del2cm(mjm1)-del2cp(mjm1)* &
!							  (r(mj)-r(mjm1))/(r(mjm1)-r(mjm2)))*xt(imt2-mnum3)*prf(j)
				else
					yt(imt1)=yt(imt1)-coef(i)*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1))*xt(imt2)*prf(j)
					yt(imt1)=yt(imt1)+coef(i)*tx(mjm1)*del2cm(mjm1)*xt(imt2-mnum3)*prf(j)
				end if
				imt1=l1+mnum2+nskpxn(i)
				imt2=l2+mnum2+nskpxn(i)
				if (m == 0) then
					yt(imt1)=yt(imt1)-coef(i)*tx(1)*(del2cm(1)+del2cp(1)-r(2)**2*del2cm(1)/(r(2)**2-r(1)**2))*xt(imt2)*prf(j)
					yt(imt1)=yt(imt1)+coef(i)*tx(1)*(del2cp(1)-r(1)**2*del2cm(1)/(r(2)**2-r(1)**2))*xt(imt2+mnum3)*prf(j)
				else
					yt(imt1)=yt(imt1)-coef(i)*tx(1)*(del2cm(1)+del2cp(1))*xt(imt2)*prf(j)
					yt(imt1)=yt(imt1)+coef(i)*tx(1)*del2cp(1)*xt(imt2+mnum3)*prf(j)
				end if
				do j=2,mjm2
					imt1=l1+mnum2+mnum3*(j-1)+nskpxn(i)
					imt2=l2+mnum2+mnum3*(j-1)+nskpxn(i)
					yt(imt1)=yt(imt1)-coef(i)*tx(j)*(del2cm(j)+del2cp(j))*xt(imt2)*prf(j)
					yt(imt1)=yt(imt1)+coef(i)*tx(j)*del2cm(j)*xt(imt2-mnum3)*prf(j)
					yt(imt1)=yt(imt1)+coef(i)*tx(j)*del2cp(j)*xt(imt2+mnum3)*prf(j)
				end do
				imt1=l1+mnum2+mnum3*mjm2+nskpxn(i)
				imt2=l2+mnum2+mnum3*mjm2+nskpxn(i)
				if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
					yt(imt1)=yt(imt1)-coef(i)*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1)-del2cp(mjm1)* &
							  (r(mj)-r(mjm2))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2))*xt(imt2)*prf(j)
					yt(imt1)=yt(imt1)+coef(i)*tx(mjm1)*(del2cm(mjm1)-del2cp(mjm1)* &
							  (r(mj)-r(mjm1))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2))*xt(imt2-mnum3)*prf(j)
!				else if (ivar == 4) then
!					yt(imt1)=yt(imt1)-coef(i)*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1)-del2cp(mjm1)* &
!							  (r(mj)-r(mjm2))/(r(mjm1)-r(mjm2)))*xt(imt2)*prf(j)
!					yt(imt1)=yt(imt1)+coef(i)*tx(mjm1)*(del2cm(mjm1)-del2cp(mjm1)* &
!							  (r(mj)-r(mjm1))/(r(mjm1)-r(mjm2)))*xt(imt2-mnum3)*prf(j)
				else
					yt(imt1)=yt(imt1)-coef(i)*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1))*xt(imt2)*prf(j)
					yt(imt1)=yt(imt1)+coef(i)*tx(mjm1)*del2cm(mjm1)*xt(imt2-mnum3)*prf(j)
				end if
			end do
		end do

	end subroutine b2lx0_dlsqr