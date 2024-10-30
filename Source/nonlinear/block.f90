	subroutine blockj(tx,itx,ieqn,ivar,ith,ir,izt,coef)

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
		integer :: ity,isv,ips,iph,itz,icn,ilp,ich,i,mnum3,l,m,n,l1,l2,l3,l1t,l2t,l3n,lp,lpp,lpn,ibnd,j,imat,il,lpq,lq,lq2,imtq

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
								imat=l1+mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
								amat(imat)=amat(imat)+xa(j)*tx(j,l3n)
							end do
						else if (ir == 1) then
							do j=1,mjm2
								imat=l1+mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
								amat(imat)=amat(imat)+wt10(j,ibnd)*xa(j)*tx(j,l3n)
								cmat(imat)=cmat(imat)+wt1m(j,ibnd)*xa(j)*tx(j,l3n)
								bmat(imat)=bmat(imat)+wt1p(j,ibnd)*xa(j)*tx(j,l3n)
							end do
							imat=l1+mnum3*(l2-1+mnum3*mjm2)+nskpn(i)
							if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
								amat(imat)=amat(imat)+wt10(mjm1,2)*xa(mjm1)*tx(mjm1,l3n)
								cmat(imat)=cmat(imat)+wt1m(mjm1,2)*xa(mjm1)*tx(mjm1,l3n)
								bmat(imat)=bmat(imat)+wt1p(mjm1,2)*xa(mjm1)*tx(mjm1,l3n)
!							else if (ivar == 4) then
!								amat(imat)=amat(imat)+wt10(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)
!								cmat(imat)=cmat(imat)+wt1m(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)
!								bmat(imat)=bmat(imat)+wt1p(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)
!								amat(imat)=amat(imat)+dc1p(mjm1)*xa(mjm1)*tx(mjm1,l3n)*(r(mj)-r(mjm2))/ &
!										      (r(mjm1)-r(mjm2))
!								cmat(imat)=cmat(imat)-dc1p(mjm1)*xa(mjm1)*tx(mjm1,l3n)*(r(mj)-r(mjm1))/ &
!										      (r(mjm1)-r(mjm2))
							else
								amat(imat)=amat(imat)+wt10(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)
								cmat(imat)=cmat(imat)+wt1m(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)
								bmat(imat)=bmat(imat)+wt1p(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)
							end if
						else
							do j=1,mjm2
								imat=l1+mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
								amat(imat)=amat(imat)+wt20(j,ibnd)*xa(j)*tx(j,l3n)
								cmat(imat)=cmat(imat)+wt2m(j,ibnd)*xa(j)*tx(j,l3n)
								bmat(imat)=bmat(imat)+wt2p(j,ibnd)*xa(j)*tx(j,l3n)
							end do
							imat=l1+mnum3*(l2-1+mnum3*mjm2)+nskpn(i)
							if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
								amat(imat)=amat(imat)+wt20(mjm1,2)*xa(mjm1)*tx(mjm1,l3n)
								cmat(imat)=cmat(imat)+wt2m(mjm1,2)*xa(mjm1)*tx(mjm1,l3n)
								bmat(imat)=bmat(imat)+wt2p(mjm1,2)*xa(mjm1)*tx(mjm1,l3n)
!							else if (ivar == 4) then
!								amat(imat)=amat(imat)+wt20(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)
!								cmat(imat)=cmat(imat)+wt2m(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)
!								bmat(imat)=bmat(imat)+wt2p(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)
!								amat(imat)=amat(imat)+dc2p(mjm1)*xa(mjm1)*tx(mjm1,l3n)*(r(mj)-r(mjm2))/ &
!										      (r(mjm1)-r(mjm2))
!								cmat(imat)=cmat(imat)-dc2p(mjm1)*xa(mjm1)*tx(mjm1,l3n)*(r(mj)-r(mjm1))/ &
!										      (r(mjm1)-r(mjm2))
							else
								amat(imat)=amat(imat)+wt20(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)
								cmat(imat)=cmat(imat)+wt2m(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)
								bmat(imat)=bmat(imat)+wt2p(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)
							end if
						end if
					end do
				end do
			end do
		end do
		do i=nc_start,nc_end
			mnum3=noeqn*mnumn(i)
			do l2t=1,mnumn(i)
				l=l2t+lnumn(i-1)
				l2=l2t+(ivar-1)*mnumn(i)
				m=mm(lln(l))
				n=nn(lln(l))
				ibnd=1
				if (isv == 1 .and. m == 0) ibnd=2
				do l1t=1,mnumn(i)
					l1=l1t+(ieqn-1)*mnumn(i)+ilp*mnum3
					lpn=nskp2n(i)+mnumn(i)*(l2t-1)+l1t
					do l3=1,leqmax
						l3n=l3
						if (itz ==  1 .and. itx ==  1) x=coef*cmapp(lpn,l3n)
						if (itz ==  1 .and. itx == -1) x=coef*cmamp(lpn,l3n)
						if (itz == -1 .and. itx == -1) x=coef*cmamm(lpn,l3n)
						if (itz == -1 .and. itx ==  1) x=coef*cmapm(lpn,l3n)
						if (l3 == leq0 .and. itx == -1 .and. l1t == l2t) then
							l1=l1t+(ieqn-1)*mnumn(i)+(1-ilp)*mnum3
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
								imat=l1+2*mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
								amat(imat)=amat(imat)+xa(j)*tx(j,l3n)
							end do
						else if (ir == 1) then
							do j=1,mjm2
								imat=l1+2*mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
								amat(imat)=amat(imat)+wt10(j,ibnd)*xa(j)*tx(j,l3n)
								cmat(imat)=cmat(imat)+wt1m(j,ibnd)*xa(j)*tx(j,l3n)
								bmat(imat)=bmat(imat)+wt1p(j,ibnd)*xa(j)*tx(j,l3n)
							end do
							imat=l1+2*mnum3*(l2-1+mnum3*mjm2)+nskpn(i)
							if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
								amat(imat)=amat(imat)+wt10(mjm1,2)*xa(mjm1)*tx(mjm1,l3n)
								cmat(imat)=cmat(imat)+wt1m(mjm1,2)*xa(mjm1)*tx(mjm1,l3n)
								bmat(imat)=bmat(imat)+wt1p(mjm1,2)*xa(mjm1)*tx(mjm1,l3n)
!							else if (ivar == 4) then
!								amat(imat)=amat(imat)+wt10(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)
!								cmat(imat)=cmat(imat)+wt1m(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)
!								bmat(imat)=bmat(imat)+wt1p(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)
!								amat(imat)=amat(imat)+dc1p(mjm1)*xa(mjm1)*tx(mjm1,l3n)*(r(mj)-r(mjm2))/ &
!										      (r(mjm1)-r(mjm2))
!								cmat(imat)=cmat(imat)-dc1p(mjm1)*xa(mjm1)*tx(mjm1,l3n)*(r(mj)-r(mjm1))/ &
!										      (r(mjm1)-r(mjm2))
							else
								amat(imat)=amat(imat)+wt10(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)
								cmat(imat)=cmat(imat)+wt1m(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)
								bmat(imat)=bmat(imat)+wt1p(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)
							end if
						else
							do j=1,mjm2
								imat=l1+2*mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
								amat(imat)=amat(imat)+wt20(j,ibnd)*xa(j)*tx(j,l3n)
								cmat(imat)=cmat(imat)+wt2m(j,ibnd)*xa(j)*tx(j,l3n)
								bmat(imat)=bmat(imat)+wt2p(j,ibnd)*xa(j)*tx(j,l3n)
							end do
							imat=l1+2*mnum3*(l2-1+mnum3*mjm2)+nskpn(i)
							if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
								amat(imat)=amat(imat)+wt20(mjm1,2)*xa(mjm1)*tx(mjm1,l3n)
								cmat(imat)=cmat(imat)+wt2m(mjm1,2)*xa(mjm1)*tx(mjm1,l3n)
								bmat(imat)=bmat(imat)+wt2p(mjm1,2)*xa(mjm1)*tx(mjm1,l3n)
!							else if (ivar == 4) then
!								amat(imat)=amat(imat)+wt20(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)
!								cmat(imat)=cmat(imat)+wt2m(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)
!								bmat(imat)=bmat(imat)+wt2p(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)
!								amat(imat)=amat(imat)+dc2p(mjm1)*xa(mjm1)*tx(mjm1,l3n)*(r(mj)-r(mjm2))/ &
!										      (r(mjm1)-r(mjm2))
!								cmat(imat)=cmat(imat)-dc2p(mjm1)*xa(mjm1)*tx(mjm1,l3n)*(r(mj)-r(mjm1))/ &
!										      (r(mjm1)-r(mjm2))
							else
								amat(imat)=amat(imat)+wt20(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)
								cmat(imat)=cmat(imat)+wt2m(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)
								bmat(imat)=bmat(imat)+wt2p(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)
							end if
						end if
					end do
				end do
			end do
		end do

	end subroutine blockj

	subroutine block0(tx,ieqn,ivar,ith,ir,izt,coef)

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
		integer :: ity,isv,ips,iph,itz,ilp,ich,i,mnum3,l,m,n,l1,l2,l2t,lp,ibnd,j,imat,il,lpq,lq,lq2,imtq

		ity= 1
		if (ivar == 2 .or. ivar == 4 .or. ivar == 6 .or. ivar == 7 .or. ivar == ivalp .or. ivar > iq) ity=-1
		isv=1
		ips=(-1)**((ith+izt+1)/2)
		iph=ity**(ith+izt)*ips
		itz=(-1)**(ith+izt)
		ilp=(1-itz)/2
		ity=ity*itz
		ich=ity
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
						imat=l1+mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
						amat(imat)=amat(imat)+xa(j)*tx(j)
					end do
				else if(ir == 1) then
					do j=1,mjm2
						imat=l1+mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
						amat(imat)=amat(imat)+wt10(j,ibnd)*xa(j)*tx(j)
						cmat(imat)=cmat(imat)+wt1m(j,ibnd)*xa(j)*tx(j)
						bmat(imat)=bmat(imat)+wt1p(j,ibnd)*xa(j)*tx(j)
					end do
					imat=l1+mnum3*(l2-1+mnum3*mjm2)+nskpn(i)
					if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
						amat(imat)=amat(imat)+wt10(mjm1,2)*xa(mjm1)*tx(mjm1)
						cmat(imat)=cmat(imat)+wt1m(mjm1,2)*xa(mjm1)*tx(mjm1)
						bmat(imat)=bmat(imat)+wt1p(mjm1,2)*xa(mjm1)*tx(mjm1)
!					else if (ivar == 4) then
!						amat(imat)=amat(imat)+wt10(mjm1,1)*xa(mjm1)*tx(mjm1)
!						cmat(imat)=cmat(imat)+wt1m(mjm1,1)*xa(mjm1)*tx(mjm1)
!						bmat(imat)=bmat(imat)+wt1p(mjm1,1)*xa(mjm1)*tx(mjm1)
!						amat(imat)=amat(imat)+dc1p(mjm1)*xa(mjm1)*tx(mjm1)*(r(mj)-r(mjm2))/ &
!								      (r(mjm1)-r(mjm2))
!						cmat(imat)=cmat(imat)-dc1p(mjm1)*xa(mjm1)*tx(mjm1)*(r(mj)-r(mjm1))/ &
!								      (r(mjm1)-r(mjm2))
					else
						amat(imat)=amat(imat)+wt10(mjm1,1)*xa(mjm1)*tx(mjm1)
						cmat(imat)=cmat(imat)+wt1m(mjm1,1)*xa(mjm1)*tx(mjm1)
						bmat(imat)=bmat(imat)+wt1p(mjm1,1)*xa(mjm1)*tx(mjm1)
					end if
				else
					do j=1,mjm2
						imat=l1+mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
						amat(imat)=amat(imat)+wt20(j,ibnd)*xa(j)*tx(j)
						cmat(imat)=cmat(imat)+wt2m(j,ibnd)*xa(j)*tx(j)
						bmat(imat)=bmat(imat)+wt2p(j,ibnd)*xa(j)*tx(j)
					end do
					imat=l1+mnum3*(l2-1+mnum3*mjm2)+nskpn(i)
					if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
						amat(imat)=amat(imat)+wt20(mjm1,2)*xa(mjm1)*tx(mjm1)
						cmat(imat)=cmat(imat)+wt2m(mjm1,2)*xa(mjm1)*tx(mjm1)
						bmat(imat)=bmat(imat)+wt2p(mjm1,2)*xa(mjm1)*tx(mjm1)
!					else if (ivar == 4) then
!						amat(imat)=amat(imat)+wt20(mjm1,1)*xa(mjm1)*tx(mjm1)
!						cmat(imat)=cmat(imat)+wt2m(mjm1,1)*xa(mjm1)*tx(mjm1)
!						bmat(imat)=bmat(imat)+wt2p(mjm1,1)*xa(mjm1)*tx(mjm1)
!						amat(imat)=amat(imat)+dc2p(mjm1)*xa(mjm1)*tx(mjm1)*(r(mj)-r(mjm2))/ &
!								      (r(mjm1)-r(mjm2))
!						cmat(imat)=cmat(imat)-dc2p(mjm1)*xa(mjm1)*tx(mjm1)*(r(mj)-r(mjm1))/ &
!								      (r(mjm1)-r(mjm2))
					else
						amat(imat)=amat(imat)+wt20(mjm1,1)*xa(mjm1)*tx(mjm1)
						cmat(imat)=cmat(imat)+wt2m(mjm1,1)*xa(mjm1)*tx(mjm1)
						bmat(imat)=bmat(imat)+wt2p(mjm1,1)*xa(mjm1)*tx(mjm1)
					end if
				end if
			end do
		end do
		do i=nc_start,nc_end
			mnum3=noeqn*mnumn(i)
			do l2t=1,mnumn(i)
				l=l2t+lnumn(i-1)
				l2=l2t+(ivar-1)*mnumn(i)
				m=mm(lln(l))
				n=nn(lln(l))
				ibnd=1
				if (isv == 1 .and. m == 0) ibnd=2
				l1=l2t+(ieqn-1)*mnumn(i)+ilp*mnum3
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
						imat=l1+2*mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
						amat(imat)=amat(imat)+xa(j)*tx(j)
					end do
				else if (ir == 1) then
					do j=1,mjm2
						imat=l1+2*mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
						amat(imat)=amat(imat)+wt10(j,ibnd)*xa(j)*tx(j)
						cmat(imat)=cmat(imat)+wt1m(j,ibnd)*xa(j)*tx(j)
						bmat(imat)=bmat(imat)+wt1p(j,ibnd)*xa(j)*tx(j)
					end do
					imat=l1+2*mnum3*(l2-1+mnum3*mjm2)+nskpn(i)
					if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
						amat(imat)=amat(imat)+wt10(mjm1,2)*xa(mjm1)*tx(mjm1)
						cmat(imat)=cmat(imat)+wt1m(mjm1,2)*xa(mjm1)*tx(mjm1)
						bmat(imat)=bmat(imat)+wt1p(mjm1,2)*xa(mjm1)*tx(mjm1)
!					else if (ivar == 4) then
!						amat(imat)=amat(imat)+wt10(mjm1,1)*xa(mjm1)*tx(mjm1)
!						cmat(imat)=cmat(imat)+wt1m(mjm1,1)*xa(mjm1)*tx(mjm1)
!						bmat(imat)=bmat(imat)+wt1p(mjm1,1)*xa(mjm1)*tx(mjm1)
!						amat(imat)=amat(imat)+dc1p(mjm1)*xa(mjm1)*tx(mjm1)*(r(mj)-r(mjm2))/ &
!								      (r(mjm1)-r(mjm2))
!						cmat(imat)=cmat(imat)-dc1p(mjm1)*xa(mjm1)*tx(mjm1)*(r(mj)-r(mjm1))/ &
!								      (r(mjm1)-r(mjm2))
					else
						amat(imat)=amat(imat)+wt10(mjm1,1)*xa(mjm1)*tx(mjm1)
						cmat(imat)=cmat(imat)+wt1m(mjm1,1)*xa(mjm1)*tx(mjm1)
						bmat(imat)=bmat(imat)+wt1p(mjm1,1)*xa(mjm1)*tx(mjm1)
					end if
				else
					do j=1,mjm2
						imat=l1+2*mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
						amat(imat)=amat(imat)+wt20(j,ibnd)*xa(j)*tx(j)
						cmat(imat)=cmat(imat)+wt2m(j,ibnd)*xa(j)*tx(j)
						bmat(imat)=bmat(imat)+wt2p(j,ibnd)*xa(j)*tx(j)
					end do
					imat=l1+2*mnum3*(l2-1+mnum3*mjm2)+nskpn(i)
					if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
						amat(imat)=amat(imat)+wt20(mjm1,2)*xa(mjm1)*tx(mjm1)
						cmat(imat)=cmat(imat)+wt2m(mjm1,2)*xa(mjm1)*tx(mjm1)
						bmat(imat)=bmat(imat)+wt2p(mjm1,2)*xa(mjm1)*tx(mjm1)
!					else if (ivar == 4) then
!						amat(imat)=amat(imat)+wt20(mjm1,1)*xa(mjm1)*tx(mjm1)
!						cmat(imat)=cmat(imat)+wt2m(mjm1,1)*xa(mjm1)*tx(mjm1)
!						bmat(imat)=bmat(imat)+wt2p(mjm1,1)*xa(mjm1)*tx(mjm1)
!						amat(imat)=amat(imat)+dc2p(mjm1)*xa(mjm1)*tx(mjm1)*(r(mj)-r(mjm2))/ &
!								      (r(mjm1)-r(mjm2))
!						cmat(imat)=cmat(imat)-dc2p(mjm1)*xa(mjm1)*tx(mjm1)*(r(mj)-r(mjm1))/ &
!								      (r(mjm1)-r(mjm2))
					else
						amat(imat)=amat(imat)+wt20(mjm1,1)*xa(mjm1)*tx(mjm1)
						cmat(imat)=cmat(imat)+wt2m(mjm1,1)*xa(mjm1)*tx(mjm1)
						bmat(imat)=bmat(imat)+wt2p(mjm1,1)*xa(mjm1)*tx(mjm1)
					end if
				end if
			end do
		end do

	end subroutine block0
	
	subroutine block_dlsq(ieqn,ivar,coef)

		use param
		use equil
		implicit none

		integer :: ieqn,ivar
		real(IDP) :: coef

		interface
			subroutine blockj(tx,itx,ieqn,ivar,ith,ir,izt,coef)
				use param
				implicit none
				integer :: itx,ieqn,ivar,ith,ir,izt
				real(IDP) :: coef
				real(IDP), dimension(0:,0:) :: tx
			end subroutine blockj
		end interface

		call blockj(lplrr,1,ieqn,ivar,0,2,0,coef)
		call blockj(lplrt,-1,ieqn,ivar,1,1,0,coef)
		call blockj(lplrz,-1,ieqn,ivar,0,1,1,coef)
		call blockj(lpltt,1,ieqn,ivar,2,0,0,coef)
		call blockj(lpltz,1,ieqn,ivar,1,0,1,coef)
		call blockj(lplzz,1,ieqn,ivar,0,0,2,coef)
		call blockj(lplr,1,ieqn,ivar,0,1,0,coef)
		call blockj(lplt,-1,ieqn,ivar,1,0,0,coef)
		call blockj(lplz,-1,ieqn,ivar,0,0,1,coef)

	end subroutine block_dlsq

	subroutine blockj_landau_grad_parallel(tx,itx,ieqn,ivar,ith,ir,izt,coef)

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
		integer :: ity,isv,ips,iph,itz,icn,ilp,ich,i,mnum3,l,m,n,l1,l2,l3,l1t,l2t,l3n,lp,lpp,lpn,ibnd,j,imat,il

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
								imat=l1+mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
								amat(imat)=amat(imat)+xa(j)*tx(j,l3n)
							end do
						else if (ir == 1) then
							do j=1,mjm1
								imat=l1+mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
								amat(imat)=amat(imat)+wt10(j,ibnd)*xa(j)*tx(j,l3n)
								cmat(imat)=cmat(imat)+wt1m(j,ibnd)*xa(j)*tx(j,l3n)
								bmat(imat)=bmat(imat)+wt1p(j,ibnd)*xa(j)*tx(j,l3n)
							end do
						else
							do j=1,mjm1
								imat=l1+mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
								amat(imat)=amat(imat)+wt20(j,ibnd)*xa(j)*tx(j,l3n)
								cmat(imat)=cmat(imat)+wt2m(j,ibnd)*xa(j)*tx(j,l3n)
								bmat(imat)=bmat(imat)+wt2p(j,ibnd)*xa(j)*tx(j,l3n)
							end do
						end if
					end do
				end do
			end do
		end do
		do i=nc_start,nc_end
			mnum3=noeqn*mnumn(i)
			do l2t=1,mnumn(i)
				l=l2t+lnumn(i-1)
				l2=l2t+(ivar-1)*mnumn(i)
				m=mm(lln(l))
				n=nn(lln(l))
				ibnd=1
				if (isv == 1 .and. m == 0) ibnd=2
				do l1t=1,mnumn(i)
					l1=l1t+(ieqn-1)*mnumn(i)+ilp*mnum3
					lpn=nskp2n(i)+mnumn(i)*(l2t-1)+l1t
					do l3=1,leqmax
						l3n=l3
						if (itz ==  1 .and. itx ==  1) x=coef*cmapp(lpn,l3n)
						if (itz ==  1 .and. itx == -1) x=coef*cmamp(lpn,l3n)
						if (itz == -1 .and. itx == -1) x=coef*cmamm(lpn,l3n)
						if (itz == -1 .and. itx ==  1) x=coef*cmapm(lpn,l3n)
						if (l3 == leq0 .and. itx == -1 .and. l1t == l2t) then
							l1=l1t+(ieqn-1)*mnumn(i)+(1-ilp)*mnum3
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
								imat=l1+2*mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
								amat(imat)=amat(imat)+xa(j)*tx(j,l3n)
							end do
						else if (ir == 1) then
							do j=1,mjm1
								imat=l1+2*mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
								amat(imat)=amat(imat)+wt10(j,ibnd)*xa(j)*tx(j,l3n)
								cmat(imat)=cmat(imat)+wt1m(j,ibnd)*xa(j)*tx(j,l3n)
								bmat(imat)=bmat(imat)+wt1p(j,ibnd)*xa(j)*tx(j,l3n)
							end do
						else
							do j=1,mjm1
								imat=l1+2*mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
								amat(imat)=amat(imat)+wt20(j,ibnd)*xa(j)*tx(j,l3n)
								cmat(imat)=cmat(imat)+wt2m(j,ibnd)*xa(j)*tx(j,l3n)
								bmat(imat)=bmat(imat)+wt2p(j,ibnd)*xa(j)*tx(j,l3n)
							end do
						end if
					end do
				end do
			end do
		end do

	end subroutine blockj_landau_grad_parallel
	
	subroutine block0_dlsq(tx,ieqn,ivar,coef)

		use param
		use var_para
		use cotrol
		use domain
		use dynamo
		implicit none

		integer :: ieqn,ivar
		real(IDP) :: coef
		real(IDP), dimension(0:) :: tx
		integer :: i,l,m,l1,l2,l2t,j,mnum3,imat

		do i=nr_start,nr_end
			mnum3=noeqn*mnumn(i)
			do l2t=1,mnumn(i)
				l=l2t+lnumn(i-1)
				l2=l2t+(ivar-1)*mnumn(i)
				m=mm(lln(l))
				l1=l2t+(ieqn-1)*mnumn(i)
				do j=1,mjm1
					imat=l1+mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
					amat(imat)=amat(imat)-coef*m*m*rinv(j)*rinv(j)*tx(j)
				end do
				imat=l1+mnum3*(l2-1)+nskpn(i)
				if (m == 0) then
					amat(imat)=amat(imat)-coef*tx(1)*(del2cm(1)+del2cp(1)-r(2)**2*del2cm(1)/(r(2)**2-r(1)**2))
					bmat(imat)=bmat(imat)+coef*tx(1)*(del2cp(1)-r(1)**2*del2cm(1)/(r(2)**2-r(1)**2))
				else
					amat(imat)=amat(imat)-coef*tx(1)*(del2cm(1)+del2cp(1))
					bmat(imat)=bmat(imat)+coef*tx(1)*del2cp(1)
				end if
				do j=2,mjm2
					imat=l1+mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
					amat(imat)=amat(imat)-coef*tx(j)*(del2cm(j)+del2cp(j))
					cmat(imat)=cmat(imat)+coef*tx(j)*del2cm(j)
					bmat(imat)=bmat(imat)+coef*tx(j)*del2cp(j)
				end do
				imat=l1+mnum3*(l2-1+mnum3*mjm2)+nskpn(i)
				if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
					amat(imat)=amat(imat)-coef*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1)-del2cp(mjm1)* &
							      (r(mj)-r(mjm2))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2))
					cmat(imat)=cmat(imat)+coef*tx(mjm1)*(del2cm(mjm1)-del2cp(mjm1)* &
							      (r(mj)-r(mjm1))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2))
!				else if (ivar == 4) then
!					amat(imat)=amat(imat)-coef*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1)-del2cp(mjm1)* &
!							      (r(mj)-r(mjm2))/(r(mjm1)-r(mjm2)))
!					cmat(imat)=cmat(imat)+coef*tx(mjm1)*(del2cm(mjm1)-del2cp(mjm1)* &
!							      (r(mj)-r(mjm1))/(r(mjm1)-r(mjm2)))
				else
					amat(imat)=amat(imat)-coef*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1))
					cmat(imat)=cmat(imat)+coef*tx(mjm1)*del2cm(mjm1)
				end if
			end do
		end do
		do i=nc_start,nc_end
			mnum3=noeqn*mnumn(i)
			do l2t=1,mnumn(i)
				l=l2t+lnumn(i-1)
				l2=l2t+(ivar-1)*mnumn(i)
				m=mm(lln(l))
				l1=l2t+(ieqn-1)*mnumn(i)
				do j=1,mjm1
					imat=l1+2*mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
					amat(imat)=amat(imat)-coef*m*m*rinv(j)*rinv(j)*tx(j)
				end do
				imat=l1+2*mnum3*(l2-1)+nskpn(i)
				if (m == 0) then
					amat(imat)=amat(imat)-coef*tx(1)*(del2cm(1)+del2cp(1)-r(2)**2*del2cm(1)/(r(2)**2-r(1)**2))
					bmat(imat)=bmat(imat)+coef*tx(1)*(del2cp(1)-r(1)**2*del2cm(1)/(r(2)**2-r(1)**2))
				else
					amat(imat)=amat(imat)-coef*tx(1)*(del2cm(1)+del2cp(1))
					bmat(imat)=bmat(imat)+coef*tx(1)*del2cp(1)
				end if
				do j=2,mjm2
					imat=l1+2*mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
					amat(imat)=amat(imat)-coef*tx(j)*(del2cm(j)+del2cp(j))
					cmat(imat)=cmat(imat)+coef*tx(j)*del2cm(j)
					bmat(imat)=bmat(imat)+coef*tx(j)*del2cp(j)
				end do
				imat=l1+2*mnum3*(l2-1+mnum3*mjm2)+nskpn(i)
				if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
					amat(imat)=amat(imat)-coef*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1)-del2cp(mjm1)* &
							      (r(mj)-r(mjm2))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2))
					cmat(imat)=cmat(imat)+coef*tx(mjm1)*(del2cm(mjm1)-del2cp(mjm1)* &
							      (r(mj)-r(mjm1))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2))
!				else if (ivar == 4) then
!					amat(imat)=amat(imat)-coef*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1)-del2cp(mjm1)* &
!							      (r(mj)-r(mjm2))/(r(mjm1)-r(mjm2)))
!					cmat(imat)=cmat(imat)+coef*tx(mjm1)*(del2cm(mjm1)-del2cp(mjm1)* &
!							      (r(mj)-r(mjm1))/(r(mjm1)-r(mjm2)))
				else
					amat(imat)=amat(imat)-coef*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1))
					cmat(imat)=cmat(imat)+coef*tx(mjm1)*del2cm(mjm1)
				end if
			end do
		end do

	end subroutine block0_dlsq

	subroutine block0_dlsqr(tx,ieqn,ivar,coef,prf)

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
		integer :: i,l,m,l1,l2,l2t,j,mnum3,imat

		do i=nr_start,nr_end
			mnum3=noeqn*mnumn(i)
			do l2t=1,mnumn(i)
				l=l2t+lnumn(i-1)
				l2=l2t+(ivar-1)*mnumn(i)
				m=mm(lln(l))
				l1=l2t+(ieqn-1)*mnumn(i)
				do j=1,mjm1
					imat=l1+mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
					amat(imat)=amat(imat)-coef(i)*m*m*rinv(j)*rinv(j)*tx(j)*prf(j)
				end do
				imat=l1+mnum3*(l2-1)+nskpn(i)
				if (m == 0) then
					amat(imat)=amat(imat)-coef(i)*tx(1)*(del2cm(1)+del2cp(1)-r(2)**2*del2cm(1)/(r(2)**2-r(1)**2))*prf(1)
					bmat(imat)=bmat(imat)+coef(i)*tx(1)*(del2cp(1)-r(1)**2*del2cm(1)/(r(2)**2-r(1)**2))*prf(1)
				else
					amat(imat)=amat(imat)-coef(i)*tx(1)*(del2cm(1)+del2cp(1))*prf(1)
					bmat(imat)=bmat(imat)+coef(i)*tx(1)*del2cp(1)*prf(1)
				end if
				do j=2,mjm2
					imat=l1+mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
					amat(imat)=amat(imat)-coef(i)*tx(j)*(del2cm(j)+del2cp(j))*prf(j)
					cmat(imat)=cmat(imat)+coef(i)*tx(j)*del2cm(j)*prf(j)
					bmat(imat)=bmat(imat)+coef(i)*tx(j)*del2cp(j)*prf(j)
				end do
				imat=l1+mnum3*(l2-1+mnum3*mjm2)+nskpn(i)
				if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
					amat(imat)=amat(imat)-coef(i)*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1)-del2cp(mjm1)* &
							      (r(mj)-r(mjm2))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2))*prf(mj)
					cmat(imat)=cmat(imat)+coef(i)*tx(mjm1)*(del2cm(mjm1)-del2cp(mjm1)* &
							      (r(mj)-r(mjm1))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2))*prf(mj)
!				else if (ivar == 4) then
!					amat(imat)=amat(imat)-coef(i)*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1)-del2cp(mjm1)* &
!							      (r(mj)-r(mjm2))/(r(mjm1)-r(mjm2)))*prf(mj)
!					cmat(imat)=cmat(imat)+coef(i)*tx(mjm1)*(del2cm(mjm1)-del2cp(mjm1)* &
!							      (r(mj)-r(mjm1))/(r(mjm1)-r(mjm2)))*prf(mj)
				else
					amat(imat)=amat(imat)-coef(i)*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1))*prf(mjm1)
					cmat(imat)=cmat(imat)+coef(i)*tx(mjm1)*del2cm(mjm1)*prf(mjm1)
				end if
			end do
		end do
		do i=nc_start,nc_end
			mnum3=noeqn*mnumn(i)
			do l2t=1,mnumn(i)
				l=l2t+lnumn(i-1)
				l2=l2t+(ivar-1)*mnumn(i)
				m=mm(lln(l))
				l1=l2t+(ieqn-1)*mnumn(i)
				do j=1,mjm1
					imat=l1+2*mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
					amat(imat)=amat(imat)-coef(i)*m*m*rinv(j)*rinv(j)*tx(j)*prf(j)
				end do
				imat=l1+2*mnum3*(l2-1)+nskpn(i)
				if (m == 0) then
					amat(imat)=amat(imat)-coef(i)*tx(1)*(del2cm(1)+del2cp(1)-r(2)**2*del2cm(1)/(r(2)**2-r(1)**2))*prf(1)
					bmat(imat)=bmat(imat)+coef(i)*tx(1)*(del2cp(1)-r(1)**2*del2cm(1)/(r(2)**2-r(1)**2))*prf(1)
				else
					amat(imat)=amat(imat)-coef(i)*tx(1)*(del2cm(1)+del2cp(1))*prf(1)
					bmat(imat)=bmat(imat)+coef(i)*tx(1)*del2cp(1)*prf(1)
				end if
				do j=2,mjm2
					imat=l1+2*mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
					amat(imat)=amat(imat)-coef(i)*tx(j)*(del2cm(j)+del2cp(j))*prf(j)
					cmat(imat)=cmat(imat)+coef(i)*tx(j)*del2cm(j)*prf(j)
					bmat(imat)=bmat(imat)+coef(i)*tx(j)*del2cp(j)*prf(j)
				end do
				imat=l1+2*mnum3*(l2-1+mnum3*mjm2)+nskpn(i)
				if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
					amat(imat)=amat(imat)-coef(i)*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1)-del2cp(mjm1)* &
							      (r(mj)-r(mjm2))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2))*prf(mj)
					cmat(imat)=cmat(imat)+coef(i)*tx(mjm1)*(del2cm(mjm1)-del2cp(mjm1)* &
							      (r(mj)-r(mjm1))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2))*prf(mj)
!				else if (ivar == 4) then
!					amat(imat)=amat(imat)-coef(i)*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1)-del2cp(mjm1)* &
!							      (r(mj)-r(mjm2))/(r(mjm1)-r(mjm2)))*prf(mj)
!					cmat(imat)=cmat(imat)+coef(i)*tx(mjm1)*(del2cm(mjm1)-del2cp(mjm1)* &
!							      (r(mj)-r(mjm1))/(r(mjm1)-r(mjm2)))*prf(mj)
				else
					amat(imat)=amat(imat)-coef(i)*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1))*prf(mjm1)
					cmat(imat)=cmat(imat)+coef(i)*tx(mjm1)*del2cm(mjm1)*prf(mjm1)
				end if
			end do
		end do

	end subroutine block0_dlsqr