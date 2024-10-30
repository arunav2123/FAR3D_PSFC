	subroutine setmod

		use param
		use processor
		use var_para
		use cotrol
		use domain
		use dynamo
		implicit none

		interface
			subroutine inputlist
			end subroutine inputlist	
		end interface	
		
		integer :: l,n,m,lh,mt,nt,ngcd,ndiv,lll,i,lm,lq,mcnt,nl,lp,lpp,lppn,lop,l1,mband,mmp,nnp,np,idel,idelp1,i_mn,ip,mjst,mjnd, &
			   iPE,mxnum,neq,nm
		integer, dimension(1) :: PE_min
		integer, dimension(:), allocatable :: nsort,nmap,nnum_PE,nrsq_PE,nrl_PE,nvalue,iflag
		integer, dimension(:,:), allocatable :: nr_PE

!		Set up the modes distribution and the mode couplings of the model  

		mmin=0
		mmax=0
		nmin=0
		nmax=0
		mmineq=0
		mmaxeq=0
		nmineq=0
		nmaxeq=0
		do l=1,lmax
			m=mm(l)
			n=nn(l)
			signl(l)=1
			if (n < 0 .or. (n == 0 .and. m < 0)) signl(l)=-1
			if (n == 0 .and. m == 0) signl(l)=0
		end do
		mmax=maxval(mm(1:lmax))
		mmin=minval(mm(1:lmax))
		nmax=maxval(nn(1:lmax))
		nmin=minval(nn(1:lmax))
		mmaxx=max(mmax,abs(mmin))		
		do l=1,leqmax
			m=mmeq(l)
			n=nneq(l)
			sgnleq(l)=1
			if (n < 0 .or. (n == 0 .and. m < 0)) sgnleq(l)=-1
			if (n == 0 .and. m == 0) sgnleq(l)=0
		end do
		mmaxeq=maxval(mmeq(1:leqmax))
		mmineq=minval(mmeq(1:leqmax))
		nmaxeq=maxval(nneq(1:leqmax))
		nmineq=minval(nneq(1:leqmax))
		mmaxxeq=max(mmaxeq,abs(mmineq))		

		allocate (ll(mmin:mmax,nmin:nmax))

		ll=0
		do l=1,lmax
			ll(mm(l),nn(l))=l
		end do
		l0=0
		if (0 >= mmin .and. 0 <= mmax .and. 0 >= nmin .and. 0 <= nmax) l0=ll(0,0)
		if (l0 == 0) then
			if (myPE == 0) write (6,'("  setmod: l0=0")')
			stop
		end if

		allocate (lleq(mmineq:mmaxeq,nmineq:nmaxeq))

		lleq=0
		do l=1,leqmax
			lleq(mmeq(l),nneq(l))=l
		end do
		leq0=0
		if (0 >= mmineq .and. 0 <= mmaxeq .and. 0 >= nmineq .and. 0 <= nmaxeq) leq0=lleq(0,0)
		if (leq0 == 0) then
			if (myPE == 0) write (6,'("  setmod: leq0=0")')
			stop
		end if

!  find prime harmonics.

		lh=0
		mnxxx=1
		do l=1,lmax

			mt=mm(l)*signl(l)
			nt=nn(l)*signl(l)

			if (mt == 0 .or. nt == 0) then
				if (nt /= 0) then
					mnxxx=max(mnxxx,nt)
					nt=1
				end if
				if (mt /= 0) then
					mnxxx=max(mnxxx,mt)
					mt=1
				end if
			else
				ngcd=1
				do ndiv=2,nt
					if(mt /= ndiv*(mt/ndiv)) cycle
					if(nt /= ndiv*(nt/ndiv)) cycle
					ngcd=ndiv
				end do
				mnxxx=max(mnxxx,ngcd)
				mt=mt/ngcd
				nt=nt/ngcd
			end if

			do lll=1,lh
				if (mt == mh(lll) .and. nt == nh(lll)) exit
			end do
			if (lh == 0 .or. lll > lh) then
				lh=lh+1
				mh(lh)=mt
				nh(lh)=nt
			end if
		end do
		lhmax=lh

		if (myPE == 0) then
			write(6,'(/" mnxxx = ",i3)') mnxxx
		end if

		lh=0
		do l=1,leqmax

			mt=mmeq(l)*sgnleq(l)
			nt=nneq(l)*sgnleq(l)

			if (mt == 0 .or. nt == 0) then
				if (nt /= 0) nt=1
				if (mt /= 0) mt=1
			else
				ngcd=1
				do ndiv=2,nt
					if(mt /= ndiv*(mt/ndiv)) cycle
					if(nt /= ndiv*(nt/ndiv)) cycle
					ngcd=ndiv
				end do
				mt=mt/ngcd
				nt=nt/ngcd
			end if

			do lll=1,lh
				if (mt == mheq(lll) .and. nt == nheq(lll)) exit
			end do
			if (lh == 0 .or. lll > lh) then
				lh=lh+1
				mheq(lh)=mt
				nheq(lh)=nt
			end if
		end do
		lheqmx=lh

		allocate (m1n(nmin:nmax))
		allocate (mrang(nmin:nmax))

		nfp=nmaxeq
		if (nmaxeq > 0) then
			do l=1,leqmax
				if (nneq(l) > 0) nfp=min(nfp,nneq(l))
			end do
			do l=1,leqmax
				if (nneq(l) /= nfp*(nneq(l)/nfp)) then
					write (6,'("  setmod: nneq(",i2,") =",i4," is not multiple of nfp =",i3)') l,nneq(l),nfp
					stop
				end if
			end do
		else
			nfp=1
		end if

		do n=nmin,nmax
			m1n(n)=mmax+1
			mrang(n)=mmin-1
			do l=1,lmax
				if (nn(l) /= n) cycle
				m1n(n)=min(m1n(n),mm(l))
				mrang(n)=max(mrang(n),mm(l))
			end do
			if (mrang(n) == mmin-1) then
				mrang(n)=0
			else
				mrang(n)=mrang(n)-m1n(n)+1
			endif
		end do
		if (nmin < 0 .and. nmax > 0) then
			if (nmin /= -nmax) then
				if (myPE == 0) write (6,'("  setmod: nmin =",i4," nmax =",i3)') nmin,nmax
				stop
			end if			
			do n=1,nmax
				if (n == (n/nfp)*nfp) cycle			
				if (mrang(n)+mrang(-n) > 0 .and. (mrang(-n) /= mrang(n) .or. m1n(n) /= -(m1n(-n)+mrang(-n)-1))) then
					if (myPE == 0) write (6,'("  setmod: n =",i4," mrang =",i3,i4," m1n =",i3,i5)') &
								   n,mrang(n),mrang(-n),m1n(n),m1n(-n)+mrang(-n)-1
					stop
				end if
			end do
		end if

		if (nocpl == 0) then

			lmax0=0
			do i=1,mrang(0)
				m=m1n(0)+i-1
				l=ll(m,0)
				if (l == 0) cycle
				lmax0=lmax0+1
			end do

			allocate (ll0(lmax0))

			lm=0
			do i=1,mrang(0)
				m=m1n(0)+i-1
				l=ll(m,0)
				if (l == 0) cycle
				lm=lm+1
				ll0(lm)=l
			end do
			
		else
		
			lmax0=1
			allocate (ll0(lmax0))
			ll0(1)=l0
		
		end if

		lmaxn=lmax-leqmax
		if (m0dy < 0) m0dy=0
		if (nonlin == 0) m0dy=0
		lmaxn=lmaxn+m0dy
		lmx=noeqn*lmaxn

		allocate (lln(lmaxn))

		nnum=0

		if (nocpl == 0) then

			do n=nmin,nmax
				if (mrang(n) == 0) cycle
				if (n < 0 .and. mrang(-n) /= 0) cycle
				nnum=nnum+1
				mcnt=0
				do i=1,mrang(n)
					m=m1n(n)+i-1
					l=ll(m,n)
					lq=0
					if (n >= nmineq .and. n <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,n)
					if (l == 0 .or. lq > m0dy) cycle
					mcnt=mcnt+1
				end do
				if (-n >= nmin .and. -n <= nmax .and. n /= 0) then
					do i=1,mrang(-n)
						m=m1n(-n)+i-1
						l=ll(m,-n)
						lq=0
						if (-n >= nmineq .and. -n <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,-n)
						if (l == 0 .or. lq > m0dy) cycle
						mcnt=mcnt+1
					end do
				end if
				if (mcnt == 0) nnum=nnum-1
			end do

			allocate (mnumn(nnum))
			allocate (lnumn(0:nnum))
			allocate (nsort(nnum))
			allocate (nvalue(nnum))

			nnum=0

			do n=nmin,nmax
				if (mrang(n) == 0) cycle
				if (n < 0 .and. mrang(-n) /= 0) cycle
				nnum=nnum+1
				mcnt=0
				do i=1,mrang(n)
					m=m1n(n)+i-1
					l=ll(m,n)
					lq=0
					if (n >= nmineq .and. n <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,n)
					if (l == 0 .or. lq > m0dy) cycle
					mcnt=mcnt+1
				end do
				if (-n >= nmin .and. -n <= nmax .and. n /= 0) then
					do i=1,mrang(-n)
						m=m1n(-n)+i-1
						l=ll(m,-n)
						lq=0
						if (-n >= nmineq .and. -n <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,-n)
						if (l == 0 .or. lq > m0dy) cycle
						mcnt=mcnt+1
					end do
				end if
				if (mcnt == 0) then 
					nnum=nnum-1
				else
					nvalue(nnum)=n
					mnumn(nnum)=mcnt
				end if
			end do

			do n=1,nnum
				nsort(n)=n
			end do
			do n=2,nnum
				nl=n
				do while (mnumn(nsort(nl)) < mnumn(nsort(nl-1)))
					np=nsort(nl)
					nsort(nl)=nsort(nl-1)
					nsort(nl-1)=np
					nl=nl-1
					if (nl == 1) exit
				end do
			end do
				
			nl=nnum/(2*numPEs)
			mxnum=2*nl+mod(nnum,2*numPEs)
			allocate (nnum_PE(numPEs))
			allocate (nrsq_PE(numPEs))
			allocate (nrl_PE(numPEs))
			allocate (nr_PE(mxnum,numPEs))
			nnum_PE=0
			nrsq_PE=0
			nrl_PE=0
			n=nnum
			if (nl > 0) then
				do i=1,nl
					do iPE=numPEs,1,-1
						nr_PE(2*i-1,iPE)=nsort(n)
						nrsq_PE(iPE)=nrsq_PE(iPE)+mnumn(nsort(n))*mnumn(nsort(n))
						nrl_PE(iPE)=nrl_PE(iPE)+mnumn(nsort(n))
						n=n-1
					end do
					do iPE=1,numPEs
						nr_PE(2*i,iPE)=nsort(n)
						nrsq_PE(iPE)=nrsq_PE(iPE)+mnumn(nsort(n))*mnumn(nsort(n))
						nrl_PE(iPE)=nrl_PE(iPE)+mnumn(nsort(n))
						n=n-1
					end do
				end do
				nnum_PE=2*nl
				do while (n > 0)
					PE_min=minloc(nrsq_PE)
					iPE=PE_min(1)
					nnum_PE(iPE)=nnum_PE(iPE)+1
					nr_PE(nnum_PE(iPE),iPE)=nsort(n)
					nrsq_PE(iPE)=nrsq_PE(iPE)+mnumn(nsort(n))*mnumn(nsort(n))
					nrl_PE(iPE)=nrl_PE(iPE)+mnumn(nsort(n))
					n=n-1
				end do
			else if(nnum < numPEs) then
				iPE=numPEs
				do i=nnum,1,-1
					nnum_PE(iPE)=1
					nr_PE(1,iPE)=nsort(i)
					nrsq_PE(iPE)=mnumn(nsort(i))*mnumn(nsort(i))
					nrl_PE(iPE)=mnumn(nsort(i))
					iPE=iPE-1
				end do
			else
				do iPE=numPEs,1,-1
					nr_PE(1,iPE)=nsort(n)
					nrsq_PE(iPE)=mnumn(nsort(n))*mnumn(nsort(n))
					nrl_PE(iPE)=mnumn(nsort(n))
					n=n-1
				end do
				nnum_PE=1
				do while (n > 0)
					PE_min=minloc(nrsq_PE)
					iPE=PE_min(1)
					nnum_PE(iPE)=nnum_PE(iPE)+1
					nr_PE(nnum_PE(iPE),iPE)=nsort(n)
					nrsq_PE(iPE)=nrsq_PE(iPE)+mnumn(nsort(n))*mnumn(nsort(n))
					nrl_PE(iPE)=nrl_PE(iPE)+mnumn(nsort(n))
					n=n-1
				end do
			end if
		
			do iPE=1,numPEs
				do i=1,nnum_PE(iPE)
					if (nvalue(nr_PE(i,iPE)) == 0) exit
				end do
				if (i <= nnum_PE(iPE)) exit
			end do
			if (iPE > numPEs) then
				if (myPE == 0) write (6,'("  setmod: No modes with n = 0")')
				iPE=1
			end if
			if (iPE > 1) then
				do i=1,max(nnum_PE(1),nnum_PE(iPE))
					n=nr_PE(i,1)
					nr_PE(i,1)=nr_PE(i,iPE)
					nr_PE(i,iPE)=n
				end do
				n=nnum_PE(1)
				nnum_PE(1)=nnum_PE(iPE)
				nnum_PE(iPE)=n
			end if

			allocate (nmap(nnum))
			allocate (n_st(0:numPEsm1))
			allocate (n_nd(0:numPEsm1))
			n=1
			n_st(0)=1
			n_nd(0)=nnum_PE(1)
			do i=nnum_PE(1),1,-1
				nmap(n)=nvalue(nr_PE(i,1))
				n=n+1
			end do
			do iPE=2,numPEs
				n_st(iPE-1)=n_nd(iPE-2)+1
				n_nd(iPE-1)=n_nd(iPE-2)+nnum_PE(iPE)
				do i=nnum_PE(iPE),1,-1
					nmap(n)=nvalue(nr_PE(i,iPE))
					n=n+1
				end do
			end do
		
			if (myPE == 0) then
				write(6,'(/,"n_st=",20i5)') (n_st(iPE),iPE=0,numPEsm1)
				write(6,'("n_nd=",20i5)') (n_nd(iPE),iPE=0,numPEsm1)
				write(6,'("nmap=",20i5)') (nmap(n),n=1,nnum)
			end if

			nnd=0
			lnumn(0)=0
			nl=0
			lp=0
			do np=1,nnum
				n=nmap(np)
				if (mrang(n) == 0) cycle
				if (n < 0 .and. mrang(-n) /= 0) cycle
				nl=nl+1
				lpp=0
				lppn=0
				do i=1,mrang(n)
					m=m1n(n)+i-1
					l=ll(m,n)
					lq=0
					if (n >= nmineq .and. n <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,n)
					if (l == 0 .or. lq > m0dy) cycle
					lp=lp+1
					lln(lp)=l
					lpp=lpp+1
				end do
				if (-n >= nmin .and. -n <= nmax .and. n /= 0) then
					do i=mrang(-n),1,-1
						m=m1n(-n)+i-1
						l=ll(m,-n)
						lq=0
						if (-n >= nmineq .and. -n <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,-n)
						if (l == 0 .or. lq > m0dy) cycle
						lp=lp+1
						lln(lp)=l
						lppn=lppn+1
					end do
				end if
				if (lpp > 0) then
					mnumn(nl)=lpp
					lnumn(nl)=lp
					if (nmin == 0 .or. nmax == 0 .or. n == 0) then
						nnd=nnd+1
					else if (lppn /= lpp) then
						if (myPE == 0) write (6,'("  setmod: n =",i4," mnumn =",i3,i4)') n,lpp,lppn
						stop
					end if
				else
					nl=nl-1
				end if
			end do
			nst=nnd+1

			if (myPE == 0) write (6,'(/"lln="/("    ",25i5))') (lln(l),l=1,lmaxn)

		else

			do n=nmin,nmax
				if (mrang(n) == 0) cycle
				do i=1,mrang(n)
					m=m1n(n)+i-1
					l=ll(m,n)
					if(l == 0) cycle
					lop=0
					if (-m >= mmin .and. -m <= mmax .and. -n >= nmin .and. -n <= nmax) lop=ll(-m,-n)
					if (signl(l) < 0 .and. lop /= 0) cycle
					lq=0
					if (n >= nmineq .and. n <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,n)
					if (lq > m0dy) cycle
					nnum=nnum+1
					mcnt=1
					if (lop /= 0 .and. (m /= 0 .or. n /= 0)) mcnt=mcnt+1
				end do
			end do

			allocate (mnumn(nnum))
			allocate (lnumn(0:nnum))

			nnd=0
			lnumn(0)=0
			nl=0
			lp=0
			do n=nmin,nmax
				if (mrang(n) == 0) cycle
				do i=1,mrang(n)
					m=m1n(n)+i-1
					l=ll(m,n)
					if (l == 0) cycle
					lop=0
					if (-m >= mmin .and. -m <= mmax .and. -n >= nmin .and. -n <= nmax) lop=ll(-m,-n)
					if (signl(l) < 0 .and. lop /= 0) cycle
					lq=0
					if (n >= nmineq .and. n <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,n)
					if (lq > m0dy) cycle
					nl=nl+1
					lp=lp+1
					lpp=1
					lln(lp)=l
					if (lop /= 0 .and. (m /= 0 .or. n /= 0)) then
						lp=lp+1
						lln(lp)=lop
					else
						nnd=nnd+1
					end if
					mnumn(nl)=lpp
					lnumn(nl)=lp
				end do
			end do
			nst=nnd+1

		end if

		allocate (lo(lmaxn))

		lo=0
		do l=1,lmaxn
			m=-mm(lln(l))
			n=-nn(lln(l))
			if (m < mmin .or. m > mmax .or. n < nmin .or. n > nmax) cycle
			do l1=1,lmaxn
				if (mm(lln(l1)) == m .and. nn(lln(l1)) == n) exit
			end do
			if (l1 > lmaxn) then
				lo(l)=0
			else
				lo(l)=l1
			end if
		end do

		mband=0
		do m=0,mmax
			if (ll(m,0) /= 0) mband=mband+1
		end do
		mxmband=2*mband-1
		do n=1,nmax
			mband=0
			do m=mmin,mmax
				if (ll(m,n) /= 0) mband=mband+1
			end do
			if (mband > mxmband) mxmband=mband
		end do

		mband=0
		do m=0,mmaxeq
			if (lleq(m,0) /= 0) mband=mband+1
		end do
		mbandeq=2*mband-1
		do n=1,nmaxeq
			mband=0
			do m=mmineq,mmaxeq
				if (lleq(m,n) /= 0) mband=mband+1
			end do
			if (mband > mbandeq) mbandeq=mband
		end do

	end subroutine setmod