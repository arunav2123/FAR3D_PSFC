	subroutine linstart

		use mpi
		use param
		use processor
		use var_para
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		integer :: i,j,l,m,n,ibnd,mcnt2,nbdim,nbxdim,nbidim,it1,it2,le,lp,ln,l1,l2,l3,l1t,l2t,l3t,mnum3,mnum6,loca,loci, &
			   ier,it,imat,imt1,imt2,nend,lskp,it3,is1,is2,is3,m1,m2,m3,mm3,mp3,n1,n2,n3,np3,nm3,sgn12,l3p, &
			   i1,ip,ierr,iPE,tag,tagst,ind,lsk,mnum2,lpn,lq
		integer, dimension(MPI_STATUS_SIZE) :: status
		real(IDP) :: epsq,oneos,betfc,betfc_f,betfc_alp,coef,omcyd,omcydalp,beteom,betiom,rmrbar,acoef,bcoef,x,aux
		real(IDP), dimension(:), allocatable :: xa
		real(IDP), dimension(:,:), allocatable :: scp2,scp3
		character*1 :: tb
		
		interface
			subroutine clgam(rgam,na,nb,nc,c1,c2)
				use param
				implicit none
				integer :: na,nb,nc
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: rgam
			end subroutine clgam
			subroutine om(tx,itx,ieqn,ivar,ith,ir,izt,coef)
				use param
				implicit none
				integer :: itx,ieqn,ivar,ith,ir,izt
				real(IDP) :: coef
				real(IDP), dimension(0:,0:) :: tx
			end subroutine om
			subroutine blockj(tx,itx,ieqn,ivar,ith,ir,izt,coef)
				use param
				implicit none
				integer :: itx,ieqn,ivar,ith,ir,izt
				real(IDP) :: coef
				real(IDP), dimension(0:,0:) :: tx
			end subroutine blockj
			subroutine blockj_landau_grad_parallel(tx,itx,ieqn,ivar,ith,ir,izt,coef)
				use param
				implicit none
				integer :: itx,ieqn,ivar,ith,ir,izt
				real(IDP) :: coef
				real(IDP), dimension(0:,0:) :: tx
			end subroutine blockj_landau_grad_parallel
			subroutine om0(tx,ieqn,ivar,ith,ir,izt,coef)
				use param
				implicit none
				integer :: ieqn,ivar,ith,ir,izt
				real(IDP) :: coef
				real(IDP), dimension(0:) :: tx
			end subroutine om0
			subroutine block0(tx,ieqn,ivar,ith,ir,izt,coef)
				use param
				implicit none
				integer :: ieqn,ivar,ith,ir,izt
				real(IDP) :: coef
				real(IDP), dimension(0:) :: tx
			end subroutine block0
			subroutine block_dlsq(ieqn,ivar,coef)
				use param
				implicit none
				integer :: ieqn,ivar
				real(IDP) :: coef
			end subroutine block_dlsq
			subroutine block0_dlsq(tx,ieqn,ivar,coef)
				use param
				implicit none
				integer :: ieqn,ivar
				real(IDP) :: coef
				real(IDP), dimension(0:) :: tx
			end subroutine block0_dlsq
			subroutine block0_dlsqr(tx,ieqn,ivar,coef,prf)
				use param
				implicit none
				integer :: ieqn,ivar
		        real(IDP), dimension(:) :: coef
		        real(IDP), dimension(0:) :: prf
				real(IDP), dimension(0:) :: tx
			end subroutine block0_dlsqr			
			subroutine dbydtheq(d,a,ltype,c1,c2,k)
				use param
				implicit none
				integer :: k,ltype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine dbydtheq
			subroutine grpareq(d,a,ltype,c1,c2)
				use param
				implicit none
				integer :: ltype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine grpareq
			subroutine dbydreq(d,a,c1,c2,k)
				use param
				implicit none
				integer :: k
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine dbydreq
			subroutine dbydzteq(d,a,ltype,c1,c2)
				use param
				implicit none
				integer :: ltype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine dbydzteq
			subroutine dbydr0(d,a,c1,c2,k)
				use param
				implicit none
				integer :: k
				real(IDP) :: c1,c2
				real(IDP), dimension(0:) :: a,d
			end subroutine dbydr0
			subroutine d2bydr20(d,a,c1,c2,k)
				use param
				implicit none
				integer :: k
				real(IDP) :: c1,c2
				real(IDP), dimension(0:) :: a,d
			end subroutine d2bydr20
			subroutine dlstar(ss,ff,itypf,wk1,wk2,wk3,c1,c2)
				use param
				use var_para
				implicit none
				integer :: itypf
				real(IDP) :: c1,c2
				real(IDP), dimension(mj_start:,0:) :: ss,ff,wk1,wk2,wk3
			end subroutine dlstar
			subroutine dlsq(ss,ff,itypf,wk1,wk2,wk3,c1,c2)
				use param
				use var_para
				implicit none
				integer :: itypf
				real(IDP) :: c1,c2
				real(IDP), dimension(mj_start:,0:) :: ss,ff,wk1,wk2,wk3
			end subroutine dlsq
			subroutine dbydth_par(d,a,ltype,c1,c2,k)
				use param
				use var_para
				implicit none
				integer :: k,ltype
				real(IDP) :: c1,c2
				real(IDP), dimension(mj_start:,0:) :: a,d
			end subroutine dbydth_par
			subroutine dbydzt_par(d,a,ltype,c1,c2)
				use param
				use var_para
				implicit none
				integer :: ltype
				real(IDP) :: c1,c2
				real(IDP), dimension(mj_start:,0:) :: a,d
			end subroutine dbydzt_par
			subroutine trnsfr(f,ftype,idir)
				use param
				use var_para
				implicit none
				integer :: ftype,idir
				real(IDP), dimension(mj_start:,0:) :: f
			end subroutine trnsfr
			subroutine decbt(m,n,a,b,c,ip,ier)
				use param
				implicit none
				integer :: m,n,ier
				integer, dimension(m,n) :: ip
				real(IDP), dimension(m,m,n) :: a, b, c
			end subroutine decbt
			subroutine cdecbt(m,n,a,b,c,ip,ier)
				use param
				implicit none
				integer :: m,n,ier
				integer, dimension(m,n) :: ip
				real(IDP), dimension(2*m,m,n) :: a, b, c
			end subroutine cdecbt
			subroutine solbt(m,n,a,b,c,y,ip)
				use param
				implicit none
				integer :: m,n
				integer, dimension(m,n) :: ip
				real(IDP), dimension(m,m,n) :: a, b, c
				real(IDP), dimension(m,n) :: y
			end subroutine solbt
			subroutine csolbt(m,n,a,b,c,y,ip)
				use param
				implicit none
				integer :: m,n
				integer, dimension(m,n) :: ip
				real(IDP), dimension(2*m,m,n) :: a, b, c
				real(IDP), dimension(2*m,n) :: y
			end subroutine csolbt
		end interface

		tb = char(9)
		epsq=eps*eps
		oneos=0.0_IDP
		if (s > 0.0_IDP) oneos=1.0_IDP/s
		betfc=bet0/(2.*epsq)
		betfc_f=LcA2*bet0_f/(2.*epsq)
		betfc_alp=LcA2alp*bet0_alp/(2.*epsq)
		omcyd = omcy
		omcydalp = omcyalp
		beteom=dpres*bet0/(2*epsq*omcyd)
		betiom=(1-dpres)*bet0/(2*epsq*omcyd)

		allocate (stdifpr(nnum))
		allocate (stdifur(nnum))
		allocate (stdifnfr(nnum))
		allocate (stdifvfr(nnum))
		allocate (stdifvr(nnum))				
               
		do i=1,nnum
			l=lnumn(i-1)+1
			n=nn(lln(l))
			ann=abs(n)
!			fctr=1.+100.*(1.+tanh(ann-10.0_IDP))
			fctr=1.0_IDP
			stdifpr(i)=stdifp*fctr
			stdifur(i)=stdifu*fctr
			stdifnfr(i)=stdifnf*fctr
			stdifvfr(i)=stdifvf*fctr		
			stdifvr(i)=stdifv*fctr		
			if (lln(l) == l0) then
				i0=i
			    stdifpr(i)=stdifp
			    stdifur(i)=stdifu
			    stdifnfr(i)=stdifnf
			    stdifvfr(i)=stdifvf		
			    stdifvr(i)=stdifv
			end if
		end do
                aux=1.0_IDP
                do j=0,mj
!                       sd7(j)=(1.+tanh(-100*(0.9-(j*aux/mj)))+0.002)/((1+tanh(-100*0.9))+0.002)
                        sd7(j)=1.0
                end do

!	edge-factor for vfova, nfeq, denseq and preq
!		if (nstep1 == 0) then
!			vfova = vfova*(1.0-exp(-(1.0-r)**2/1.e-3))
!			acoef=nfeq(925)-37.5*(nfeq(924)-nfeq(926))/6.0
!			bcoef=500.0*(nfeq(924)-nfeq(926))/(6*0.075**5)
!			do j=926,mj
!				nfeq(j)=acoef+bcoef*(1.0-r(j))**6
!			end do
!			acoef=denseq(925)-37.5*(denseq(924)-denseq(926))/6.0
!			bcoef=500.0*(denseq(924)-denseq(926))/(6*0.075**5)
!			do j=926,mj
!				denseq(j)=acoef+bcoef*(1.0-r(j))**6
!			end do
!			acoef=preq(925)-37.5*(preq(924)-preq(926))/6.0
!			bcoef=500.0*(preq(924)-preq(926))/(6*0.075**5)
!			do j=926,mj
!				preq(j)=acoef+bcoef*(1.0-r(j))**6
!			end do
!		end if

		call dbydr0(qqinvp,qqinv,0.0_IDP,1.0_IDP,0)
		call dbydr0(denseqr,denseq,0.0_IDP,1.0_IDP,0)
		call dbydr0(dnfeqdr,nfeq,0.0_IDP,1.0_IDP,0)

!  define viscosities

!		sd7=1.0_IDP+49.5*(1.0+tanh((r-0.90)/0.01))
!		sd7=1.0
!		call dbydr0(sd6,sd7,0.0_IDP,1.0_IDP,1)
!		do l=1,leqmax
!			lplrr(:,l)=sd7*lplrr(:,l)
!			lplrt(:,l)=sd7*lplrt(:,l)
!			lpltt(:,l)=sd7*lpltt(:,l)
!			lplr(:,l)=sd7*lplr(:,l)+sd6*gttoj(:,l)
!			lplt(:,l)=sd7*lplt(:,l)-sd6*grtoj(:,l)
!		end do
!		lplrr=0.0
!		lplrt=0.0
!		lpltt=0.0
!		lplr=0.0
!		lplt=0.0
!		lplrr(:,leq0)=sd7
!		lpltt(:,leq0)=sd7
!		lplr(:,leq0)=rinv*sd7+sd6

		vfova2 = vfova*vfova
		if (alpha_on == 1) then
			valphaova2=valphaova*valphaova
			call dbydr0(dnalpeqdr,nalpeq,0.0_IDP,1.0_IDP,0)
		end if

		allocate (nr_st(0:numPEsm1))
		allocate (nr_nd(0:numPEsm1))
		allocate (nc_st(0:numPEsm1))
		allocate (nc_nd(0:numPEsm1))
		allocate (lmaxPE(0:numPEsm1))
		
		do ip=0,numPEsm1
			nr_st(ip)=1
			nr_nd(ip)=0
			nc_st(ip)=nnum+1
			nc_nd(ip)=nnum
		end do
		if (nnd == 0) then
			do ip=0,numPEsm1
				nc_st(ip)=n_st(ip)
				nc_nd(ip)=n_nd(ip)
			end do
		else if (nnd == nnum) then
			do ip=0,numPEsm1
				nr_st(ip)=n_st(ip)
				nr_nd(ip)=n_nd(ip)
			end do
		else
			do ip=0,numPEsm1
				if (n_st(ip) < nst) then
					nr_st(ip)=n_st(ip)
					if (n_nd(ip) > nnd) then
						nr_nd(ip)=nnd
						nc_st(ip)=nst
						nc_nd(ip)=n_nd(ip)
					else
						nr_nd(ip)=n_nd(ip)
					end if
				else
					nc_st(ip)=n_st(ip)
					nc_nd(ip)=n_nd(ip)
				end if
			end do
		end if

		allocate (nskpn(nnum))
		allocate (nskp2n(nnum))
		allocate (nskpxn(nnum))
		allocate (nskpin(nnum))

		nskp2n(1)=0
		nskpxn(1)=0
		nend=min(nnd,nnum-1)
		do i=1,nend
			nskpxn(i+1)=nskpxn(i)+mnumn(i)*noeqn*mj
			nskp2n(i+1)=nskp2n(i)+mnumn(i)*mnumn(i)
		end do
		do i=nst,nnum-1
			nskpxn(i+1)=nskpxn(i)+2*mnumn(i)*noeqn*mj
			nskp2n(i+1)=nskp2n(i)+mnumn(i)*mnumn(i)
		end do
		nbxdim=lmaxn*noeqn*mj
		mcnt2=nskp2n(nnum)+mnumn(nnum)*mnumn(nnum)

		do ip=0,numPEsm1
			nend=min(nr_nd(ip),nc_nd(ip)-1)
			if (nend > 0) then
				nskpn(nr_st(ip))=0
				nskpin(nr_st(ip))=0
			else
				nskpn(nc_st(ip))=0
				nskpin(nc_st(ip))=0
			end if
			do i=nr_st(ip),nend
				nskpn(i+1)=nskpn(i)+mnumn(i)*mnumn(i)*noeqn*noeqn*mjm1
				nskpin(i+1)=nskpin(i)+mnumn(i)*noeqn*mjm1
			end do
			do i=nc_st(ip),nc_nd(ip)-1
				nskpn(i+1)=nskpn(i)+2*mnumn(i)*mnumn(i)*noeqn*noeqn*mjm1
				nskpin(i+1)=nskpin(i)+mnumn(i)*noeqn*mjm1
			end do
		end do

		if (nc_st(myPE) > nnum) then
			nbdim=nskpn(nr_nd(myPE))+mnumn(nr_nd(myPE))*mnumn(nr_nd(myPE))*noeqn*noeqn*mjm1
			nbidim=nskpin(nr_nd(myPE))+mnumn(nr_nd(myPE))*noeqn*mjm1
		else
			nbdim=nskpn(nc_nd(myPE))+2*mnumn(nc_nd(myPE))*mnumn(nc_nd(myPE))*noeqn*noeqn*mjm1
			nbidim=nskpin(nc_nd(myPE))+mnumn(nc_nd(myPE))*noeqn*mjm1
		end if

		if (myPE == 0) then
			write (6,'("mcnt2=",i12)') mcnt2
			write (6,'("nbdim=",i12)') nbdim
			write (6,'("nbidim=",i12)') nbidim
			write (6,'("nbxdim=",i9)') nbxdim
		else
			write (6,'("iPE=",i3," nbdim=",i12)') myPE,nbdim
			write (6,'("iPE=",i3," nbidim=",i12)') myPE,nbidim
		end if
		
		nr_start=nr_st(myPE)
		nr_end=nr_nd(myPE)
		nc_start=nc_st(myPE)
		nc_end=nc_nd(myPE)

		lmaxPE=0
		do ip=0,numPEsm1
			do i=nr_st(ip),nr_nd(ip)
				lmaxPE(ip)=lmaxPE(ip)+mnumn(i)
			end do
			do i=nc_st(ip),nc_nd(ip)
				lmaxPE(ip)=lmaxPE(ip)+2*mnumn(i)
			end do
		end do
		
		allocate (scp(lmaxPE(myPE),0:mj))
		allocate (scu(maxval(lmaxPE)*maxval(mj_dl),0:numPEsm1))
		
		if (myPE == 0) then
			write(6,'(/"nnd =",i4," nst =",i4/)') nnd,nst
			do ip=0,numPEsm1
				write(6,'("iPE =",i3," nr =",2i6," nc =",2i6)') ip,nr_st(ip),nr_nd(ip),nc_st(ip),nc_nd(ip)
			end do
			do ip=0,numPEsm1
				write(6,'("iPE =",i3," lmaxPE =",i4," mj_br =",i4," mj_inc =",i4," mj_st =",i4," mj_dl =",i4)') &
					 ip,lmaxPE(ip),mj_br(ip),mj_inc(ip),mj_st(ip),mj_dl(ip)
			end do
		end if

		if (myPE < numPEsm1) then
			allocate (cmamm(nskp2n(n_st(myPE))+1:nskp2n(n_nd(myPE)+1),leqdim))
			allocate (cmamp(nskp2n(n_st(myPE))+1:nskp2n(n_nd(myPE)+1),leqdim))
			allocate (cmapm(nskp2n(n_st(myPE))+1:nskp2n(n_nd(myPE)+1),leqdim))
			allocate (cmapp(nskp2n(n_st(myPE))+1:nskp2n(n_nd(myPE)+1),leqdim))
			allocate (amat(nbdim))
			allocate (bmat(nbdim))
			allocate (cmat(nbdim))
			allocate (ipc(nbidim))
			allocate (xt(nskpxn(n_st(myPE))+1:nskpxn(n_nd(myPE)+1)))
			allocate (yt(nskpxn(n_st(myPE))+1:nskpxn(n_nd(myPE)+1)))
		else
			allocate (cmamm(nskp2n(n_st(myPE))+1:mcnt2,leqdim))
			allocate (cmamp(nskp2n(n_st(myPE))+1:mcnt2,leqdim))
			allocate (cmapm(nskp2n(n_st(myPE))+1:mcnt2,leqdim))
			allocate (cmapp(nskp2n(n_st(myPE))+1:mcnt2,leqdim))
			allocate (amat(nbdim))
			allocate (bmat(nbdim))
			allocate (cmat(nbdim))
			allocate (ipc(nbidim))
			allocate (xt(nskpxn(n_st(myPE))+1:nbxdim))
			allocate (yt(nskpxn(n_st(myPE))+1:nbxdim))
		end if

		cmamm=0.
		cmamp=0.
		cmapm=0.
		cmapp=0.

		amat=0.
		bmat=0.
		cmat=0.

		ipc=0
		xt=0.
		yt=0.

!  calculate derivative weights to be used in blockj, block0, b2lx and b2lx0

		wt1m(1,1)=0.
		wt2m(1,1)=0.
		wt1m(1,2)=0.
		wt2m(1,2)=0.
		wt10(1,1)=-dc1m(1)-dc1p(1)
		wt20(1,1)=-dc2m(1)-dc2p(1)
		wt10(1,2)=r(2)**2*dc1m(1)/(r(2)**2-r(1)**2)-dc1m(1)-dc1p(1)
		wt20(1,2)=r(2)**2*dc2m(1)/(r(2)**2-r(1)**2)-dc2m(1)-dc2p(1)
		wt1p(1,1)=dc1p(1)
		wt2p(1,1)=dc2p(1)
		wt1p(1,2)=-r(1)**2*dc1m(1)/(r(2)**2-r(1)**2)+dc1p(1)
		wt2p(1,2)=-r(1)**2*dc2m(1)/(r(2)**2-r(1)**2)+dc2p(1)
		wt1m(2:mjm1,1)=dc1m(2:mjm1)
		wt2m(2:mjm1,1)=dc2m(2:mjm1)
		wt10(2:mjm1,1)=-dc1m(2:mjm1)-dc1p(2:mjm1)
		wt20(2:mjm1,1)=-dc2m(2:mjm1)-dc2p(2:mjm1)
		wt1p(2:mjm2,1)=dc1p(2:mjm2)
		wt2p(2:mjm2,1)=dc2p(2:mjm2)
		wt1p(mjm1,1)=0.
		wt2p(mjm1,1)=0.
		wt1m(2:mjm2,2)=dc1m(2:mjm2)
		wt2m(2:mjm2,2)=dc2m(2:mjm2)
		wt1m(mjm1,2)=dc1m(mjm1)-dc1p(mjm1)*(r(mj)-r(mjm1))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2)
		wt2m(mjm1,2)=dc2m(mjm1)-dc2p(mjm1)*(r(mj)-r(mjm1))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2)
		wt10(2:mjm2,2)=-dc1m(2:mjm2)-dc1p(2:mjm2)
		wt20(2:mjm2,2)=-dc2m(2:mjm2)-dc2p(2:mjm2)
		wt10(mjm1,2)=-dc1m(mjm1)-dc1p(mjm1)+dc1p(mjm1)*(r(mj)-r(mjm2))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2)
		wt20(mjm1,2)=-dc2m(mjm1)-dc2p(mjm1)+dc2p(mjm1)*(r(mj)-r(mjm2))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2)
		wt1p(2:mjm2,2)=dc1p(2:mjm2)
		wt2p(2:mjm2,2)=dc2p(2:mjm2)
		wt1p(mjm1,2)=0.
		wt2p(mjm1,2)=0.
	
!  	set source

		allocate (src(0:mj))
		do j=0,mj
			rmrbar = (r(j)-rsrc)/wsrc
			src(j)=0.
			if (abs(rmrbar) < 15.0_IDP) then
				src(j)=asrc*exp(-rmrbar**2)
			else
				src(j)=0.
			end if
		end do

!  read start or continue values into proper arrays

		call dlstar(uzt,phi,-1,sc1,sc2,sc3,0.0_IDP,1.0_IDP)
		if (myPE == numPEsm1) uzt(mj,:)=0.0

		do ind=1,nvar
			select case (ind)
				case (1)
					call trnsfr(psi,1,2)
				case (2)
					call trnsfr(phi,-1,2)
				case (3)
					call trnsfr(pr,1,2)
				case (4)
					call trnsfr(uzt,-1,2)
				case (5)
					call trnsfr(nf,1,2)
				case (6)
					call trnsfr(vprlf,-1,2)
				case (7)
					call trnsfr(vthprlf,-1,2)
				case (8)
					call trnsfr(nalp,1,2)
				case (9)
					call trnsfr(vprlalp,-1,2)
			end select
                	
			l1=0
			do i=nr_start,nr_end
				mnum3=noeqn*mnumn(i)
				lsk=(ind-1)*mnumn(i)+nskpxn(i)
				do l1t=1,mnumn(i)
					l1=l1+1
					lskp=l1t+lsk
					do j=1,mj
						imt1=lskp+mnum3*(j-1)
						xt(imt1)=scp(l1,j)
					end do
				end do
			end do
			do i=nc_start,nc_end
				mnum2=noeqn*mnumn(i)
				mnum3=2*noeqn*mnumn(i)
				lsk=(ind-1)*mnumn(i)+nskpxn(i)
				do l1t=1,mnumn(i)
					l1=l1+1
					lskp=l1t+lsk
					do j=1,mj
						imt1=lskp+mnum3*(j-1)
						xt(imt1)=scp(l1,j)
						xt(imt1+mnum2)=scp(l1+mnumn(i),j)
					end do
				end do
			end do
		end do
	
!  calculate "c" coefficients to be used in blockj, b2lx and b2lx0

		allocate (llno(lmax))
		do l=1,lmax
			do l1=1,lmaxn
				if(lln(l1) == l) exit
			end do
			if (l1 > lmaxn) then
				llno(l)=0
			else
				llno(l)=l1
			end if
		end do

		do it1=-1,1,2

			do it2=-1,1,2

				it3=it1*it2

				do le=1,leqmax
					m1=mmeq(le)*sgnleq(le)
					n1=nneq(le)*sgnleq(le)
					is1=it1*sgnleq(le)

					do i=n_st(myPE),n_nd(myPE)
						do l2t=1,mnumn(i)
							l2=l2t+lnumn(i-1)
							l=lln(l2)
							m2=mm(l)*signl(l)
							n2=nn(l)*signl(l)
							is2=it2*signl(l)

							if (is1 == 0) then
								lp=nskp2n(i)+mnumn(i)*(l2t-1)+l2t
								if (it1 == 1) then
									if (it2 == 1) then
										cmapp(lp,le)=1.0_IDP
									else
										cmapm(lp,le)=1.0_IDP
									end if
								else
									l3=lo(l2)
									if (l3 > 0) then
										l3t=l3-lnumn(i-1)
										if (l3t <= mnumn(i)) then
											lp=nskp2n(i)+mnumn(i)*(l2t-1)+l3t
											if (it2 == 1) then
												cmamp(lp,le)=1.0_IDP
											else
												cmamm(lp,le)=1.0_IDP
											end if
										end if
									end if
								end if

							else if (is2 == 0) then
								m3=it2*mmeq(le)
								n3=it2*nneq(le)
								if (m3 >= mmin .and. m3 <= mmax .and. n3 >= nmin .and. n3 <= nmax) then
									l3p=ll(m3,n3)
									if (l3p > 0) then
										l3=llno(l3p)
										if (l3 > 0) then
											l3t=l3-lnumn(i-1)
											if (l3t <= mnumn(i)) then
												lp=nskp2n(i)+mnumn(i)*(l2t-1)+l3t
												if (it2 == 1) then
													if (it1 == 1) then
														cmapp(lp,le)=1.0_IDP
													else
														cmamp(lp,le)=1.0_IDP
													end if
												else
													if (it1 == 1) then
														cmapm(lp,le)=1.0_IDP
													else
														cmamm(lp,le)=1.0_IDP
													end if
												end if
											end if
										end if
									end if
								end if

							else
								is3=is1*is2

								if (n1 < n2) then
									sgn12=-1
									mm3=m2-m1
									nm3=n2-n1
								else if (n1 > n2) then
									sgn12=1
									mm3=m1-m2
									nm3=n1-n2
								else
									sgn12=sign(1,m1-m2)
									mm3=abs(m1-m2)
									nm3=0
								end if
								m3=is3*it3*mm3
								n3=is3*it3*nm3
								l3p=0
								if (m3 >= mmin .and. m3 <= mmax .and. n3 >= nmin .and. n3 <= nmax) l3p=ll(m3,n3)
								l3=0
								if (l3p > 0) l3=llno(l3p)
								if (l3 > 0) then
									l3t=l3-lnumn(i-1)
									lp=nskp2n(i)+mnumn(i)*(l2t-1)+l3t
									coef=0.5_IDP
									if (is3 == -1) then
										if (is1*sgn12 == 1) coef=-0.5_IDP
										if (mm3 == 0 .and. nm3 == 0) coef=0.0_IDP
									end if
									if (it2 == 1) then
										if (it1 == 1) then
											cmapp(lp,le)=coef
										else
											cmamp(lp,le)=coef
										end if
									else
										if (it1 == 1) then
											cmapm(lp,le)=coef
										else
											cmamm(lp,le)=coef
										end if
									end if
								end if

								mp3=m1+m2
								np3=n1+n2
								m3=is3*it3*mp3
								n3=is3*it3*np3
								l3p=0
								if (m3 >= mmin .and. m3 <= mmax .and. n3 >= nmin .and. n3 <= nmax) l3p=ll(m3,n3)
								l3=0
								if (l3p > 0) l3=llno(l3p)
								if (l3 > 0) then
									l3t=l3-lnumn(i-1)
									lp=nskp2n(i)+mnumn(i)*(l2t-1)+l3t
									coef=0.5_IDP
									if (is1 == -1 .and. is2 == -1) coef=-0.5_IDP
									if (it2 == 1) then
										if (it1 == 1) then
											cmapp(lp,le)=coef
										else
											cmamp(lp,le)=coef
										end if
									else
										if (it1 == 1) then
											cmapm(lp,le)=coef
										else
											cmamm(lp,le)=coef
										end if
									end if
								end if

							end if
						end do
					end do
				end do
			end do
		end do

!   Ion FLR effects auxiliary variable

		if (iflr_on == 1) then
			call dlsq(sc1,psi,1,sc2,sc3,sc4,0.0_IDP,1.0_IDP)
			call trnsfr(sc1,1,2)
			l1=0
			do i=nr_start,nr_end
				mnum3=noeqn*mnumn(i)
				lsk=(iq-1)*mnumn(i)+nskpxn(i)
				do l1t=1,mnumn(i)
					l1=l1+1
					lskp=l1t+lsk
					do j=1,mj
						imt1=lskp+mnum3*(j-1)
						xt(imt1)=scp(l1,j)
					end do
				end do
			end do
			do i=nc_start,nc_end
				mnum2=noeqn*mnumn(i)
				mnum3=2*noeqn*mnumn(i)
				lsk=(iq-1)*mnumn(i)+nskpxn(i)
				do l1t=1,mnumn(i)
					l1=l1+1
					lskp=l1t+lsk
					do j=1,mj
						imt1=lskp+mnum3*(j-1)
						xt(imt1)=scp(l1,j)
						xt(imt1+mnum2)=scp(l1+mnumn(i),j)
					end do
				end do
			end do
		end if

!   EP FLR effects auxiliary variables

		if (epflr_on == 1) then

			call dlsq(sc1,phi,-1,sc2,sc3,sc4,0.0_IDP,1.0_IDP)
			call dbydzt_par(sc2,psi,1,0.0_IDP,1.0_IDP)
			call dbydth_par(sc4,psi,1,0.0_IDP,1.0_IDP,0)
			do l=1,lmax
				sc3(:,l)=r(mj_start:mj_end)*sc4(:,l)
			end do

			allocate (xa(0:mj))
			allocate (scp2(lmaxPE(myPE),0:mj),scp3(lmaxPE(myPE),0:mj))
			call trnsfr(sc2,-1,2)
			scp2=scp
			call trnsfr(sc3,-1,2)
			scp3=scp
			call trnsfr(sc1,-1,2)

			lsk=0
			do i=nr_start,nr_end
				allocate (amatw(mnumn(i),mnumn(i),mjm1),bmatw(mnumn(i),mnumn(i),mjm1),cmatw(mnumn(i),mnumn(i),mjm1))
				allocate (xw(mnumn(i),mjm1))
				allocate (ipcw(mnumn(i),mjm1))

				amatw=0.0_IDP
				bmatw=0.0_IDP
				cmatw=0.0_IDP
				mnum3=noeqn*mnumn(i)

				do l2t=1,mnumn(i)
					l2=l2t+lnumn(i-1)
					m=mm(lln(l2))
					n=nn(lln(l2))
					ibnd=1
					if (m == 0) ibnd=2
					do l3t=1,mnumn(i)
						do lq=1,leqmax
							lp=nskp2n(i)+mnumn(i)*(l2t-1)+l3t
							x=cmapm(lp,lq)
							if (x == 0.0_IDP) cycle
							xa=x
							do j=1,mjm2
								amatw(l3t,l2t,j)=amatw(l3t,l2t,j)+wt20(j,ibnd)*xa(j)*lplrr(j,lq)
								cmatw(l3t,l2t,j)=cmatw(l3t,l2t,j)+wt2m(j,ibnd)*xa(j)*lplrr(j,lq)
								bmatw(l3t,l2t,j)=bmatw(l3t,l2t,j)+wt2p(j,ibnd)*xa(j)*lplrr(j,lq)
							end do
							amatw(l3t,l2t,mjm1)=amatw(l3t,l2t,mjm1)+wt20(mjm1,1)*xa(mjm1)*lplrr(mjm1,lq)
							cmatw(l3t,l2t,mjm1)=cmatw(l3t,l2t,mjm1)+wt2m(mjm1,1)*xa(mjm1)*lplrr(mjm1,lq)
							bmatw(l3t,l2t,mjm1)=bmatw(l3t,l2t,mjm1)+wt2p(mjm1,1)*xa(mjm1)*lplrr(mjm1,lq)
							do j=1,mjm2
								amatw(l3t,l2t,j)=amatw(l3t,l2t,j)+wt10(j,ibnd)*xa(j)*lplr(j,lq)
								cmatw(l3t,l2t,j)=cmatw(l3t,l2t,j)+wt1m(j,ibnd)*xa(j)*lplr(j,lq)
								bmatw(l3t,l2t,j)=bmatw(l3t,l2t,j)+wt1p(j,ibnd)*xa(j)*lplr(j,lq)
							end do
							amatw(l3t,l2t,mjm1)=amatw(l3t,l2t,mjm1)+wt10(mjm1,1)*xa(mjm1)*lplr(mjm1,lq)
							cmatw(l3t,l2t,mjm1)=cmatw(l3t,l2t,mjm1)+wt1m(mjm1,1)*xa(mjm1)*lplr(mjm1,lq)
							bmatw(l3t,l2t,mjm1)=bmatw(l3t,l2t,mjm1)+wt1p(mjm1,1)*xa(mjm1)*lplr(mjm1,lq)
							xa=-x*m*m*rinv*rinv
							do j=1,mjm1
								amatw(l3t,l2t,j)=amatw(l3t,l2t,j)+xa(j)*lpltt(j,lq)
							end do
							xa=-x*m*n*rinv
							do j=1,mjm1
								amatw(l3t,l2t,j)=amatw(l3t,l2t,j)+xa(j)*lpltz(j,lq)
							end do
							xa=-x*n*n
							do j=1,mjm1
								amatw(l3t,l2t,j)=amatw(l3t,l2t,j)+xa(j)*lplzz(j,lq)
							end do
						end do
						do lq=1,leqmax
							lp=nskp2n(i)+mnumn(i)*(l2t-1)+l3t
							x=cmamp(lp,lq)
							if (x == 0.0_IDP) cycle
							xa=x*m*rinv
							do j=1,mjm2
								amatw(l3t,l2t,j)=amatw(l3t,l2t,j)+wt10(j,ibnd)*xa(j)*lplrt(j,lq)
								cmatw(l3t,l2t,j)=cmatw(l3t,l2t,j)+wt1m(j,ibnd)*xa(j)*lplrt(j,lq)
								bmatw(l3t,l2t,j)=bmatw(l3t,l2t,j)+wt1p(j,ibnd)*xa(j)*lplrt(j,lq)
							end do
							amatw(l3t,l2t,mjm1)=amatw(l3t,l2t,mjm1)+wt10(mjm1,1)*xa(mjm1)*lplrt(mjm1,lq)
							cmatw(l3t,l2t,mjm1)=cmatw(l3t,l2t,mjm1)+wt1m(mjm1,1)*xa(mjm1)*lplrt(mjm1,lq)
							bmatw(l3t,l2t,mjm1)=bmatw(l3t,l2t,mjm1)+wt1p(mjm1,1)*xa(mjm1)*lplrt(mjm1,lq)
							do j=1,mjm1
								amatw(l3t,l2t,j)=amatw(l3t,l2t,j)+xa(j)*lplt(j,lq)
							end do
							xa=x*n
							do j=1,mjm2
								amatw(l3t,l2t,j)=amatw(l3t,l2t,j)+wt10(j,ibnd)*xa(j)*lplrz(j,lq)
								cmatw(l3t,l2t,j)=cmatw(l3t,l2t,j)+wt1m(j,ibnd)*xa(j)*lplrz(j,lq)
								bmatw(l3t,l2t,j)=bmatw(l3t,l2t,j)+wt1p(j,ibnd)*xa(j)*lplrz(j,lq)
							end do
							amatw(l3t,l2t,mjm1)=amatw(l3t,l2t,mjm1)+wt10(mjm1,1)*xa(mjm1)*lplrz(mjm1,lq)
							cmatw(l3t,l2t,mjm1)=cmatw(l3t,l2t,mjm1)+wt1m(mjm1,1)*xa(mjm1)*lplrz(mjm1,lq)
							bmatw(l3t,l2t,mjm1)=bmatw(l3t,l2t,mjm1)+wt1p(mjm1,1)*xa(mjm1)*lplrz(mjm1,lq)
							do j=1,mjm1
								amatw(l3t,l2t,j)=amatw(l3t,l2t,j)+xa(j)*lplz(j,lq)
							end do
						end do
					end do
				end do

				if (alpha_on == 1 .and. r_epflralp > 0.0) then
					allocate (amatwalp(mnumn(i),mnumn(i),mjm1),bmatwalp(mnumn(i),mnumn(i),mjm1),cmatwalp(mnumn(i),mnumn(i),mjm1))
					allocate (ipcwalp(mnumn(i),mjm1))
					coef=-r_epflralp*r_epflralp
					amatwalp=coef*amatw
					bmatwalp=coef*bmatw
					cmatwalp=coef*cmatw
					do l=1,mnumn(i)
						do j=1,mjm1
							amatwalp(l,l,j)=1.0+amatwalp(l,l,j)
						end do
					end do

					bmatwalp(:,:,mjm1)=0.0_IDP
					cmatwalp(:,:,1)=0.0_IDP
					call decbt(mnumn(i),mjm1,amatwalp,bmatwalp,cmatwalp,ipcwalp,ier)
					if (ier /= 0) stop 18

					do j=1,mjm1
						xw(:,j)=coef*scp(lsk+1:lsk+mnumn(i),j)
					end do
					call solbt(mnumn(i),mjm1,amatwalp,bmatwalp,cmatwalp,xw,ipcwalp)
					do l=1,mnumn(i)
						lskp=l+(iwa-1)*mnumn(i)+nskpxn(i)
						do j=1,mjm1
							imt1=lskp+mnum3*(j-1)
							xt(imt1)=xw(l,j)
						end do
					end do

					do j=1,mjm1
						xw(:,j)=scp2(lsk+1:lsk+mnumn(i),j)
					end do
					call solbt(mnumn(i),mjm1,amatwalp,bmatwalp,cmatwalp,xw,ipcwalp)
					do l=1,mnumn(i)
						lskp=l+(ix1a-1)*mnumn(i)+nskpxn(i)
						do j=1,mjm1
							imt1=lskp+mnum3*(j-1)
							xt(imt1)=xw(l,j)
						end do
					end do

					do j=1,mjm1
						xw(:,j)=scp3(lsk+1:lsk+mnumn(i),j)
					end do
					call solbt(mnumn(i),mjm1,amatwalp,bmatwalp,cmatwalp,xw,ipcwalp)
					do l=1,mnumn(i)
						lskp=l+(ix2a-1)*mnumn(i)+nskpxn(i)
						do j=1,mjm1
							imt1=lskp+mnum3*(j-1)
							xt(imt1)=xw(l,j)
						end do
					end do
					deallocate (amatwalp,bmatwalp,cmatwalp,ipcwalp)
				end if

				if (r_epflr > 0.0) then
					coef=-r_epflr*r_epflr
					amatw=coef*amatw
					bmatw=coef*bmatw
					cmatw=coef*cmatw
					do l=1,mnumn(i)
						do j=1,mjm1
							amatw(l,l,j)=1.0+amatw(l,l,j)
						end do
					end do

					bmatw(:,:,mjm1)=0.0_IDP
					cmatw(:,:,1)=0.0_IDP

					call decbt(mnumn(i),mjm1,amatw,bmatw,cmatw,ipcw,ier)
					if (ier /= 0) stop 17

					do j=1,mjm1
						xw(:,j)=coef*scp(lsk+1:lsk+mnumn(i),j)
					end do
					call solbt(mnumn(i),mjm1,amatw,bmatw,cmatw,xw,ipcw)
					do l=1,mnumn(i)
						lskp=l+(iw-1)*mnumn(i)+nskpxn(i)
						do j=1,mjm1
							imt1=lskp+mnum3*(j-1)
							xt(imt1)=xw(l,j)
						end do
					end do

					do j=1,mjm1
						xw(:,j)=scp2(lsk+1:lsk+mnumn(i),j)
					end do
					call solbt(mnumn(i),mjm1,amatw,bmatw,cmatw,xw,ipcw)
					do l=1,mnumn(i)
						lskp=l+(ix1-1)*mnumn(i)+nskpxn(i)
						do j=1,mjm1
							imt1=lskp+mnum3*(j-1)
							xt(imt1)=xw(l,j)
						end do
					end do

					do j=1,mjm1
						xw(:,j)=scp3(lsk+1:lsk+mnumn(i),j)
					end do
					call solbt(mnumn(i),mjm1,amatw,bmatw,cmatw,xw,ipcw)
					do l=1,mnumn(i)
						lskp=l+(ix2-1)*mnumn(i)+nskpxn(i)
						do j=1,mjm1
							imt1=lskp+mnum3*(j-1)
							xt(imt1)=xw(l,j)
						end do
					end do
				end if
				lsk=lsk+mnumn(i)
				deallocate (amatw,bmatw,cmatw,xw,ipcw)
			end do

			do i=nc_start,nc_end
				allocate (amatw(2*mnumn(i),mnumn(i),mjm1),bmatw(2*mnumn(i),mnumn(i),mjm1),cmatw(2*mnumn(i),mnumn(i),mjm1))
				allocate (xw(2*mnumn(i),mjm1))
				allocate (ipcw(mnumn(i),mjm1))

				amatw=0.0_IDP
				bmatw=0.0_IDP
				cmatw=0.0_IDP
				mnum2=noeqn*mnumn(i)
				mnum3=2*noeqn*mnumn(i)

				do l2t=1,mnumn(i)
					l2=l2t+lnumn(i-1)
					m=mm(lln(l2))
					n=nn(lln(l2))
					ibnd=1
					if (m == 0) ibnd=2
					do l3t=1,mnumn(i)
						do lq=1,leqmax
							lp=nskp2n(i)+mnumn(i)*(l2t-1)+l3t
							x=cmapm(lp,lq)
							if (x == 0.0_IDP) cycle
							xa=x
							do j=1,mjm2
								amatw(l3t,l2t,j)=amatw(l3t,l2t,j)+wt20(j,ibnd)*xa(j)*lplrr(j,lq)
								cmatw(l3t,l2t,j)=cmatw(l3t,l2t,j)+wt2m(j,ibnd)*xa(j)*lplrr(j,lq)
								bmatw(l3t,l2t,j)=bmatw(l3t,l2t,j)+wt2p(j,ibnd)*xa(j)*lplrr(j,lq)
							end do
							amatw(l3t,l2t,mjm1)=amatw(l3t,l2t,mjm1)+wt20(mjm1,1)*xa(mjm1)*lplrr(mjm1,lq)
							cmatw(l3t,l2t,mjm1)=cmatw(l3t,l2t,mjm1)+wt2m(mjm1,1)*xa(mjm1)*lplrr(mjm1,lq)
							bmatw(l3t,l2t,mjm1)=bmatw(l3t,l2t,mjm1)+wt2p(mjm1,1)*xa(mjm1)*lplrr(mjm1,lq)
							do j=1,mjm2
								amatw(l3t,l2t,j)=amatw(l3t,l2t,j)+wt10(j,ibnd)*xa(j)*lplr(j,lq)
								cmatw(l3t,l2t,j)=cmatw(l3t,l2t,j)+wt1m(j,ibnd)*xa(j)*lplr(j,lq)
								bmatw(l3t,l2t,j)=bmatw(l3t,l2t,j)+wt1p(j,ibnd)*xa(j)*lplr(j,lq)
							end do
							amatw(l3t,l2t,mjm1)=amatw(l3t,l2t,mjm1)+wt10(mjm1,1)*xa(mjm1)*lplr(mjm1,lq)
							cmatw(l3t,l2t,mjm1)=cmatw(l3t,l2t,mjm1)+wt1m(mjm1,1)*xa(mjm1)*lplr(mjm1,lq)
							bmatw(l3t,l2t,mjm1)=bmatw(l3t,l2t,mjm1)+wt1p(mjm1,1)*xa(mjm1)*lplr(mjm1,lq)
							xa=-x*m*m*rinv*rinv
							do j=1,mjm1
								amatw(l3t,l2t,j)=amatw(l3t,l2t,j)+xa(j)*lpltt(j,lq)
							end do
							xa=-x*m*n*rinv
							do j=1,mjm1
								amatw(l3t,l2t,j)=amatw(l3t,l2t,j)+xa(j)*lpltz(j,lq)
							end do
							xa=-x*n*n
							do j=1,mjm1
								amatw(l3t,l2t,j)=amatw(l3t,l2t,j)+xa(j)*lplzz(j,lq)
							end do
						end do
						do lq=1,leqmax
							lp=nskp2n(i)+mnumn(i)*(l2t-1)+l3t
							x=cmamp(lp,lq)
							if (x == 0.0_IDP) cycle
							xa=x*m*rinv
							do j=1,mjm2
								amatw(l3t,l2t,j)=amatw(l3t,l2t,j)+wt10(j,ibnd)*xa(j)*lplrt(j,lq)
								cmatw(l3t,l2t,j)=cmatw(l3t,l2t,j)+wt1m(j,ibnd)*xa(j)*lplrt(j,lq)
								bmatw(l3t,l2t,j)=bmatw(l3t,l2t,j)+wt1p(j,ibnd)*xa(j)*lplrt(j,lq)
							end do
							amatw(l3t,l2t,mjm1)=amatw(l3t,l2t,mjm1)+wt10(mjm1,1)*xa(mjm1)*lplrt(mjm1,lq)
							cmatw(l3t,l2t,mjm1)=cmatw(l3t,l2t,mjm1)+wt1m(mjm1,1)*xa(mjm1)*lplrt(mjm1,lq)
							bmatw(l3t,l2t,mjm1)=bmatw(l3t,l2t,mjm1)+wt1p(mjm1,1)*xa(mjm1)*lplrt(mjm1,lq)
							do j=1,mjm1
								amatw(l3t,l2t,j)=amatw(l3t,l2t,j)+xa(j)*lplt(j,lq)
							end do
							xa=x*n
							do j=1,mjm2
								amatw(l3t,l2t,j)=amatw(l3t,l2t,j)+wt10(j,ibnd)*xa(j)*lplrz(j,lq)
								cmatw(l3t,l2t,j)=cmatw(l3t,l2t,j)+wt1m(j,ibnd)*xa(j)*lplrz(j,lq)
								bmatw(l3t,l2t,j)=bmatw(l3t,l2t,j)+wt1p(j,ibnd)*xa(j)*lplrz(j,lq)
							end do
							amatw(l3t,l2t,mjm1)=amatw(l3t,l2t,mjm1)+wt10(mjm1,1)*xa(mjm1)*lplrz(mjm1,lq)
							cmatw(l3t,l2t,mjm1)=cmatw(l3t,l2t,mjm1)+wt1m(mjm1,1)*xa(mjm1)*lplrz(mjm1,lq)
							bmatw(l3t,l2t,mjm1)=bmatw(l3t,l2t,mjm1)+wt1p(mjm1,1)*xa(mjm1)*lplrz(mjm1,lq)
							do j=1,mjm1
								amatw(l3t,l2t,j)=amatw(l3t,l2t,j)+xa(j)*lplz(j,lq)
							end do
						end do
					end do
				end do

				if (alpha_on == 1 .and. r_epflralp > 0.0) then
					allocate (amatwalp(2*mnumn(i),mnumn(i),mjm1),bmatwalp(2*mnumn(i),mnumn(i),mjm1),cmatwalp(2*mnumn(i),mnumn(i),mjm1))
					allocate (ipcwalp(mnumn(i),mjm1))
					coef=-r_epflralp*r_epflralp
					amatwalp=coef*amatw
					bmatwalp=coef*bmatw
					cmatwalp=coef*cmatw
					do l=1,mnumn(i)
						do j=1,mjm1
							amatwalp(l,l,j)=1.0+amatwalp(l,l,j)
						end do
					end do

					bmatwalp(:,:,mjm1)=0.0_IDP
					cmatwalp(:,:,1)=0.0_IDP

					call cdecbt(mnumn(i),mjm1,amatwalp,bmatwalp,cmatwalp,ipcwalp,ier)
					if (ier /= 0) stop 18

					do j=1,mjm1
						xw(:,j)=coef*scp(lsk+1:lsk+2*mnumn(i),j)
					end do
					call csolbt(mnumn(i),mjm1,amatwalp,bmatwalp,cmatwalp,xw,ipcwalp)
					do l=1,mnumn(i)
						lskp=l+(iwa-1)*mnumn(i)+nskpxn(i)
						do j=1,mjm1
							imt1=lskp+mnum3*(j-1)
							xt(imt1)=xw(l,j)
							xt(imt1+mnum2)=xw(l+mnumn(i),j)
						end do
					end do

					do j=1,mjm1
						xw(:,j)=scp2(lsk+1:lsk+2*mnumn(i),j)
					end do
					call csolbt(mnumn(i),mjm1,amatwalp,bmatwalp,cmatwalp,xw,ipcwalp)
					do l=1,mnumn(i)
						lskp=l+(ix1a-1)*mnumn(i)+nskpxn(i)
						do j=1,mjm1
							imt1=lskp+mnum3*(j-1)
							xt(imt1)=xw(l,j)
							xt(imt1+mnum2)=xw(l+mnumn(i),j)
						end do
					end do

					do j=1,mjm1
						xw(:,j)=scp3(lsk+1:lsk+2*mnumn(i),j)
					end do
					call csolbt(mnumn(i),mjm1,amatwalp,bmatwalp,cmatwalp,xw,ipcwalp)
					do l=1,mnumn(i)
						lskp=l+(ix2a-1)*mnumn(i)+nskpxn(i)
						do j=1,mjm1
							imt1=lskp+mnum3*(j-1)
							xt(imt1)=xw(l,j)
							xt(imt1+mnum2)=xw(l+mnumn(i),j)
						end do
					end do
					deallocate (amatwalp,bmatwalp,cmatwalp,ipcwalp)
				end if

				if (r_epflr > 0.0) then
					coef=-r_epflr*r_epflr
					amatw=coef*amatw
					bmatw=coef*bmatw
					cmatw=coef*cmatw
					do l=1,mnumn(i)
						do j=1,mjm1
							amatw(l,l,j)=1.0+amatw(l,l,j)
						end do
					end do

					bmatw(:,:,mjm1)=0.0_IDP
					cmatw(:,:,1)=0.0_IDP

					call cdecbt(mnumn(i),mjm1,amatw,bmatw,cmatw,ipcw,ier)
					if (ier /= 0) stop 17

					do j=1,mjm1
						xw(:,j)=coef*scp(lsk+1:lsk+2*mnumn(i),j)
					end do
					call csolbt(mnumn(i),mjm1,amatw,bmatw,cmatw,xw,ipcw)
					do l=1,mnumn(i)
						lskp=l+(iw-1)*mnumn(i)+nskpxn(i)
						do j=1,mjm1
							imt1=lskp+mnum3*(j-1)
							xt(imt1)=xw(l,j)
							xt(imt1+mnum2)=xw(l+mnumn(i),j)
						end do
					end do

					do j=1,mjm1
						xw(:,j)=scp2(lsk+1:lsk+2*mnumn(i),j)
					end do
					call csolbt(mnumn(i),mjm1,amatw,bmatw,cmatw,xw,ipcw)
					do l=1,mnumn(i)
						lskp=l+(ix1-1)*mnumn(i)+nskpxn(i)
						do j=1,mjm1
							imt1=lskp+mnum3*(j-1)
							xt(imt1)=xw(l,j)
							xt(imt1+mnum2)=xw(l+mnumn(i),j)
						end do
					end do

					do j=1,mjm1
						xw(:,j)=scp3(lsk+1:lsk+2*mnumn(i),j)
					end do
					call csolbt(mnumn(i),mjm1,amatw,bmatw,cmatw,xw,ipcw)
					do l=1,mnumn(i)
						lskp=l+(ix2-1)*mnumn(i)+nskpxn(i)
						do j=1,mjm1
							imt1=lskp+mnum3*(j-1)
							xt(imt1)=xw(l,j)
							xt(imt1+mnum2)=xw(l+mnumn(i),j)
						end do
					end do
				end if
				lsk=lsk+2*mnumn(i)
				deallocate (amatw,bmatw,cmatw,xw,ipcw)
			end do
			deallocate (xa,scp2,scp3)

		end if

!  put in r.h.s. of equations

!  psi equation

		sd1=1.0_IDP
		call om0(sd1,1,2,0,0,0,1.0_IDP)
		call clgam(sceq2,1,1,1,0.0_IDP,-1.0_IDP)
		do l=1,leqmax
			sceq1(:,l)=eta*feq*sceq2(:,l)
		end do
		call blockj(sceq1,-1,1,1,1,0,0,oneos)
		call clgam(sceq2,1,2,1,0.0_IDP,1.0_IDP)
		do l=1,leqmax
			sceq1(:,l)=eta*feq*sceq2(:,l)
		end do
		call blockj(sceq1,1,1,1,0,1,0,oneos)
		do l=1,leqmax
			sceq1(:,l)=eta*feq*grroj(:,l)
		end do
		call blockj(sceq1,1,1,1,2,0,0,oneos)
		do l=1,leqmax
			sceq1(:,l)=-2.*eta*feq*grtoj(:,l)
		end do
		call blockj(sceq1,-1,1,1,1,1,0,oneos)
		do l=1,leqmax
			sceq1(:,l)=eta*feq*gttoj(:,l)
		end do
		call blockj(sceq1,1,1,1,0,2,0,oneos)
		do l=1,leqmax
			sceq2(:,l)=rinv*sceq1(:,l)
		end do
		call blockj(sceq2,1,1,1,0,1,0,oneos)

!   Ion FLR effects

		if (iflr_on == 1) then

			if (ext_prof == 1) then
				do l=1,leqmax
					sceq5(:,l)=1.2533*iflr*iflr*vAlfven*vAlfven*bmod(:,l)/(vtherm_elecP*(feq-qqinv*cureq))
				end do
			else
				do l=1,leqmax
					sceq5(:,l)=1.2533*iflr*iflr*bmod(:,l)/(denseq*vtherm_elc*(feq-qqinv*cureq))
				end do
			end if
			call blockj_landau_grad_parallel(sceq5,1,1,iq,0,0,0,1.0_IDP)

!   Ion FLR effects auxiliary equation

			call block_dlsq(iq,1,1.0_IDP)
			sd1=1.0_IDP
			call block0(sd1,iq,iq,0,0,0,-1.0_IDP)

		end if

!   Two fluid terms

		if (twofl_on == 1) then
			sd1=beteom/denseq
			call om0(sd1,1,3,0,0,0,-1.0_IDP)
		end if

!  u-zeta equation

		call dbydreq(sceq1,sqg,0.0_IDP,-1.0_IDP,0)
		call blockj(sceq1,1,2,3,1,0,0,betfc)
		call dbydtheq(sceq2,sqg,1,0.0_IDP,1.0_IDP,0)
		call blockj(sceq2,-1,2,3,0,1,0,betfc)
		write(*,'(/"bet0 = ",e12.4,3x,"betf0_f = ",e12.4)') bet0, bet0_f
		write(*,'("eps = <a>/<R> = ",e12.4,3x,"S = tau_R/tau_A = ",e12.4)') eps, s
		write(*,'("bet0/eps**2 = ",e12.4,3x,"bet0_f/eps**2 = ",e12.4)') betfc, betfc_f

!  fast ion coupling

		call blockj(sceq1,1,2,5,1,0,0,betfc_f)
		call blockj(sceq2,-1,2,5,0,1,0,betfc_f)

		if (alpha_on == 1) then
			call blockj(sceq1,1,2,8,1,0,0,betfc_alp)
			call blockj(sceq2,-1,2,8,0,1,0,betfc_alp)
		end if

!  Shared equilibrium toroidal flow velocity for u-zeta equation

		call block0(vzt_eq,2,4,0,0,1,-1.0_IDP)

		call clgam(sceq2,1,1,1,0.0_IDP,-1.0_IDP)
		call om(sceq2,-1,2,1,1,0,0,1.0_IDP)
		call grpareq(sceq1,sceq2,-1,0.0_IDP,1.0_IDP)
		call blockj(sceq1,1,2,1,1,0,0,1.0_IDP)
		call clgam(sceq2,1,2,1,0.0_IDP,1.0_IDP)
		call om(sceq2,1,2,1,0,1,0,1.0_IDP)
		call grpareq(sceq1,sceq2,1,0.0_IDP,1.0_IDP)
		call blockj(sceq1,-1,2,1,0,1,0,1.0_IDP)
		call grpareq(sceq1,grtoj,-1,0.0_IDP,-2.0_IDP)
		call blockj(sceq1,1,2,1,1,1,0,1.0_IDP)
		call grpareq(sceq1,gttoj,1,0.0_IDP,1.0_IDP)
		call blockj(sceq1,-1,2,1,0,2,0,1.0_IDP)
		do l=1,leqmax
			sceq2(:,l)=rinv*sceq1(:,l)
		end do
		call blockj(sceq2,-1,2,1,0,1,0,1.0_IDP)
		call grpareq(sceq1,grroj,1,0.0_IDP,1.0_IDP)
		call blockj(sceq1,-1,2,1,2,0,0,1.0_IDP)
		do l=1,leqmax
			sceq1(:,l)=-grtoj(:,l)
		end do
		call om(sceq1,-1,2,1,1,1,0,2.0_IDP*1.0_IDP)
		do l=1,leqmax
			sceq1(:,l)=grroj(:,l)
		end do
		call om(sceq1,1,2,1,2,0,0,1.0_IDP)
		do l=1,leqmax
			sceq1(:,l)=gttoj(:,l)
		end do
		call om(sceq1,1,2,1,0,2,0,1.0_IDP)
		do l=1,leqmax
			sceq2(:,l)=rinv*sceq1(:,l)
		end do
		call om(sceq2,1,2,1,0,1,0,1.0_IDP)
		call dbydr0(sd1,cureq,0.0_IDP,1.0_IDP,0)
		sd2=rinv*sd1
		call d2bydr20(sd2,cureq,1.0_IDP,-1.0_IDP,0)
		sd1=rinv*sd2/epsq
		call block0(sd1,2,1,1,0,0,1.0_IDP)
		call dbydtheq(sceq2,bst,-1,0.0_IDP,-1.0_IDP,0)
		do l=1,leqmax
			sceq1(:,l)=r*sceq2(:,l)/epsq
		end do
		call dbydreq(sceq2,sceq1,0.0_IDP,-1.0_IDP,0)
		call blockj(sceq2,1,2,1,1,0,0,1.0_IDP)
		call dbydtheq(sceq2,sceq1,1,0.0_IDP,1.0_IDP,0)
		call blockj(sceq2,-1,2,1,0,1,0,1.0_IDP)

!   Ion FLR effects

		if (iflr_on == 1) then
			coef=omegar*iflr*iflr
			call block_dlsq(2,4,coef)
		end if

!   Electron-ion Landau damping

		if (ieldamp_on == 1) then   
			coef=-(1-dpres)*bet0/(2.*epsq*omegar)
			call blockj(eildrr,1,2,2,0,2,0,coef)
			call blockj(eildrt,-1,2,2,1,1,0,coef)
			call blockj(eildrz,-1,2,2,0,1,1,coef)
			call blockj(eildtt,1,2,2,2,0,0,coef)
			call blockj(eildtz,1,2,2,1,0,1,coef)
			call blockj(eildzz,1,2,2,0,0,2,coef)
			call blockj(eildr,1,2,2,0,1,0,coef)
			call blockj(eildt,-1,2,2,1,0,0,coef)
			call blockj(eildz,-1,2,2,0,0,1,coef)
		end if

!   Two fluid terms

		if (twofl_on == 1) then
			call dbydr0(sd1,preq,0.0_IDP,1.0_IDP,0)
			sd2=feq-qqinv*cureq
			do l=1,leqmax
				sceq1(:,l)=feq*sd1*jbgrr(:,l)/sd2
				sceq2(:,l)=feq*sd1*jbgrt(:,l)/sd2
				sceq3(:,l)=feq*sd1*jbgtt(:,l)/sd2
			end do
			call blockj(sceq1,1,2,2,3,0,0,-betiom)
			call blockj(sceq2,-1,2,2,2,1,0,2*betiom)
			call blockj(sceq3,1,2,2,1,2,0,-betiom)
			call dbydtheq(sceq4,sceq1,1,0.0_IDP,1.0_IDP,3)
			call dbydreq(sceq4,sceq2,1.0_IDP,-1.0_IDP,3)
			do l=1,leqmax
				sceq4(:,l)=sceq4(:,l)+rinv*sceq2(:,l)
			end do
			call blockj(sceq4,-1,2,2,2,0,0,-betiom)
			call dbydtheq(sceq4,sceq2,-1,0.0_IDP,1.0_IDP,3)
			call dbydreq(sceq4,sceq3,1.0_IDP,-1.0_IDP,3)
			call blockj(sceq4,1,2,2,1,1,0,betiom)
			do l=1,leqmax
				sceq1(:,l)=rinv*cureq*sd1*jbgrr(:,l)/sd2
				sceq2(:,l)=rinv*cureq*sd1*jbgrt(:,l)/sd2
				sceq3(:,l)=rinv*cureq*sd1*jbgtt(:,l)/sd2
			end do
			call blockj(sceq1,1,2,2,2,0,1,betiom)
			call blockj(sceq2,-1,2,2,1,1,1,-2*betiom)
			call blockj(sceq3,1,2,2,0,2,1,betiom)
			call dbydtheq(sceq4,sceq1,1,0.0_IDP,1.0_IDP,4)
			call dbydreq(sceq4,sceq2,1.0_IDP,-1.0_IDP,4)
			call blockj(sceq4,-1,2,2,1,0,1,betiom)
			call dbydtheq(sceq4,sceq2,-1,0.0_IDP,1.0_IDP,4)
			call dbydreq(sceq4,sceq3,1.0_IDP,-1.0_IDP,4)
			do l=1,leqmax
				sceq4(:,l)=sceq4(:,l)-rinv*sceq3(:,l)
			end do
			call blockj(sceq4,1,2,2,0,1,1,-betiom)
			do l=1,leqmax
				sceq1(:,l)=rinv*cureq*sd1*dgrrz(:,l)
				sceq2(:,l)=rinv*cureq*sd1*dgrtz(:,l)
				sceq3(:,l)=rinv*cureq*sd1*dgttz(:,l)
			end do
			call blockj(sceq1,-1,2,2,2,0,0,0.5*betiom)
			call blockj(sceq2,1,2,2,1,1,0,-betiom)
			call blockj(sceq3,-1,2,2,0,2,0,0.5*betiom)
			sd3=rinv*feq*sd1
			sd3(0)=2*feq(0)*(preq(1)-preq(0))/(r(1)*r(1))
			do l=1,leqmax
				sceq1(:,l)=sd3*dgrrt(:,l)-rinv*cureq*sd1*dgrrz(:,l)
				sceq2(:,l)=feq*sd1*dgttr(:,l)+2*sd3*jbgtt(:,l)/sd2-rinv*cureq*sd1*dgrtz(:,l)
			end do
			call dbydtheq(sceq3,sceq1,-1,0.0_IDP,1.0_IDP,2)
			call dbydreq(sceq3,sceq2,1.0_IDP,-1.0_IDP,2)
			call dbydr0(sd4,qqinv,0.0_IDP,1.0_IDP,0)
			call dbydr0(sd5,cureq,0.0_IDP,1.0_IDP,0)
			call dbydreq(sceq2,jsq,0.0_IDP,0.5_IDP,0)
			do l=1,leqmax
				sceq1(:,l)=cureq*sd1*(qqinv*(rinv*dgrtt(:,l)-dgttr(:,l))-(sd4+2*rinv*qqinv)*jbgtt(:,l)/sd2+ &
						      rinv*(dbsjtbj(:,l)-rinv*(cureq*sceq2(:,l)+sd5*jsq(:,l))/sd2)/(eps*eps))
			end do
			call dbydreq(sceq3,sceq1,1.0_IDP,-1.0_IDP,4)
			call blockj(sceq3,1,2,2,1,0,0,-0.5*betiom)
			call dbydtheq(sceq3,sceq1,1,0.0_IDP,1.0_IDP,4)
			do l=1,leqmax
				sceq1(:,l)=sd1*(feq*dgttt(:,l)-cureq*dgttz(:,l))
				sceq2(:,l)=2*sd3*(dgrtt(:,l)-jbgtt(:,l)/sd2)-sd1*(feq*dgttr(:,l)+rinv*cureq*dgrtz(:,l))
			end do
			call dbydreq(sceq5,sceq1,0.0_IDP,1.0_IDP,3)
			do l=1,leqmax
				sceq3(:,l)=sceq3(:,l)+rinv*sceq5(:,l)
			end do
			call dbydtheq(sceq3,sceq2,1,1.0_IDP,-1.0_IDP,2)
			call blockj(sceq3,-1,2,2,0,1,0,-0.5*betiom)
		end if

! diffusion term added
		sd1=1.0_IDP
!		call block0_dlsq(sd1,2,4,stdifu)
                call block0_dlsqr(sd1,2,4,stdifur,sd7)

! p equation

		call dbydr0(sd1,preq,0.0_IDP,1.0_IDP,0)
		call block0(sd1,3,2,1,0,0,1.0_IDP)

!		do l=1,leqmax
!			sceq1(:,l)=preq*djroj(:,l)
!		end do
!		call blockj(sceq1,1,3,2,1,0,0,gamma)
!		do l=1,leqmax
!			sceq1(:,l)=-preq*djtoj(:,l)
!		end do
!		call blockj(sceq1,-1,3,2,0,1,0,gamma)

		sd2=feq/(feq-qqinv*cureq)
		sd3=cureq/(feq-qqinv*cureq)
		do l=1,leqmax
			sceq1(:,l)=-preq*(sd2*djtoj(:,l)-rinv*sd3*djzoj(:,l))
		end do
		call blockj(sceq1,-1,3,2,0,1,0,gamma)
		do l=1,leqmax
			sceq1(:,l)=-preq*(r*dbsjzoj(:,l)-sd2*djroj(:,l))
		end do
		call blockj(sceq1,1,3,2,1,0,0,gamma)
		do l=1,leqmax
			sceq1(:,l)=-preq*(rinv*sd3*djroj(:,l)-dbsjtoj(:,l))
		end do
		call blockj(sceq1,1,3,2,0,0,1,gamma)
		call dbydr0(sd4,sd2,0.0_IDP,1.0_IDP,0)
		sd2=preq*sd4
		call block0(sd2,3,2,1,0,0,gamma)
		call dbydr0(sd4,sd3,0.0_IDP,1.0_IDP,0)
		sd2=-rinv*preq*sd4
		call block0(sd2,3,2,0,0,1,gamma)

!  parallel thermal velocity term		

		do l=1,leqmax
			sceq1(:,l)=preq*bmod(:,l)/(feq-qqinv*cureq)
		end do
		call om(sceq1,1,3,7,0,0,0,-gamma)
		call grpareq(sceq2,bmod,1,0.0_IDP,1.0_IDP)
		do l=1,leqmax
			sceq1(:,l)=preq*sceq2(:,l)/(feq-qqinv*cureq)
		end do
		call blockj(sceq1,-1,3,7,0,0,0,gamma)

!  Shared equilibrium toroidal flow velocity for pressure equation

		call block0(vzt_eq,3,3,0,0,1,-1.0_IDP)

!   Two fluid terms

		if (twofl_on == 1) then
			sd2=gamma*betiom*preq/denseq
			sd3=sd2/(feq-qqinv*cureq)
 	 		call dbydr0(sd4,feq,0.0_IDP,1.0_IDP,0)
			sd4=sd3*sd4
			call block0(sd4,3,3,1,0,0,1.0_IDP)
 	 		call dbydr0(sd4,cureq,0.0_IDP,1.0_IDP,0)
			sd4=rinv*sd3*sd4
			call block0(sd4,3,3,0,0,1,-1.0_IDP)
			call dbydzteq(sceq1,bst,-1,0.0_IDP,1.0_IDP)
			call dbydtheq(sceq2,bst,-1,0.0_IDP,1.0_IDP,0)
			do l=1,leqmax
				sceq1(:,l)=r*sd3*sceq1(:,l)
				sceq2(:,l)=r*sd3*sceq2(:,l)
			end do
			call blockj(sceq1,1,3,3,1,0,0,-1.0_IDP)
			call blockj(sceq2,1,3,3,0,0,1,1.0_IDP)
			do l=1,leqmax
				sceq1(:,l)=sd2*omdr(:,l)
				sceq2(:,l)=sd2*omdt(:,l)
				sceq3(:,l)=sd2*omdz(:,l)
			end do
			call blockj(sceq1,-1,3,3,0,1,0,-2.0_IDP)
			call blockj(sceq2,1,3,3,1,0,0,-2.0_IDP)
			call blockj(sceq3,1,3,3,0,0,1,-2.0_IDP)
			do l=1,leqmax
				sceq1(:,l)=sd1*sd3*grtoj(:,l)
				sceq2(:,l)=sd1*sd3*gttoj(:,l)
			end do
			call om(sceq1,-1,3,1,1,0,0,-epsq)
			call om(sceq2,1,3,1,0,1,0,epsq)
			sd4=rinv*sd1*sd3*cureq
			call block0(sd4,3,1,1,1,0,-1.0_IDP)
			do l=1,leqmax
				sceq1(:,l)=r*sd1*sd3*bst(:,l)
			end do
			call blockj(sceq1,-1,3,1,2,0,0,1.0_IDP)
			do l=1,leqmax
				sceq1(:,l)=sd1*sd2*dgrtp(:,l)
				sceq2(:,l)=sd1*sd2*dgttp(:,l)
			end do
			call blockj(sceq1,1,3,1,1,0,0,-epsq)
			call blockj(sceq2,-1,3,1,0,1,0,epsq)
			do l=1,leqmax
				sceq1(:,l)=rinv*sd1*sd3*cureq*djtoj(:,l)
				sceq2(:,l)=sd1*sd2*dbsjtoj(:,l)
			end do
			call blockj(sceq1,-1,3,1,0,1,0,-1.0_IDP)
			call blockj(sceq2,1,3,1,1,0,0,1.0_IDP)
		end if

! diffusion term added
		sd1=1.0_IDP
!		call block0_dlsq(sd1,3,3,stdifp)
		call block0_dlsqr(sd1,3,3,stdifpr,sd7)

!  u-zeta expression

		do l=1,leqmax
			sceq2(:,l)=denseq*jbgrt(:,l)
		end do
		call blockj(sceq2,-1,4,2,1,1,0,-2.0_IDP)
		do l=1,leqmax
			sceq2(:,l)=denseq*jbgrr(:,l)
		end do
		call blockj(sceq2,1,4,2,2,0,0,1.0_IDP)
		do l=1,leqmax
			sceq2(:,l)=denseq*jbgtt(:,l)
		end do
		call blockj(sceq2,1,4,2,0,2,0,1.0_IDP)
		do l=1,leqmax
			sceq1(:,l)=rinv*sceq2(:,l)
		end do
		call blockj(sceq1,1,4,2,0,1,0,1.0_IDP)
		call clgam(sceq2,1,1,2,0.0_IDP,-1.0_IDP)
		do l=1,leqmax
			sceq1(:,l)=denseq*sceq2(:,l)-denseqr*jbgrt(:,l)
		end do
		call blockj(sceq1,-1,4,2,1,0,0,1.0_IDP)
		call clgam(sceq2,1,2,2,0.0_IDP,1.0_IDP)
		do l=1,leqmax
			sceq1(:,l)=denseq*sceq2(:,l)+denseqr*jbgtt(:,l)
		end do
		call blockj(sceq1,1,4,2,0,1,0,1.0_IDP)
		sd1=-1.0_IDP
		call block0(sd1,4,4,0,0,0,1.0_IDP)
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	NBI particle effects	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		

!  Load 1st Omega-d terms in fast ion density equation:
		write(*,'(/"vfova2 = ",e12.4)') vfova2(1)

        if(Trapped_on .eq. 0) then 
		do l=1,leqmax
			sceq1(:,l)=vfova2*omdr(:,l)/(epsq*omcyd)
			sceq2(:,l)=vfova2*omdt(:,l)/(epsq*omcyd)
			sceq3(:,l)=vfova2*omdz(:,l)/(epsq*omcyd)
		end do
	else
		do l=1,leqmax
			sceq1(:,l)=vfova*omdrprp(:,l)/(epsq*omcyd)
			sceq2(:,l)=vfova*omdtprp(:,l)/(epsq*omcyd)
			sceq3(:,l)=vfova*omdzprp(:,l)/(epsq*omcyd)
		end do
	end if		
		call blockj(sceq1,-1,5,5,0,1,0,-1.0_IDP)
		call blockj(sceq2,1,5,5,1,0,0,-1.0_IDP)
		call blockj(sceq3,1,5,5,0,0,1,-1.0_IDP)

!  Shared equilibrium toroidal flow velocity for fast ion density equation

		call block0(vzt_eq,5,5,0,0,1,-1.0_IDP)

!  Load 1st Omega-d terms in fast ion parallel velocity equation:

		call blockj(sceq1,-1,6,6,0,1,0,-1.0_IDP)
		call blockj(sceq2,1,6,6,1,0,0,-1.0_IDP)
		call blockj(sceq3,1,6,6,0,0,1,-1.0_IDP)

!  Shared equilibrium toroidal flow velocity for fast ion parallel velocity equation

		call block0(vzt_eq,6,6,0,0,1,-1.0_IDP)

!  Load 2nd Omega-d terms in fast ion density equation:
		write(*,'("nfeq = ",e12.4)') nfeq(1)

        if(Trapped_on .eq. 0) then 
		do l=1,leqmax
		       sceq1(:,l)=nfeq(:)*omdr(:,l)
		       sceq2(:,l)=nfeq(:)*omdt(:,l)
		       sceq3(:,l)=nfeq(:)*omdz(:,l)
		end do
	else
		do l=1,leqmax
		       sceq1(:,l)=nfeq(:)*omdrprp(:,l)/vfova
		       sceq2(:,l)=nfeq(:)*omdtprp(:,l)/vfova
		       sceq3(:,l)=nfeq(:)*omdzprp(:,l)/vfova
		end do
	end if	

		call blockj(sceq1,-1,5,2,0,1,0,-1.0_IDP)
		call blockj(sceq2,1,5,2,1,0,0,-1.0_IDP)
		call blockj(sceq3,1,5,2,0,0,1,-1.0_IDP)

! diffusion term added
		sd1=1.0_IDP
!		call block0_dlsq(sd1,5,5,stdifnf)
		call block0_dlsqr(sd1,5,5,stdifnfr,sd7)

!  Load remaining terms in fast ion density equation

!    Parallel gradient term

		do l=1,leqmax
			sceq1(:,l)=nfeq*bmod(:,l)/(feq-qqinv*cureq)
		end do
		
		call om(sceq1,1,5,6,0,0,0,-1.0_IDP)

!    Omega* term

		sd1=dnfeqdr*rinv*cureq/(feq-qqinv*cureq)
		sd2=dnfeqdr*feq/(feq-qqinv*cureq)
		
		call block0(sd1,5,2,0,0,1,-1.0_IDP)
		call block0(sd2,5,2,1,0,0,1.0_IDP)

!   EP FLR effects

		if (epflr_on == 1 .and. r_epflr > 0.0) then
			sd1=-dnfeqdr*rinv*cureq/(feq-qqinv*cureq)
			sd2=-dnfeqdr*feq/(feq-qqinv*cureq)
			call block0(sd1,5,iw,0,0,1,-1.0_IDP)
			call block0(sd2,5,iw,1,0,0,1.0_IDP)
			sd1=epsq*omcyd*omegar*nfeq/vfova2
			call block0(sd1,5,iw,0,0,0,1.0_IDP)
		end if

!  Load remaining terms in fast ion parallel velocity equation

!    Landau closure term

		do l=1,leqmax
			sceq1(:,l)=1.414213*LcA1*vfova*bmod(:,l)/(feq-qqinv*cureq)
		end do
		call blockj_landau_grad_parallel(sceq1,1,6,6,0,0,0,1.0_IDP)

!    Parallel gradient term

		do l=1,leqmax
			sceq1(:,l)=2*LcA0*vfova2*bmod(:,l)/(nfeq*(feq-qqinv*cureq))
		end do
		call om(sceq1,1,6,5,0,0,0,-1.0_IDP)

!    Omega* term

		if(epflr_on .eq. 0 .or. r_epflr == 0.0) then 

			sd1=vfova2*dnfeqdr*rinv*cureq/(nfeq*(feq-qqinv*cureq))
			sd2=vfova2*dnfeqdr*feq/(nfeq*(feq-qqinv*cureq))
			call block0(sd1,6,1,0,0,1,-1.0_IDP)
			call block0(sd2,6,1,1,0,0,1.0_IDP)
		 
		else 

!   EP FLR effects

			sd1=vfova2*dnfeqdr*rinv*cureq/(nfeq*(feq-qqinv*cureq))
			sd2=vfova2*dnfeqdr*rinv*feq/(nfeq*(feq-qqinv*cureq))
			call block0(sd1,6,ix1,0,0,0,-1.0_IDP)
			call block0(sd2,6,ix2,0,0,0,1.0_IDP)

!   EP FLR effects auxiliary equations

			sd1=1.0_IDP
			call block0(sd1,iw,iw,0,0,0,1.0_IDP)
			coef=-r_epflr*r_epflr
			call block_dlsq(iw,iw,coef)
			call block_dlsq(iw,2,-coef)

			call block0(sd1,ix1,ix1,0,0,0,1.0_IDP)
			call block_dlsq(ix1,ix1,coef)
			call block0(sd1,ix1,1,0,0,1,-1.0_IDP)

			call block0(sd1,ix2,ix2,0,0,0,1.0_IDP)
			call block_dlsq(ix2,ix2,coef)
			call block0(r,ix2,1,1,0,0,-1.0_IDP)
		 
		end if 

! diffusion term added
		sd1=1.0_IDP
!		call block0_dlsq(sd1,6,6,stdifvf)
		call block0_dlsqr(sd1,6,6,stdifvfr,sd7)

!  End of fast ion moment equations

!  Load terms of the thermal moment of energetic particles

!  Pressure gradient term

		do l=1,leqmax
			sceq1(:,l)=bet0*bmod(:,l)/(2.*denseq*(feq-qqinv*cureq))
		end do
		call om(sceq1,1,7,3,0,0,0,-1.0_IDP)
		
!  Magnetic field perturbation term		
		
		call dbydr0(sd1,preq,0.0_IDP,1.0_IDP,0)		
		do l=1,leqmax
			sceq2(:,l)=sceq1(:,l)*sd1
		end do		
		call blockj(sceq2,1,7,1,1,0,0,1.0_IDP)		
		
!  Shared equilibrium toroidal flow velocity for the thermal moment of energetic particles equation

		call block0(vzt_eq,7,7,0,0,1,-1.0_IDP)

! diffusion term added
		sd1=1.0_IDP
!		call block0_dlsq(sd1,7,7,stdifv)
		call block0_dlsqr(sd1,7,7,stdifvr,sd7)
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	Alpha particle effects	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

		if (alpha_on == 1) then

!  Load 1st Omega-d terms in fast ion density equation:

			do l=1,leqmax
				sceq1(:,l)=valphaova2*omdr(:,l)/(epsq*omcydalp)
				sceq2(:,l)=valphaova2*omdt(:,l)/(epsq*omcydalp)
				sceq3(:,l)=valphaova2*omdz(:,l)/(epsq*omcydalp)
			end do
		
			call blockj(sceq1,-1,8,8,0,1,0,-1.0_IDP)
			call blockj(sceq2,1,8,8,1,0,0,-1.0_IDP)
			call blockj(sceq3,1,8,8,0,0,1,-1.0_IDP)

!  Shared equilibrium toroidal flow velocity for fast ion density equation

			call block0(vzt_eq,8,8,0,0,1,-1.0_IDP)

!  Load 1st Omega-d terms in fast ion parallel velocity equation:

			call blockj(sceq1,-1,9,9,0,1,0,-1.0_IDP)
			call blockj(sceq2,1,9,9,1,0,0,-1.0_IDP)
			call blockj(sceq3,1,9,9,0,0,1,-1.0_IDP)

!  Shared equilibrium toroidal flow velocity for fast ion parallel velocity equation

			call block0(vzt_eq,9,9,0,0,1,-1.0_IDP)

!  Load 2nd Omega-d terms in fast ion density equation:

			do l=1,leqmax
				sceq1(:,l)=nalpeq(:)*omdr(:,l)
				sceq2(:,l)=nalpeq(:)*omdt(:,l)
				sceq3(:,l)=nalpeq(:)*omdz(:,l)
			end do

			call blockj(sceq1,-1,8,2,0,1,0,-1.0_IDP)
			call blockj(sceq2,1,8,2,1,0,0,-1.0_IDP)
			call blockj(sceq3,1,8,2,0,0,1,-1.0_IDP)

! diffusion term added
			sd1=1.0_IDP
			call block0_dlsq(sd1,8,8,stdifnalp)

!  Load remaining terms in fast ion density equation

!    Parallel gradient term

			do l=1,leqmax
				sceq1(:,l)=nalpeq*bmod(:,l)/(feq-qqinv*cureq)
			end do
		
			call om(sceq1,1,8,9,0,0,0,-1.0_IDP)

!    Omega* term

			sd1=dnalpeqdr*rinv*cureq/(feq-qqinv*cureq)
			sd2=dnalpeqdr*feq/(feq-qqinv*cureq)
		
			call block0(sd1,8,2,0,0,1,-1.0_IDP)
			call block0(sd2,8,2,1,0,0,1.0_IDP)
		 

!   EP FLR effects

			if (epflr_on == 1 .and. r_epflralp > 0.0) then
				sd1=-dnalpeqdr*rinv*cureq/(feq-qqinv*cureq)
				sd2=-dnalpeqdr*feq/(feq-qqinv*cureq)
				call block0(sd1,8,iwa,0,0,1,-1.0_IDP)
				call block0(sd2,8,iwa,1,0,0,1.0_IDP)
				sd1=epsq*omcydalp*omegar*nalpeq/valphaova2
				call block0(sd1,8,iwa,0,0,0,1.0_IDP)
			end if

!  Load remaining terms in fast ion parallel velocity equation

!    Landau closure term

			do l=1,leqmax
				sceq1(:,l)=1.414213*LcA1alp*valphaova*bmod(:,l)/(feq-qqinv*cureq)
			end do
			call blockj_landau_grad_parallel(sceq1,1,9,9,0,0,0,1.0_IDP)

!    Parallel gradient term

			do l=1,leqmax
				sceq1(:,l)=2*LcA0alp*valphaova2*bmod(:,l)/(nalpeq*(feq-qqinv*cureq))
			end do
			call om(sceq1,1,9,8,0,0,0,-1.0_IDP)

!    Omega* term

			if(epflr_on .eq. 0 .or. r_epflralp == 0.0) then 

				sd1=valphaova2*dnalpeqdr*rinv*cureq/(nalpeq*(feq-qqinv*cureq))
				sd2=valphaova2*dnalpeqdr*feq/(nalpeq*(feq-qqinv*cureq))
				call block0(sd1,9,1,0,0,1,-1.0_IDP)
				call block0(sd2,9,1,1,0,0,1.0_IDP)

			else

!   EP FLR effects

				sd1=valphaova2*dnalpeqdr*rinv*cureq/(nalpeq*(feq-qqinv*cureq))
				sd2=valphaova2*dnalpeqdr*rinv*feq/(nalpeq*(feq-qqinv*cureq))
				call block0(sd1,9,ix1a,0,0,0,-1.0_IDP)
				call block0(sd2,9,ix2a,0,0,0,1.0_IDP)

!   EP FLR effects auxiliary equations

				sd1=1.0_IDP
				call block0(sd1,iwa,iwa,0,0,0,1.0_IDP)
				coef=-r_epflralp*r_epflralp
				call block_dlsq(iwa,iwa,coef)
				call block_dlsq(iwa,2,-coef)

				call block0(sd1,ix1a,ix1a,0,0,0,1.0_IDP)
				call block_dlsq(ix1a,ix1a,coef)
				call block0(sd1,ix1a,1,0,0,1,-1.0_IDP)

				call block0(sd1,ix2a,ix2a,0,0,0,1.0_IDP)
				call block_dlsq(ix2a,ix2a,coef)
				call block0(r,ix2a,1,1,0,0,-1.0_IDP)
		 
			end if 

! diffusion term added
			sd1=1.0_IDP
			call block0_dlsq(sd1,9,9,stdifvalp)

!  End of alpha particles moment equations

		end if

		amat=-dtd2*amat
		bmat=-dtd2*bmat
		cmat=-dtd2*cmat

!  put in l.h.s. of equation

!  psi equation

		sd1=1.0_IDP
		call block0(sd1,1,1,0,0,0,1.0_IDP)

!  u-zeta equation

		call block0(sd1,2,4,0,0,0,1.0_IDP)

! p equation

		call block0(sd1,3,3,0,0,0,1.0_IDP)

! fast ion density moment equation

		call block0(sd1,5,5,0,0,0,1.0_IDP)

! fast ion parallel velocity moment equation

		call block0(sd1,6,6,0,0,0,1.0_IDP)
		
! thermal moment equation

		call block0(sd1,7,7,0,0,0,1.0_IDP)		

		if (alpha_on == 1) then

! fast ion density moment equation

			call block0(sd1,8,8,0,0,0,1.0_IDP)

! fast ion parallel velocity moment equation

			call block0(sd1,9,9,0,0,0,1.0_IDP)
		
		end if

		yt=0.0_IDP
		do i=nr_start,nr_end
			mnum3=noeqn*mnumn(i)
			do l1t=1,mnum3
				do l2t=1,mnum3
					do j=1,mjm1
						imt1=l1t+mnum3*(j-1)+nskpxn(i)
						imt2=l2t+mnum3*(j-1)+nskpxn(i)
						imat=l1t+mnum3*(l2t-1+mnum3*(j-1))+nskpn(i)
						if (j > 1) yt(imt1)=yt(imt1)+cmat(imat)*xt(imt2-mnum3)
						yt(imt1)=yt(imt1)+amat(imat)*xt(imt2)
						yt(imt1)=yt(imt1)+bmat(imat)*xt(imt2+mnum3)
					end do
				end do
			end do
		end do
		do i=nc_start,nc_end
			mnum3=noeqn*mnumn(i)
			do l1t=1,mnum3
				do l2t=1,mnum3
					do j=1,mjm1
						imt1=l1t+2*mnum3*(j-1)+nskpxn(i)
						imt2=l2t+2*mnum3*(j-1)+nskpxn(i)
						imat=l1t+2*mnum3*(l2t-1+mnum3*(j-1))+nskpn(i)
						if (j > 1) then
							yt(imt1)=yt(imt1)+cmat(imat)*xt(imt2-2*mnum3)-cmat(imat+mnum3)*xt(imt2-mnum3)
							yt(imt1+mnum3)=yt(imt1+mnum3)+cmat(imat)*xt(imt2-mnum3)+cmat(imat+mnum3)*xt(imt2-2*mnum3)
						end if
						yt(imt1)=yt(imt1)+amat(imat)*xt(imt2)-amat(imat+mnum3)*xt(imt2+mnum3)
						yt(imt1+mnum3)=yt(imt1+mnum3)+amat(imat)*xt(imt2+mnum3)+amat(imat+mnum3)*xt(imt2)
						yt(imt1)=yt(imt1)+bmat(imat)*xt(imt2+2*mnum3)-bmat(imat+mnum3)*xt(imt2+3*mnum3)
						yt(imt1+mnum3)=yt(imt1+mnum3)+bmat(imat)*xt(imt2+3*mnum3)+bmat(imat+mnum3)*xt(imt2+2*mnum3)
					end do
				end do
			end do
		end do

!  u-dlpersq(phi)=0

		do i=nr_start,nr_end
			mnum3=noeqn*mnumn(i)
			lskp=3*mnumn(i)
			do l1=1,mnumn(i)
				do j=1,mjm1
					l1t=l1+lskp+mnum3*(j-1)+nskpxn(i)
					yt(l1t)=0.
				end do
			end do
		end do
		do i=nc_start,nc_end
			mnum3=noeqn*mnumn(i)
			lskp=3*mnumn(i)
			do l1=1,mnumn(i)
				do j=1,mjm1
					l1t=l1+lskp+2*mnum3*(j-1)+nskpxn(i)
					yt(l1t)=0.
					yt(l1t+mnum3)=0.
				end do
			end do
		end do

!   Ion FLR effects auxiliary equation

		if (iflr_on == 1) then
			do i=nr_start,nr_end
				mnum3=noeqn*mnumn(i)
				lskp=(iq-1)*mnumn(i)
				do l1=1,mnumn(i)
					do j=1,mjm1
						l1t=l1+lskp+mnum3*(j-1)+nskpxn(i)
						yt(l1t)=0.
					end do
				end do
			end do
			do i=nc_start,nc_end
				mnum3=noeqn*mnumn(i)
				lskp=(iq-1)*mnumn(i)
				do l1=1,mnumn(i)
					do j=1,mjm1
						l1t=l1+lskp+2*mnum3*(j-1)+nskpxn(i)
						yt(l1t)=0.
						yt(l1t+mnum3)=0.
					end do
				end do
			end do
		end if

!   EP FLR effects auxiliary equations

		if (epflr_on == 1) then
			if (r_epflr > 0.0) then
				do i=nr_start,nr_end
					mnum3=noeqn*mnumn(i)
					lskp=(iw-1)*mnumn(i)
					do l1=1,3*mnumn(i)
						do j=1,mjm1
							l1t=l1+lskp+mnum3*(j-1)+nskpxn(i)
							yt(l1t)=0.
						end do
					end do
				end do
				do i=nc_start,nc_end
					mnum3=noeqn*mnumn(i)
					lskp=(iw-1)*mnumn(i)
					do l1=1,3*mnumn(i)
						do j=1,mjm1
							l1t=l1+lskp+2*mnum3*(j-1)+nskpxn(i)
							yt(l1t)=0.
							yt(l1t+mnum3)=0.
						end do
					end do
				end do
			end if
			if (alpha_on == 1 .and. r_epflralp > 0.0) then
				do i=nr_start,nr_end
					mnum3=noeqn*mnumn(i)
					lskp=(iwa-1)*mnumn(i)
					do l1=1,3*mnumn(i)
						do j=1,mjm1
							l1t=l1+lskp+mnum3*(j-1)+nskpxn(i)
							yt(l1t)=0.
						end do
					end do
				end do
				do i=nc_start,nc_end
					mnum3=noeqn*mnumn(i)
					lskp=(iwa-1)*mnumn(i)
					do l1=1,3*mnumn(i)
						do j=1,mjm1
							l1t=l1+lskp+2*mnum3*(j-1)+nskpxn(i)
							yt(l1t)=0.
							yt(l1t+mnum3)=0.
						end do
					end do
				end do
			end if
		end if

		do i=nr_start,nr_end
			mnum6=noeqn*noeqn*mnumn(i)*mnumn(i)
			do l1=1,mnum6
				l1t=l1+nskpn(i)
				l2t=l1+mnum6*(mj-2)+nskpn(i)
				bmat(l2t)=0.
				cmat(l1t)=0.
			end do
		end do
		do i=nc_start,nc_end
			mnum6=2*noeqn*noeqn*mnumn(i)*mnumn(i)
			do l1=1,mnum6
				l1t=l1+nskpn(i)
				l2t=l1+mnum6*(mj-2)+nskpn(i)
				bmat(l2t)=0.
				cmat(l1t)=0.
			end do
		end do

		do it=nr_start,nr_end
			loca=nskpn(it)+1
			loci=nskpin(it)+1
			mnum3=noeqn*mnumn(it)
			call decbt(mnum3,mjm1,amat(loca:),bmat(loca:),cmat(loca:),ipc(loci:),ier)
			if (ier /= 0) then
				write(6,'("ier=",i5," ntor=",i5)') ier,it
				stop 19
			end if
		end do
		do it=nc_start,nc_end
			loca=nskpn(it)+1
			loci=nskpin(it)+1
			mnum3=noeqn*mnumn(it)
			call cdecbt(mnum3,mjm1,amat(loca:),bmat(loca:),cmat(loca:),ipc(loci:),ier)
			if (ier /= 0) then
				write(6,'("ier=",i5," ntor=",i5)') ier,it
				stop 19
			end if
		end do

	end subroutine linstart