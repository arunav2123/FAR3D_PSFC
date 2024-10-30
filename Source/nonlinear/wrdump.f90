	subroutine wrdump

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		implicit none

		integer :: i,l,j
		character(len=8) :: confil

!		The data of the run is stored in case there is a continuation run

		confil="fs"//numrun(1)//numrun(2)//numrun(3)
                if (nonlin == 1) then
                  open(unit=7,file=confil,status='new',form='unformatted')
                else
		  open(unit=7,file=confil,status='unknown',form='unformatted')
                end if

		write(7) ihist,numrun,numruno,numruns,numeq,numvac,(numhist(i),i=1,ihist),maxstp,nstep, &
			 ndump,nprint,lplots,ndiag,jmn,jmx,jint,itime,dt0,nstep1,nonlin,m0dy,nocpl,ifreezd,ifreezu, &
			 Adens,Bdens,LcA0,LcA1,LcA2,LcA3,omcy,dpres,ext_prof,epflr_on,r_epflr,alpha_on,iflr_on,iflr,twofl_on,ieldamp_on, &
			 omegar,bet0_alp,Adensalp,Bdensalp,LcA0alp,LcA1alp,LcA2alp,LcA3alp,omcyalp,r_epflralp, &
                         namelist_on,Trapped_on,ext_prof_name,ext_prof_len,Auto_grid_on,Edge_on,edge_p, &
                         betath_factor,spe1,spe2,omcyb,rbound,EP_dens_on,Alpha_dens_on,EP_vel_on,Alpha_vel_on, &
                         q_prof_on,Eq_vel_on,Eq_velp_on,Eq_Presseq_on,Eq_Presstot_on,deltaq,deltaiota
		write(7) mj,lmax,leqmax,mjm1
		write(7) (r(j),j=0,mj),(mm(l),l=1,lmax),(nn(l),l=1,lmax),(mmeq(l),l=1,leqmax),(nneq(l),l=1,leqmax),(rinv(j),j=0,mj), &		
			 (dc1m(j),j=1,mj),(dc1p(j),j=1,mj),(dc2m(j),j=1,mj),(dc2p(j),j=1,mj),(del2cm(j),j=1,mj),(del2cp(j),j=1,mj)
		write(7) ni,nis,ne,delta,rc,fti,fte
		write(7) (qq(j),j=0,mj),(qqinv(j),j=0,mj),(psieq(j),j=0,mj),(chieq(j),j=0,mj),(preq(j),j=0,mj),(feq(j),j=0,mj), &
			 (cureq(j),j=0,mj),(denseq(j),j=0,mj),(teeq(j),j=0,mj),(nfeq(j),j=0,mj),(vfova(j),j=0,mj),(vzt_eq(j),j=0,mj), &
			 (tieq(j),j=0,mj),(vtherm_ion(j),j=1,mj)
                
		if (alpha_on == 1) write(7) (nalpeq(j),j=0,mj),(valphaova(j),j=0,mj)
		write(7) ngeneq,ndevice(1),ndevice(2),ieq,eps,bet0,bet0_f,nqprof,q0,qlamb,ro,(rs(l),l=1,lmax),leq,qh,rsbar,qlambpp, &
			 ndat,(xr(j),j=1,ndat),(xdat(j),j=1,ndat),pop,popp,pol2,ltwo,betaj,(pdat(j),j=1,ndat),ntok,betat,betap

		write(7) ((sqgi(j,l),j=0,mj),l=1,leqmax),((sqg(j,l),j=0,mj),l=1,leqmax),((bst(j,l),j=0,mj),l=1,leqmax), &
				 ((grr(j,l),j=0,mj),l=1,leqmax),((grt(j,l),j=0,mj),l=1,leqmax),((gtt(j,l),j=0,mj),l=1,leqmax), &
				 ((grroj(j,l),j=0,mj),l=1,leqmax),((grtoj(j,l),j=0,mj),l=1,leqmax),((gttoj(j,l),j=0,mj),l=1,leqmax), &	 
				 ((jbgrr(j,l),j=0,mj),l=1,leqmax),((jbgrt(j,l),j=0,mj),l=1,leqmax),((jbgtt(j,l),j=0,mj),l=1,leqmax), &
				 ((lplrr(j,l),j=0,mj),l=1,leqmax),((lplrt(j,l),j=0,mj),l=1,leqmax),((lplrz(j,l),j=0,mj),l=1,leqmax), &
				 ((lpltt(j,l),j=0,mj),l=1,leqmax),((lpltz(j,l),j=0,mj),l=1,leqmax),((lplzz(j,l),j=0,mj),l=1,leqmax), &
				 ((lplr(j,l),j=0,mj),l=1,leqmax),((lplt(j,l),j=0,mj),l=1,leqmax),((lplz(j,l),j=0,mj),l=1,leqmax), &	
				 ((djroj(j,l),j=0,mj),l=1,leqmax),((djtoj(j,l),j=0,mj),l=1,leqmax),((djzoj(j,l),j=0,mj),l=1,leqmax), &				 
				 ((omdr(j,l),j=0,mj),l=1,leqmax),((omdt(j,l),j=0,mj),l=1,leqmax),((omdz(j,l),j=0,mj),l=1,leqmax), &	
				 ((bmod(j,l),j=0,mj),l=1,leqmax),((dbsjtoj(j,l),j=0,mj),l=1,leqmax),((dbsjzoj(j,l),j=0,mj),l=1,leqmax)

		if (ieldamp_on == 1) write(7) ((eildr(j,l),j=0,mj),l=1,leqmax),((eildt(j,l),j=0,mj),l=1,leqmax), &
					      ((eildz(j,l),j=0,mj),l=1,leqmax),((eildrr(j,l),j=0,mj),l=1,leqmax), &
					      ((eildrt(j,l),j=0,mj),l=1,leqmax),((eildrz(j,l),j=0,mj),l=1,leqmax), &
					      ((eildtt(j,l),j=0,mj),l=1,leqmax),((eildtz(j,l),j=0,mj),l=1,leqmax), &
					      ((eildzz(j,l),j=0,mj),l=1,leqmax)

                if(Trapped_on .eq. 1) then
		  write(7)  ((omdrprp(j,l),j=0,mj),l=1,leqmax),((omdtprp(j,l),j=0,mj),l=1,leqmax),((omdzprp(j,l),j=0,mj),l=1,leqmax)
		end if	

		write(7) (eta(j),j=0,mj),dt,time
		
		write(7) ((psi(j,l),j=0,mj),l=1,lmax)
		write(7) ((phi(j,l),j=0,mj),l=1,lmax)
		write(7) ((pr(j,l),j=0,mj),l=1,lmax)
		write(7) ((nf(j,l),j=0,mj),l=1,lmax)
		write(7) ((vprlf(j,l),j=0,mj),l=1,lmax)
		write(7) ((vthprlf(j,l),j=0,mj),l=1,lmax)
		if (alpha_on == 1) then
			write(7) ((nalp(j,l),j=0,mj),l=1,lmax)
			write(7) ((vprlalp(j,l),j=0,mj),l=1,lmax)
		end if

		write(7) etascl,reta,eta0,etalmb,ietaeq,(etadat(j),j=1,ndat),chipr,ichipr,rchi,stdifp,stdifu,stdifnf,stdifv,stdifvf, &
			 stdifnalp,stdifvalp,s,gamma,xle,ipert,(widthi(l),l=1,lmax),(gammai(l),l=1,lmax),pertscl

		rewind(7)
		close(unit=7)
		write (6,'(/"  wrdump: have written fs",2a2,a1)') (numrun(i),i=1,3)

	end subroutine wrdump