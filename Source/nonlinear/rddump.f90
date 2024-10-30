	subroutine rddump

		use mpi
		use param
		use processor
		use var_para
		use cotrol
		use domain
		use equil
		use dynamo
		implicit none

		integer :: i,l,j,ihisto
		real(IDP) :: pertsclo,dum
		integer, dimension(ldim) :: lmap
		character(len=8) :: confil		

!		The data of the previous run is read to continue the simulation	

		read(8) ihisto,numruno,numrunp,numruns,numeq,numvac,(numhist(i),i=1,ihisto),maxstp,nstep, &
			ndump,nprint,lplots,ndiag,jmn,jmx,jint,itime,dt0,nstep1,nonlin,m0dy,nocpl,ifreezd,ifreezu, &
			Adens,Bdens,LcA0,LcA1,LcA2,LcA3,omcy,dpres,ext_prof,epflr_on,r_epflr,alpha_on,iflr_on,iflr,twofl_on,ieldamp_on, &
			omegar,bet0_alp,Adensalp,Bdensalp,LcA0alp,LcA1alp,LcA2alp,LcA3alp,omcyalp,r_epflralp, &
                         namelist_on,Trapped_on,ext_prof_name,ext_prof_len,Auto_grid_on,Edge_on,edge_p, &
                         betath_factor,spe1,spe2,omcyb,rbound,EP_dens_on,Alpha_dens_on,EP_vel_on,Alpha_vel_on, &
                         q_prof_on,Eq_vel_on,Eq_velp_on,Eq_Presseq_on,Eq_Presstot_on,deltaq,deltaiota
		read(8) mj,lmax,leqmax,mjm1
		read(8) (r(j),j=0,mj),(mm(l),l=1,lmax),(nn(l),l=1,lmax),(mmeq(l),l=1,leqmax),(nneq(l),l=1,leqmax),(rinv(j),j=0,mj), &		
			(dc1m(j),j=1,mj),(dc1p(j),j=1,mj),(dc2m(j),j=1,mj),(dc2p(j),j=1,mj),(del2cm(j),j=1,mj),(del2cp(j),j=1,mj)	

		allocate (sqgi(0:mj,0:leqmax),sqg(0:mj,0:leqmax),bst(0:mj,0:leqmax),grr(0:mj,0:leqmax),grt(0:mj,0:leqmax), &
			  gtt(0:mj,0:leqmax),grroj(0:mj,0:leqmax),grtoj(0:mj,0:leqmax),gttoj(0:mj,0:leqmax),jbgrr(0:mj,0:leqmax), &
			  jbgrt(0:mj,0:leqmax),jbgtt(0:mj,0:leqmax),lplrr(0:mj,0:leqmax),lplrt(0:mj,0:leqmax),lplrz(0:mj,0:leqmax), &
			  lpltt(0:mj,0:leqmax),lpltz(0:mj,0:leqmax),lplzz(0:mj,0:leqmax),lplr(0:mj,0:leqmax),lplt(0:mj,0:leqmax), &
			  lplz(0:mj,0:leqmax),djroj(0:mj,0:leqmax),djtoj(0:mj,0:leqmax),djzoj(0:mj,0:leqmax),omdr(0:mj,0:leqmax), &
			  omdt(0:mj,0:leqmax),omdz(0:mj,0:leqmax),bmod(0:mj,0:leqmax),dbsjtoj(0:mj,0:leqmax),dbsjzoj(0:mj,0:leqmax))

                if(Trapped_on .eq. 1) then	 
		   allocate (omdrprp(0:mj,0:leqmax))
		   allocate (omdtprp(0:mj,0:leqmax))
		   allocate (omdzprp(0:mj,0:leqmax))	
		end if	
		
		sqgi=0.0_IDP
		sqg=0.0_IDP
		bst=0.0_IDP
		grr=0.0_IDP
		grt=0.0_IDP
		gtt=0.0_IDP
		grroj=0.0_IDP
		grtoj=0.0_IDP
		gttoj=0.0_IDP
		jbgrr=0.0_IDP
		jbgrt=0.0_IDP
		jbgtt=0.0_IDP
		lplrr=0.0_IDP
		lplrt=0.0_IDP
		lplrz=0.0_IDP
		lpltt=0.0_IDP
		lpltz=0.0_IDP
		lplzz=0.0_IDP
		lplr=0.0_IDP
		lplt=0.0_IDP
		lplz=0.0_IDP
		djroj=0.0_IDP
		djtoj=0.0_IDP
		djzoj=0.0_IDP
		omdr=0.0_IDP
		omdt=0.0_IDP
		omdz=0.0_IDP	
		bmod=0.0_IDP
		dbsjtoj=0.0_IDP
		dbsjzoj=0.0_IDP

                if(Trapped_on .eq. 1) then
		   omdrprp=0.0_IDP
		   omdtprp=0.0_IDP
		   omdzprp=0.0_IDP
		end if	

			
		read(8) ni,nis,ne,delta,rc,fti,fte
		read(8) (qq(j),j=0,mj),(qqinv(j),j=0,mj),(psieq(j),j=0,mj),(chieq(j),j=0,mj),(preq(j),j=0,mj),(feq(j),j=0,mj), &
			(cureq(j),j=0,mj),(denseq(j),j=0,mj),(teeq(j),j=0,mj),(nfeq(j),j=0,mj),(vfova(j),j=0,mj),(vzt_eq(j),j=0,mj), &
			(tieq(j),j=0,mj),(vtherm_ion(j),j=1,mj)
		if (alpha_on == 1) read(8) (nalpeq(j),j=0,mj),(valphaova(j),j=0,mj)
		read(8) ngeneq,ndevice(1),ndevice(2),ieq,eps,bet0,bet0_f,nqprof,q0,qlamb,ro,(rs(l),l=1,lmax),leq,qh,rsbar,qlambpp, &
			ndat,(xr(j),j=1,ndat),(xdat(j),j=1,ndat),pop,popp,pol2,ltwo,betaj,(pdat(j),j=1,ndat),ntok,betat,betap
		read(8) ((sqgi(j,l),j=0,mj),l=1,leqmax),((sqg(j,l),j=0,mj),l=1,leqmax),((bst(j,l),j=0,mj),l=1,leqmax), &
			((grr(j,l),j=0,mj),l=1,leqmax),((grt(j,l),j=0,mj),l=1,leqmax),((gtt(j,l),j=0,mj),l=1,leqmax), &
			((grroj(j,l),j=0,mj),l=1,leqmax),((grtoj(j,l),j=0,mj),l=1,leqmax),((gttoj(j,l),j=0,mj),l=1,leqmax), &
			((jbgrr(j,l),j=0,mj),l=1,leqmax),((jbgrt(j,l),j=0,mj),l=1,leqmax),((jbgtt(j,l),j=0,mj),l=1,leqmax), &
			((lplrr(j,l),j=0,mj),l=1,leqmax),((lplrt(j,l),j=0,mj),l=1,leqmax),((lplrz(j,l),j=0,mj),l=1,leqmax), &
			((lpltt(j,l),j=0,mj),l=1,leqmax),((lpltz(j,l),j=0,mj),l=1,leqmax),((lplzz(j,l),j=0,mj),l=1,leqmax), &
			((lplr(j,l),j=0,mj),l=1,leqmax),((lplt(j,l),j=0,mj),l=1,leqmax),((lplz(j,l),j=0,mj),l=1,leqmax), &
			((djroj(j,l),j=0,mj),l=1,leqmax),((djtoj(j,l),j=0,mj),l=1,leqmax),((djzoj(j,l),j=0,mj),l=1,leqmax), &
			((omdr(j,l),j=0,mj),l=1,leqmax),((omdt(j,l),j=0,mj),l=1,leqmax),((omdz(j,l),j=0,mj),l=1,leqmax), &
			((bmod(j,l),j=0,mj),l=1,leqmax),((dbsjtoj(j,l),j=0,mj),l=1,leqmax),((dbsjzoj(j,l),j=0,mj),l=1,leqmax)
		if (ieldamp_on == 1) then
			allocate (eildr(0:mj,0:leqmax),eildt(0:mj,0:leqmax),eildz(0:mj,0:leqmax),eildrr(0:mj,0:leqmax), &
				  eildrt(0:mj,0:leqmax),eildrz(0:mj,0:leqmax),eildtt(0:mj,0:leqmax),eildtz(0:mj,0:leqmax), &
				  eildzz(0:mj,0:leqmax))
			read(8) ((eildr(j,l),j=0,mj),l=1,leqmax),((eildt(j,l),j=0,mj),l=1,leqmax), &
				((eildz(j,l),j=0,mj),l=1,leqmax),((eildrr(j,l),j=0,mj),l=1,leqmax), &
				((eildrt(j,l),j=0,mj),l=1,leqmax),((eildrz(j,l),j=0,mj),l=1,leqmax), &
				((eildtt(j,l),j=0,mj),l=1,leqmax),((eildtz(j,l),j=0,mj),l=1,leqmax), &
				((eildzz(j,l),j=0,mj),l=1,leqmax)
		end if

                if(Trapped_on .eq. 1) then
		  read(8)  ((omdrprp(j,l),j=0,mj),l=1,leqmax),((omdtprp(j,l),j=0,mj),l=1,leqmax),((omdzprp(j,l),j=0,mj),l=1,leqmax)
		end if	

		read(8) (eta(j),j=0,mj),dt,time

		if (myPE == 0) then
			read(8) ((psi(j,l),j=0,mj),l=1,lmax)
			read(8) ((phi(j,l),j=0,mj),l=1,lmax)
			read(8) ((pr(j,l),j=0,mj),l=1,lmax)
			read(8) ((nf(j,l),j=0,mj),l=1,lmax)
			read(8) ((vprlf(j,l),j=0,mj),l=1,lmax)
			read(8) ((vthprlf(j,l),j=0,mj),l=1,lmax)
			if (alpha_on == 1) then
				read(8) ((nalp(j,l),j=0,mj),l=1,lmax)
				read(8) ((vprlalp(j,l),j=0,mj),l=1,lmax)
			end if
		else
			read(8) ((dum,j=0,mj_start-1),(psi(j,l),j=mj_start,mj_end),(dum,j=mj_end+1,mj),l=1,lmax)
			read(8) ((dum,j=0,mj_start-1),(phi(j,l),j=mj_start,mj_end),(dum,j=mj_end+1,mj),l=1,lmax)
			read(8) ((dum,j=0,mj_start-1),(pr(j,l),j=mj_start,mj_end),(dum,j=mj_end+1,mj),l=1,lmax)
			read(8) ((dum,j=0,mj_start-1),(nf(j,l),j=mj_start,mj_end),(dum,j=mj_end+1,mj),l=1,lmax)
			read(8) ((dum,j=0,mj_start-1),(vprlf(j,l),j=mj_start,mj_end),(dum,j=mj_end+1,mj),l=1,lmax)
			read(8) ((dum,j=0,mj_start-1),(vthprlf(j,l),j=mj_start,mj_end),(dum,j=mj_end+1,mj),l=1,lmax)
			if (alpha_on == 1) then
				read(8) ((dum,j=0,mj_start-1),(nalp(j,l),j=mj_start,mj_end),(dum,j=mj_end+1,mj),l=1,lmax)
				read(8) ((dum,j=0,mj_start-1),(vprlalp(j,l),j=mj_start,mj_end),(dum,j=mj_end+1,mj),l=1,lmax)
			end if
		end if
		
		read(8) etascl,reta,eta0,etalmb,ietaeq,(etadat(j),j=1,ndat),chipr,ichipr,rchi,stdifp,stdifu,stdifnf,stdifv,stdifvf, &
			stdifnalp,stdifvalp,s,gamma,xle,ipert

		rewind(8)
		close(unit=8)
			
	end subroutine rddump