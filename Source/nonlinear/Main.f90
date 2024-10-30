!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Open-source MIT LINCENSE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!    FAR3d gyrofluid code ver 2.0 (2 EP real nonlinear version)                !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!    Copyright (C) 2024  D. Spong, L. Garcia, J. Varela and Y. Ghai            !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                  						  !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!    This program is free software: you can redistribute it and/or modify      !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!    it under the terms of the open-source MIT license.                        !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!										  !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!    This program is distributed in the hope that it will be useful,           !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!    but WITHOUT ANY WARRANTY; without even the implied warranty of            !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.   		          !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	program far

		use mpi
		use param
		use processor
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		namelist/nam_par/maxstp,ndump,ndiag,nprint,lplots,itime,dt0,nonlin,mj,lmax,leqmax,jmn,jmx,jint,ni,nis,ne, &
			     delta,rc,fti,fte,stdifp,stdifu,stdifv,ngeneq,ieq,eps,epst,bet0,nqprof,q0,qlamb,ro,rs,leq,qh,rsbar,qlambpp, &
			     ndat,pop,popp,pol2,ltwo,betaj,ntok,etascl,reta,eta0,etalmb,ietaeq,chipr,ichipr,rchi, &
			     s,gamma,ipert,pertscl,m0dy,nocpl,ifreezd,ifreezu,xle,wsrc,rsrc,asrc,omcy,bet0_f, &
			     stdifnf,stdifvf,ext_prof,epflr_on,r_epflr,alpha_on,iflr_on,iflr,LcA0,LcA1,LcA2,LcA3,Adens,Bdens, &
			     twofl_on,dpres,bet0_alp,Adensalp,Bdensalp,LcA0alp,LcA1alp,LcA2alp,LcA3alp,omcyalp,r_epflralp,stdifnalp,stdifvalp, &
			     omegar,ieldamp_on,ndis,nwdth,omcyb,rbound,trapped_on,namelist_on,betath_factor,spe1,spe2,EP_dens_on, &
                             EP_vel_on,Alpha_dens_on,Alpha_vel_on,DIIID_u,Eq_vel_on,Eq_velp_on,q_prof_on,deltaq,deltaiota,Eq_Presseq_on, &
                             Eq_Presstot_on,Edge_on,edge_p,Auto_grid_on,ext_prof_len,nopsievol_on,noprevol_on,nonfevol_on,nonalpevol_on
		namelist/nam_arr/mm,nn,mmeq,nneq,xr,xdat,pdat,etadat,widthi,gammai,cnep,ctep,cvep,cnfp,cvfp,cnfpalp,cvfpalp,eqvt,eqvp

		character(len=132) :: char5
		character(len=10) :: timew,datew,zone_h
		integer, dimension(8) :: values_s,values_e
		character(len=16) :: confil
		character(len=5) :: numeqc
		character(len=1) :: t
		character(len=10) :: hifrq_out
		integer :: i,l,j,iend,l2t,mnum,ierror
		real(IDP) :: d,dk,scmn,scmx,sclmn,sclmx,scnorm,dum

		interface
			subroutine inputlist
			end subroutine inputlist
			subroutine dfault
			end subroutine dfault
			subroutine inital
			end subroutine inital
			subroutine resume
			end subroutine resume
			subroutine linstart
			end subroutine linstart
			subroutine linstep
			end subroutine linstep
			subroutine lincheck
			end subroutine lincheck
			subroutine output
			end subroutine output
			subroutine wrdump
			end subroutine wrdump
			subroutine endrun
			end subroutine endrun			
			subroutine numinc
			end subroutine numinc
			subroutine hifrq
			end subroutine hifrq
			subroutine trnsfr0(f,ftype)
				use param
				use var_para
				implicit none
				integer :: ftype
				real(IDP), dimension(mj_start:,0:) :: f
			end subroutine trnsfr0
			subroutine energy(i)
				implicit none
				integer :: i
			end subroutine energy
			subroutine elapsed_time(values_s,values_e)
				implicit none
				integer, dimension(8) :: values_s,values_e
			end subroutine elapsed_time
		end interface

!		Input read	 
                                                   
		print *, "====================/ WELCOME TO \===================="	
		print *, "======================================================"			
		print *, "===============================  ______========_======"		
		print *, "=====_____====______====_____===|_____ |=======||====="	
		print *, "====| ____|==| ____ |==| ___ |========||=======||====="	
		print *, "====||____===||====||==||===||===_____||=======||====="	
		print *, "====| ____|==||____||==||___||==|_____ |== ____||====="			
		print *, "====||=======| ____ |==||=\\==========||==| ___ |====="	
		print *, "====||=======||====||==||==\\====_____||==||___||====="	
		print *, "====||=======||====||==||== \\==|______|==|_____|====="	
		print *, "======================================================"	
		print *, "======================\ ver1.0 /======================"		

!		Initialize MPI
		call MPI_INIT(ierror)
!		Find number of PEs
		call MPI_COMM_SIZE(MPI_COMM_WORLD, numPEs, ierror)
!		Find the ID of this PE
		call MPI_COMM_RANK(MPI_COMM_WORLD, myPE  , ierror)
		numPEsm1 = numPEs-1
		
		write(0,'(" ====> Checking input list ... ")')

		call inputlist	
                nstres_com = nstres

		lmax=ldim
		leqmax=leqdim	
		mj=jdim

!	        Date and time	
	
		  call date_and_time(datew,timew,zone_h,values_s)	
		  write(6,'(/" time = ",a2,":",a2,":",a2,"       date = ",a2,"/",a2,"/",a4/)') timew(1:2),timew(3:4),timew(5:6), &
									     datew(7:8),datew(5:6),datew(1:4)

!		Identify the number of equations in the model	

		noeqn=7
		ivalp=0
		iq=0
		iw=0
		ix1=0
		ix2=0
		iwa=0
		ix1a=0
		ix2a=0
		if (alpha_on == 1) then
			noeqn=noeqn+2
			ivalp=noeqn
		end if
		nvar=noeqn
		if (iflr_on == 1) noeqn=noeqn+1
		iq=noeqn
		if (epflr_on == 1) then
			if (r_epflr > 0.0) then
				iw=noeqn+1
				ix1=noeqn+2
				noeqn=noeqn+3
				ix2=noeqn
			end if
			if (alpha_on == 1 .and. r_epflralp > 0.0) then
				iwa=noeqn+1
				ix1a=noeqn+2
				noeqn=noeqn+3
				ix2a=noeqn
			end if
		end if

		if (ext_prof == 0) ieldamp_on = 0		

!		Initialize the perturbed variables

		psi=0.0_IDP
		phi=0.0_IDP
		pr=0.0_IDP
		nf=0.0_IDP
		vprlf=0.0_IDP
		vthprlf=0.0_IDP
		if (alpha_on == 1) then
			nalp=0.0_IDP
			vprlalp=0.0_IDP
		end if

!		nstres indicates if the run is a new run or a continuation
!		If it is a new run the inital subroutine is called
!		If it is a continuation the resume subroutine is called

		if (nstres == 0) call inital

		if (nstres == 0) write(0,'(" ====> Preparing new run ...")')	

		if (nstres /= 0) call resume

		if (nstres /= 0) write(0,'(" ====> Preparing run continuation ...")')	

!		Modification of the equlibria thermal beta		
		if (betath_factor /= 1) write(0,'(" WARNING THERMAL BETA FACTOR ACTIVE = ",1pe13.6)') betath_factor 
		
!		Modification of the safety factor / iota profiles
		if (deltaq /= 0) write(0,'(" WARNING SAFETY FACTOR PROFILE DISPLACEMET ACTIVE = ",1pe13.6)') deltaq
		if (deltaiota /= 0) write(0,'(" WARNING IOTA PROFILE DISPLACEMET ACTIVE = ",1pe13.6)') deltaiota

		if (itime == 1) then
			dk=0.
			do l=1,lmax
				do j=0,mj
					d=nn(l)-mm(l)*qqinv(j)
					dk=max(dk,abs(d))
				end do
			end do
			dt=dt0*2./dk
			dtd2=dt*.5
		else
			dt=dt0
			dtd2=dt*.5
		end if
		if (myPE == 0) then

			confil="fs"//numrun(1)//numrun(2)
			numhist(ihist)=confil

!		The model parameter are included in farprt output file

			call output

!		Generation of high frequency data output

                if(nonlin .ne. 0 .and. myPE == 0) then
	         hifrq_out="hifrq_"//numrun(1)//numrun(2)
                 open(unit=82,file=hifrq_out,status='unknown',form='formatted')
                 if(nstres_com .eq. 0) write(82, &
		 '("time",a1,"bth_r.2",a1,"bth_r.4",a1,"bth_r.6",a1,"bth_r.8", &
		  a1,"bth_r1",a1,"avg(phi)",a1,"avg(psi)",a1,"avg(denf)",a1,"avg(vprlf)")') t,t,t,t,t,t,t,t,t 
                 close(unit=82)
                endif

!		The model parameter are saved in fs####z file	

		        if (nstres == 0) call wrdump
			call numinc
		end if

!		From here the model is advanced in time	

		write(0,'(" ====> Time stepping begins ...")')  

		do i = 1,10
		  nstep_count(i)=maxstp*i/10
		end do
		if (nstres /= 0) nstep_count(:)=nstep_count(:)+nstep		
		i=1

		if (maxstp > 0) then

			nstep1=nstep

!		Subroutine linstart creates the tridiagonal matrix where the right and left side of the
!		model equations are added  

			call linstart

			iend = 0
			do while (nstep < nstep1+maxstp)

				nstep=nstep+1

				if (nstep_count(i) == nstep) then		
					write(0,'(" ====> Time step = ",i8)') nstep
					i=i+1
				end if	

				time=time+dt
				if (nstep == nstep1+maxstp .and. mod(nstep,nprint) /= 0) iend=1

!		Subroutine listep calculates perturvation variable time advance
!		Subroutine energy calculates the radial and poloidal component of the magentic and velocity fields
!		as well as kinetic and magnetic energy of the system.

				if (mod(nstep,nprint) == 0 .or. iend /= 0) call energy(1)
				call linstep
				if (mod(nstep,nprint) == 0 .or. iend /= 0) call energy(2)

				call trnsfr0(psi,1)
				call trnsfr0(phi,-1)
				call trnsfr0(nf,1)
				call trnsfr0(vprlf,-1)

				if (mod(nstep,ndump) == 0 .and. nstep /= nstep1+maxstp) then
					call trnsfr0(pr,1)
					call trnsfr0(vthprlf,-1)
					if (alpha_on == 1) then
						call trnsfr0(nalp,1)
						call trnsfr0(vprlalp,-1)
					end if
!		Subroutine wrdump writes an output file to continue the run	
					if (myPE == 0) then
						call wrdump
						call numinc
					end if
				end if
!               Write high-frequency (every time step) data:
			      if(nonlin .ne. 0) then
			        if(nstep .gt. nstep1+1 .and. myPE == 0)  &
			         open(unit=82,file=hifrq_out,position='append',status='unknown')
			        call hifreq
			        if (myPE == 0) close(unit=82)
			      endif

			end do

!               stepping done.
!		Subroutine lincheck calculates the growth rate and the frequency of the instability	

			if (nonlin == 0) call lincheck

			call trnsfr0(psi,1)
			call trnsfr0(phi,-1)
			call trnsfr0(pr,1)
			call trnsfr0(nf,1)
			call trnsfr0(vprlf,-1)
			call trnsfr0(vthprlf,-1)
			if (alpha_on == 1) then
				call trnsfr0(nalp,1)
				call trnsfr0(vprlalp,-1)
			end if

			numrun(3)="z"
			if (myPE == 0) call wrdump
			call endrun

		end if

		write(0,'(" ====> Simulation DONE !! ")')

		if (myPE == 0) then
			call date_and_time(datew,timew,zone_h,values_e)
			write(6,'(/" time = ",a2,":",a2,":",a2,"       date = ",a2,"/",a2,"/",a4/)') timew(1:2),timew(3:4),timew(5:6), &
												     datew(7:8),datew(5:6),datew(1:4)

			call elapsed_time(values_s,values_e)
		end if

		CALL MPI_FINALIZE(ierror)

	end program far