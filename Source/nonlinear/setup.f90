	subroutine inital

		use param
		use processor
		use var_para
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none
		
		integer :: i,j,l
		real(IDP) :: qlambw,coeff

		interface
			subroutine setmod
			end subroutine setmod
			subroutine mmlims
			end subroutine mmlims
			subroutine grid
			end subroutine grid
			subroutine seteq
			end subroutine seteq
			subroutine etachi
			end subroutine etachi
			subroutine pert
			end subroutine pert
		end interface

!		This subroutine is called if the simulation is new, not a continuation	
		numruns=numrun

		rewind(5)
		if (rc == 0.) rc=rs(1)
		if (delta >= rc .or. delta >= (1.-rc)) delta=min(rc,1.-rc)

!		Subroutine setmod set up the modes of the model
		call setmod

!		n value in the complex exponential representation.
		call mmlims

!		Subroutine grid set up the time independent geometric scale factors
		call grid

!		Subroutine seteq set up	the equilibria
		call seteq

!		Subroutine etachi set up the magnetic diffusivity of the model	
		call etachi

!		Subroutine pert set up the perturbation in the equilibria
		call pert

	end subroutine inital

	subroutine resume

		use param
		use processor
		use var_para
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

		integer :: i,j,l
		character(len=5) :: numeqc
		real(IDP) :: dum

		interface
			subroutine inputlist_s
			end subroutine inputlist_s
			subroutine rddump
			end subroutine rddump
			subroutine setmod
			end subroutine setmod
			subroutine mmlims
			end subroutine mmlims
			subroutine etachi
			end subroutine etachi
			subroutine pert
			end subroutine pert
		end interface

!		This subroutine is called if the simulation is a continuation run		

!		Subroutine rddump reads the data required to continue the run		
		call rddump
                call inputlist_s

!		Subroutine setmod set up the modes of the model  		
		call setmod
		
!		n value in the complex exponential representation.	 			
		call mmlims

!		Change the thermal beta
		bet0=bet0*betath_factor

!		Displacment included to the safety factor / iota
		qq=(1.0_IDP/qqinv) + deltaq
		qqinv = (1.0_IDP/qq) + deltaiota	

!		Subroutine ae_profiles adds the external profiles in the model			
                if(ext_prof .eq. 1) then	
                         call ae_profiles
		end if
		
!		Subroutine etachi set up the magnetic diffusivity of the model			
		call etachi
		
!		Subroutine pert set up the perturbation in the equilibria			
		call pert

	end subroutine resume