	subroutine inputlist

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
			 
		character*1 :: cdum0
		character(len=132) :: char5	
		character(len=1) :: t	
		integer :: i
		character(len=16) :: confil	
		real(IDP) :: widthix,dum
		integer :: idel,i_mn,idelp1,ip,mjst,mjnd,iidum	

!		Read the input files
	
		open (unit=5,file="Input_Model",status="old",form="formatted")
		
		read(5,'(a1)') cdum0										 
		read(5,'(a1)') cdum0
		read(5,*) namelist_on
		read(5,'(a1)') cdum0
		read(5,*) nstres
		read(5,'(a1)') cdum0
		read(5,'(2a2)') (numrun(i),i=1,2)
		read(5,'(a1)') cdum0
		read(5,'(2a2,a1)') (numruno(i),i=1,3)
		read(5,'(a1)') cdum0
		read(5,'(a5)') numvac	
                read(5,'(a1)') cdum0
                read(5,*) nonlin	
                read(5,'(a1)') cdum0
                read(5,*) ngeneq
                read(5,'(a1)') cdum0
                read(5,'(a40)') eq_name			
                read(5,'(a1)') cdum0					
                read(5,*) maxstp	
                read(5,'(a1)') cdum0
                read(5,*) dt0	
                read(5,'(a1)') cdum0		
                read(5,*) ldim		
                read(5,'(a1)') cdum0
                read(5,*) leqdim	
                read(5,'(a1)') cdum0	
                read(5,*) jdim	
		read(5,'(a1)') cdum0	
		read(5,*) ext_prof	
		read(5,'(a1)') cdum0	
		read(5,'(a40)') ext_prof_name	
		read(5,'(a1)') cdum0	
		read(5,*) ext_prof_len	
           
                rewind(5)
                close (5)
 
		 allocate (sgnleq(leqdim))		!  Vector with the sign of the equilibrium modes
		 allocate (mmeq(leqdim))		!  Vector with the equilibrium polidal mode numbers
		 allocate (nneq(leqdim))		!  Vector with the equilibrium toroidal mode numbers
		 allocate (mheq(leqdim))		!  Vector with the helicities equilibrium polidal mode numbers
		 allocate (nheq(leqdim))		!  Vector with the helicities equilibrium toroidal mode numbers	
		 allocate (mm(ldim))			!  Vector with the dynamic polidal mode numbers
		 allocate (nn(ldim))			!  Vector with the dynamic toroidal mode numbers
		 allocate (mh(ldim))			!  Vector with the helicities dynamic polidal mode numbers
		 allocate (nh(ldim))			!  Vector with the helicities dynamic toroidal mode numbers	
		 allocate (signl(ldim))			!  Vector with the sign of the dynamic modes
		 allocate (jsl(ldim))	
		 allocate (rs(ldim))		
		 allocate (widthi(ldim))		!  Perturbation width
		 allocate (gammai(ldim))		!  Perturbation parameter
		 allocate (dc1m(jdim))	                !  Derivate weigths to be used in grid subroutine
		 allocate (dc1p(jdim))	
		 allocate (dc2m(jdim))	
		 allocate (dc2p(jdim))	
		 allocate (del2cm(jdim))	
		 allocate (del2cp(jdim))		
		 allocate (r(0:jdim))			!  Normalized minor radius
		 allocate (rinv(0:jdim))		!  Inverse of the normalized minor radius
		 allocate (wt1m(jdim,2))		!  Derivative weights to be used in blockj, block0, b2lx and b2lx0 subroutines
		 allocate (wt10(jdim,2))		
		 allocate (wt1p(jdim,2))		
		 allocate (wt2m(jdim,2))		
		 allocate (wt20(jdim,2))		
		 allocate (wt2p(jdim,2))	
		 allocate (eta(0:jdim))			!  Normalized resistivity (to the value in the axis)
		 allocate (etann(0:jdim))		!  resistivity (m^3 Kg / (s C^2))
		 allocate (sd1(0:jdim))			!  Dummy vectors
		 allocate (sd2(0:jdim))		
		 allocate (sd3(0:jdim))		
		 allocate (sd4(0:jdim))		
		 allocate (sd5(0:jdim))		
		 allocate (sd6(0:jdim))	
		 allocate (sd7(0:jdim))		
		
		 allocate (qq(0:jdim))			!  Safety factor
		 allocate (qqinv(0:jdim))		!  Inverse of the safety factor
		 allocate (qqinvp(0:jdim))		!  Radial derivate of the safety factor inverse
		 allocate (denseq(0:jdim))		!  Normalized thermal plasma density
		 allocate (denseqr(0:jdim))		!  Radial derivate of the normalized thermal plasma density
		 allocate (psieq(0:jdim))		!  Normalized equilibrium poloidal flux
		 allocate (chieq(0:jdim))		!  Normalized equilibrium toroidal flux
		 allocate (preq(0:jdim))		!  Normalized equilibrium thermal plasma pressure
		 allocate (feq(0:jdim))			!  Normalized equilibrium current density
		 allocate (cureq(0:jdim))		!  Normalized equilibrium electric currect
		 allocate (teeq(0:jdim))		!  Normalized equilibrium thermal electron temperature
		 allocate (tieq(0:jdim))		!  Normalized equilibrium thermal ion temperature
		 allocate (nfeq(0:jdim))		!  Normalized equilibrium fast particle density
		 allocate (dnfeqdr(0:jdim))		!  Radial derivate of the normalized fast particle density
		 allocate (vfova(0:jdim))		!  Fast particle thermal velocity normalized to the Alfven velocity
		 allocate (vfova2(0:jdim))		!  vfova^2
		 allocate (vtherm_ion(0:jdim))	        !  User defined normalized thermal ion velocity
		 allocate (vzt_eq(0:jdim))	        !  User defined equilibrium thermal toroidal velocity
		 allocate (vth_eq(0:jdim))	        !  User defined equilibrium thermal poloidal velocity
		 allocate (vth_eq1(0:jdim))	        !  User defined equilibrium thermal poloidal velocity
		 allocate (nalpeq(0:jdim))		!  Normalized equilibrium density of the second fast particle species	
		 allocate (valphaova(0:jdim))      	!  Fast particle thermal velocity normalized to the Alfven velocity of the second fast particle species
		 allocate (valphaova2(0:jdim))	        !  vfalphaova^2
		 allocate (dnalpeqdr(0:jdim))	        !  Radial derivate of the normalized density of the second fast particle species
		 allocate (vtherm_elc(0:jdim))

		 allocate (mj_br(0:numPEsm1))
		 allocate (mj_inc(0:numPEsm1))
		 allocate (mj_st(0:numPEsm1))
		 allocate (mj_dl(0:numPEsm1))

                call dfault
	
		if (namelist_on == 0) then	

                 open (unit=5,file="Input_Model",status="old",form="formatted")

		 read(5,'(a1)') cdum0										 
		 read(5,'(a1)') cdum0
		 read(5,*) namelist_on
		 read(5,'(a1)') cdum0
		 read(5,*) nstres
		 read(5,'(a1)') cdum0
		 read(5,'(2a2)') (numrun(i),i=1,2)
		 read(5,'(a1)') cdum0
		 read(5,'(2a2,a1)') (numruno(i),i=1,3)
		 rewind(5)

		 if (myPE == 0) then
                   confil="farprt"//numrun(1)//numrun(2)
		   open (unit=6,file=confil)
		
		    write(6,'("===================================================================================================")')			
		    write(6,'("========================================================______========_============================")')		
		    write(6,'("============================_____====______====_____===|_____ |=======||===========================")')	
		    write(6,'("===========================| ____|==| ____ |==| ___ |========||=======||===========================")')	
		    write(6,'("===========================||____===||====||==||===||===_____||=======||===========================")')	
		    write(6,'("===========================| ____|==||____||==||___||==|_____ |== ____||===========================")')			
		    write(6,'("===========================||=======| ____ |==||=\\==========||==| ___ |===========================")')	
		    write(6,'("===========================||=======||====||==||==\\====_____||==||___||===========================")')	
		    write(6,'("===========================||=======||====||==||== \\==|______|==|_____|=======VER2.0==============")')	
		    write(6,'("===================================================================================================")')			

		    write(6,'(" copy of Input_Model:")')
		    
		    do
			read(5,'(a)',end=20) char5
			write(6,*) char5
		    end do
20    rewind(5)	
		 end if

	         close(5)
                  
		 if (nstres /= 0) then

			confil="fs"//numruno(1)//numruno(2)//numruno(3)
			open(unit=8,file=confil,status='unknown',form='unformatted')

			read(8) ihist
			rewind(8)

		 end if	

		 ihist=ihist+1
		 allocate (numhist(ihist))

		 if (nstres /= 0) then

			read(8) iidum,numrunp,numrunp,numrunp,iidum,numvac,(numhist(i),i=1,ihist-1)
			rewind(8)

		 end if

		 open (unit=5,file="Input_Model",status="old",form="formatted")		

		 read(5,'(a1)') cdum0										 
		 read(5,'(a1)') cdum0
		 read(5,*) namelist_on
		 read(5,'(a1)') cdum0
		 read(5,*) nstres
		 read(5,'(a1)') cdum0
		 read(5,'(2a2)') (numrun(i),i=1,2)
		 read(5,'(a1)') cdum0
		 read(5,'(2a2,a1)') (numruno(i),i=1,3)
		 read(5,'(a1)') cdum0
		 read(5,'(a5)') numvac	
		 read(5,'(a1)') cdum0
		 read(5,*) nonlin	
		 read(5,'(a1)') cdum0
		 read(5,*) ngeneq	
		 read(5,'(a1)') cdum0
		 read(5,'(a40)') eq_name			
		 read(5,'(a1)') cdum0	
		 read(5,*) maxstp	
		 read(5,'(a1)') cdum0
		 read(5,*) dt0	
		 read(5,'(a1)') cdum0		
		 read(5,*) ldim		
		 read(5,'(a1)') cdum0
		 read(5,*) leqdim	
		 read(5,'(a1)') cdum0	
		 read(5,*) jdim	
		 read(5,'(a1)') cdum0	
		 read(5,*) ext_prof	
		 read(5,'(a1)') cdum0	
		 read(5,'(a40)') ext_prof_name	
		 read(5,'(a1)') cdum0	
		 read(5,*) ext_prof_len		
		 read(5,'(a1)') cdum0	
		 read(5,*) iflr_on	
		 read(5,'(a1)') cdum0		
		 read(5,*) epflr_on	
		 read(5,'(a1)') cdum0	
		 read(5,*) ieldamp_on	
		 read(5,'(a1)') cdum0	
		 read(5,*) twofl_on	
		 read(5,'(a1)') cdum0
		 read(5,*) alpha_on	
		 read(5,'(a1)') cdum0	
		 read(5,*) Trapped_on	
		 read(5,'(a1)') cdum0	
		 read(5,*) matrix_out
		 read(5,'(a1)') cdum0
		 read(5,*) m0dy
		 read(5,'(a1)') cdum0
		 read(5,*) nopsievol_on
		 read(5,'(a1)') cdum0
		 read(5,*) noprevol_on
		 read(5,'(a1)') cdum0
		 read(5,*) nonfevol_on
		 read(5,'(a1)') cdum0
		 read(5,*) nonalpevol_on
		 read(5,'(a1)') cdum0	
					
		 read(5,'(a1)') cdum0										 
		 read(5,'(a1)') cdum0
		 read(5,'(a1)') cdum0		
		 read(5,*) (mm(i),i=1,ldim)
		 read(5,'(a1)') cdum0		
		 read(5,*) (nn(i),i=1,ldim)				
		 read(5,'(a1)') cdum0		
		 read(5,*) (mmeq(i),i=1,leqdim)				
		 read(5,'(a1)') cdum0		
		 read(5,*) (nneq(i),i=1,leqdim)	  
		 read(5,'(a1)') cdum0	
		 read(5,'(a1)') cdum0		
		 read(5,*) ipert	
		 read(5,'(a1)') cdum0		
		 read(5,*) widthix
		 read(5,'(a1)') cdum0		
		 read(5,*) Auto_grid_on		
		 read(5,'(a1)') cdum0		
		 read(5,*) ni
		 read(5,'(a1)') cdum0		
		 read(5,*) nis
		 read(5,'(a1)') cdum0		
		 read(5,*) ne
		 read(5,'(a1)') cdum0		
		 read(5,*) delta
		 read(5,'(a1)') cdum0		
		 read(5,*) rc	
		 read(5,'(a1)') cdum0
		 read(5,*) Edge_on	
		 read(5,'(a1)') cdum0
		 read(5,*) edge_p	
		 read(5,'(a1)') cdum0			
		 read(5,'(a1)') cdum0		
		 read(5,*) gamma			
		 read(5,'(a1)') cdum0		
		 read(5,*) s	
		 read(5,'(a1)') cdum0		
		 read(5,*) betath_factor
		 read(5,'(a1)') cdum0		
		 read(5,*) ietaeq
		 read(5,'(a1)') cdum0		
		 read(5,*) spe1		
		 read(5,'(a1)') cdum0		
		 read(5,*) bet0_f
		 read(5,'(a1)') cdum0		
		 read(5,*) spe2		
		 read(5,'(a1)') cdum0		
		 read(5,*) bet0_alp		
		 read(5,'(a1)') cdum0		
		 read(5,*) omcy			
		 read(5,'(a1)') cdum0		
	 	 read(5,*) omcyb			
		 read(5,'(a1)') cdum0
		 read(5,*) rbound			
		 read(5,'(a1)') cdum0		
		 read(5,*) omcyalp		
		 read(5,'(a1)') cdum0		
		 read(5,*) itime	
		 read(5,'(a1)') cdum0		
		 read(5,*) dpres		
		 read(5,'(a1)') cdum0
		 read(5,'(a1)') cdum0		
		 read(5,*) stdifp			
		 read(5,'(a1)') cdum0		
		 read(5,*) stdifu	
		 read(5,'(a1)') cdum0		
		 read(5,*) stdifv	
		 read(5,'(a1)') cdum0	
		 read(5,*) stdifnf		
		 read(5,'(a1)') cdum0		
		 read(5,*) stdifvf	
		 read(5,'(a1)') cdum0		
		 read(5,*) stdifnalp	
		 read(5,'(a1)') cdum0		
		 read(5,*) stdifvalp
		 read(5,'(a1)') cdum0
		 read(5,'(a1)') cdum0		
		 read(5,*) LcA0			
		 read(5,'(a1)') cdum0		
		 read(5,*) LcA1	
		 read(5,'(a1)') cdum0		
		 read(5,*) LcA2		
		 read(5,'(a1)') cdum0		
		 read(5,*) LcA3	
		 read(5,'(a1)') cdum0		
		 read(5,*) LcA0alp	
		 read(5,'(a1)') cdum0		
		 read(5,*) LcA1alp
		 read(5,'(a1)') cdum0		
		 read(5,*) LcA2alp	
		 read(5,'(a1)') cdum0		
		 read(5,*) LcA3alp		
		 read(5,'(a1)') cdum0
		 read(5,'(a1)') cdum0		
		 read(5,*) omegar			
		 read(5,'(a1)') cdum0		
		 read(5,*) iflr	
		 read(5,'(a1)') cdum0		
		 read(5,*) r_epflr		
		 read(5,'(a1)') cdum0		
		 read(5,*) r_epflralp	
		 read(5,'(a1)') cdum0
		 read(5,'(a1)') cdum0
		 read(5,*) lplots			
		 read(5,'(a1)') cdum0		
		 read(5,*) nprint	
		 read(5,'(a1)') cdum0		
		 read(5,*) ndump	
		 read(5,'(a1)') cdum0			
		 read(5,'(a1)') cdum0
		 read(5,*) DIIID_u	
		
		 read(5,'(a1)') cdum0
		 read(5,'(a1)') cdum0
		 read(5,'(a1)') cdum0		
		 read(5,*) EP_dens_on			
		 read(5,'(a1)') cdum0		
		 read(5,*) Adens			
		 read(5,'(a1)') cdum0		
		 read(5,*) Bdens	
		 read(5,'(a1)') cdum0		
		 read(5,*) Alpha_dens_on		
		 read(5,'(a1)') cdum0		
		 read(5,*) Adensalp		
		 read(5,'(a1)') cdum0		
		 read(5,*) Bdensalp		
		 read(5,'(a1)') cdum0		
		 read(5,*) EP_vel_on	
		 read(5,'(a1)') cdum0	
		 read(5,*) Alpha_vel_on	
		 read(5,'(a1)') cdum0		
		 read(5,*) q_prof_on
		 read(5,'(a1)') cdum0		
		 read(5,*) Eq_vel_on
		 read(5,'(a1)') cdum0	
		 read(5,*) Eq_velp_on				
		 read(5,'(a1)') cdum0	
		 read(5,*) Eq_Presseq_on			
		 read(5,'(a1)') cdum0
		 read(5,*) Eq_Presstot_on			
		 read(5,'(a1)') cdum0
		 read(5,*) deltaq			
		 read(5,'(a1)') cdum0	
		 read(5,*) deltaiota	
		 read(5,'(a1)') cdum0	
		 read(5,*) etascl
		 read(5,'(a1)') cdum0	
		 read(5,*) eta0
		 read(5,'(a1)') cdum0	
		 read(5,*) reta
		 read(5,'(a1)') cdum0	
		 read(5,*) etalmb		
		 read(5,'(a1)') cdum0			
		 read(5,*) (cnep(i),i=0,10)		
		 read(5,'(a1)') cdum0		
		 read(5,*) (ctep(i),i=0,10)	
		 read(5,'(a1)') cdum0	
		 read(5,*) (cnfp(i),i=0,10)			
		 read(5,'(a1)') cdum0		
		 read(5,*) (cvep(i),i=0,10)
		 read(5,'(a1)') cdum0		
		 read(5,*) (cvfp(i),i=0,10)	
		 read(5,'(a1)') cdum0		
		 read(5,*) (cnfpalp(i),i=0,10)	
		 read(5,'(a1)') cdum0		
		 read(5,*) (cvfpalp(i),i=0,10)	
		 read(5,'(a1)') cdum0		
		 read(5,*) (eqvt(i),i=0,10)
		 read(5,'(a1)') cdum0		
		 read(5,*) (eqvp(i),i=0,10)	

		 widthi(:)=widthix

		 if (Auto_grid_on == 1) then

		  ni =  (jdim/2) - 1
		  nis = (jdim/4) + 1
		  ne =  jdim/4
                  
		 end if	

		 rewind(5)
		 close(5)
		
                else

		 open (unit=11,file="farin",status="old")
		
		 read(11,'(i1,2x,2a2,2x,2a2,a1)') nstres,(numrun(i),i=1,2),(numruno(i),i=1,3)
		 rewind(11)
                 
		 t=char(9)

		 if (myPE == 0) then

			confil="farprt"//numrun(1)//numrun(2)
			open (unit=6,file=confil)

		        write(6,'("===================================================================================================")')			
		        write(6,'("========================================================______========_============================")')		
		        write(6,'("============================_____====______====_____===|_____ |=======||===========================")')	
		        write(6,'("===========================| ____|==| ____ |==| ___ |========||=======||===========================")')	
		        write(6,'("===========================||____===||====||==||===||===_____||=======||===========================")')	
		        write(6,'("===========================| ____|==||____||==||___||==|_____ |== ____||===========================")')			
		        write(6,'("===========================||=======| ____ |==||=\\==========||==| ___ |===========================")')	
		        write(6,'("===========================||=======||====||==||==\\====_____||==||___||===========================")')	
		        write(6,'("===========================||=======||====||==||== \\==|______|==|_____|=======VER2.0==============")')	
		        write(6,'("===================================================================================================")')	

			write(6,'(" copy of in:")')

			do while (.true.)
				read(11,'(a)',end=40) char5
				write(6,*) char5
			end do
40    			rewind(11)

		 end if

		 call dfault

		 read(11,'(i1,2x,2a2,2x,2a2,a1)') nstres,(numrun(i),i=1,2),(numruno(i),i=1,3)

		 if (nstres /= 0) then

			confil="fs"//numruno(1)//numruno(2)//numruno(3)
			open(unit=8,file=confil,status='unknown',form='unformatted')

			read(8) ihist
			rewind(8)

		 end if

		 ihist=ihist+1
		 allocate (numhist(ihist))

		 if (nstres /= 0) then

			read(8) idum,numrunp,numrunp,numrunp,idum,numvac,(numhist(i),i=1,ihist-1)
			rewind(8)

		 end if

		 read(11,nam_par)
		 read(11,nam_arr)
                 rewind(11)
                 close(11)
               end if

		 mjm1=mj-1
		 idel=mjm1/numPEs
		 idelp1=idel+1
		 i_mn=numPEs*idelp1-mjm1
		 if (i_mn > 1) i_mn=i_mn-1
		 do ip=0,numPEsm1
			if (ip < i_mn) then
				mjst=ip*idel
				mjnd=(ip+1)*idel+1
			else
				mjst=ip*idelp1-i_mn
				mjnd=(ip+1)*idelp1-i_mn+1
			end if
			if (ip == numPEsm1) mjnd=mj
			mj_br(ip)=mjst+1
			if (ip == 0) mj_br(ip)=mjst
			mj_inc(ip)=mjnd-mj_br(ip)
			if (ip == numPEsm1) mj_inc(ip)=mj_inc(ip)+1
			mj_st(ip)=mjst
			mj_dl(ip)=mjnd-mjst+1
			if (ip == myPE) mj_start=mjst
			if (ip == myPE) mj_end=mjnd
		 end do

                  if (myPE == 0) then		
		   allocate (psi(0:jdim,0:ldim))	!  Normalized fluctuating poloidal flux
		   allocate (phi(0:jdim,0:ldim))	!  Normalized fluctuating electrostatic potential
		   allocate (pr(0:jdim,0:ldim))		!  Normalized fluctuating pressure
		   allocate (uzt(0:jdim,0:ldim))	!  Normalized fluctuating vorticity (toroidal component)
		   allocate (nf(0:jdim,0:ldim))		!  Normalized fluctuating fast particle density
		   allocate (vprlf(0:jdim,0:ldim))	!  Normalized fluctuating fast particle parallel velocity
		   allocate (vthprlf(0:jdim,0:ldim))	!  Normalized fluctuating parallel thermal velocity
	           allocate (nalp(0:jdim,0:ldim))	!  Normalized fluctuating density of the second fast particle species
		   allocate (vprlalp(0:jdim,0:ldim))	!  Normalized fluctuating parallel velocity of the second fast particle species
		 else
		   allocate (psi(mj_start:mj_end,0:ldim))	!  Normalized fluctuating poloidal flux
		   allocate (phi(mj_start:mj_end,0:ldim))	!  Normalized fluctuating electrostatic potential
		   allocate (pr(mj_start:mj_end,0:ldim))	!  Normalized fluctuating pressure
		   allocate (uzt(mj_start:mj_end,0:ldim))	!  Normalized fluctuating vorticity (toroidal component)
		   allocate (nf(mj_start:mj_end,0:ldim))	!  Normalized fluctuating fast particle density
		   allocate (vprlf(mj_start:mj_end,0:ldim))	!  Normalized fluctuating fast particle parallel velocity
		   allocate (vthprlf(mj_start:mj_end,0:ldim))	!  Normalized fluctuating parallel thermal velocity
	           allocate (nalp(mj_start:mj_end,0:ldim))	!  Normalized fluctuating density of the second fast particle species
		   allocate (vprlalp(mj_start:mj_end,0:ldim))	!  Normalized fluctuating parallel velocity of the second fast particle species
                 end if

		 allocate (psip(mj_start:mj_end,0:ldim))	        !  Normalized fluctuating poloidal flux
		 allocate (prp(mj_start:mj_end,0:ldim))	        !  Normalized fluctuating pressure
		 allocate (uztp(mj_start:mj_end,0:ldim))	        !  Normalized fluctuating vorticity (toroidal component)
		 allocate (nfpp(mj_start:mj_end,0:ldim))	        !  Normalized fluctuating fast particle density
		 allocate (vprlfp(mj_start:mj_end,0:ldim))	!  Normalized fluctuating fast particle parallel velocity
		 allocate (vthprlfp(mj_start:mj_end,0:ldim))	!  Normalized fluctuating parallel thermal velocity
	         allocate (nalpp(mj_start:mj_end,0:ldim) )	!  Normalized fluctuating density of the second fast particle species
		 allocate (vprlalpp(mj_start:mj_end,0:ldim))	!  Normalized fluctuating parallel velocity of the second fast particle species

		 allocate (psi_nl(mj_start:mj_end,0:ldim,3))	    !  Normalized fluctuating poloidal flux
		 allocate (pr_nl(mj_start:mj_end,0:ldim,3))	    !  Normalized fluctuating pressure
		 allocate (uzt_nl(mj_start:mj_end,0:ldim,3))	    !  Normalized fluctuating vorticity (toroidal component)
		 allocate (nf_nl(mj_start:mj_end,0:ldim,3))	    !  Normalized fluctuating fast particle density
		 allocate (vprlf_nl(mj_start:mj_end,0:ldim,3))	    !  Normalized fluctuating fast particle parallel velocity
		 allocate (vthprlf_nl(mj_start:mj_end,0:ldim,3))      !  Normalized fluctuating parallel thermal velocity
	         allocate (nalp_nl(mj_start:mj_end,0:ldim,3))	    !  Normalized fluctuating density of the second fast particle species
		 allocate (vprlalp_nl(mj_start:mj_end,0:ldim,3))      !  Normalized fluctuating parallel velocity of the second fast particle species
		
		 allocate (sc1(mj_start:mj_end,0:ldim))		!  dummy variable for dynamic modes
		 allocate (sc2(mj_start:mj_end,0:ldim))	
		 allocate (sc3(mj_start:mj_end,0:ldim))	
		 allocate (sc4(mj_start:mj_end,0:ldim))	
		 allocate (sc5(mj_start:mj_end,0:ldim))	
		 allocate (sc6(mj_start:mj_end,0:ldim))	
		 allocate (sc7(mj_start:mj_end,0:ldim))	
		 allocate (sc8(mj_start:mj_end,0:ldim))		
		 allocate (sc9(mj_start:mj_end,0:ldim))	
		 allocate (sc10(mj_start:mj_end,0:ldim))	
		 allocate (sc11(mj_start:mj_end,0:ldim))
		 allocate (sc12(mj_start:mj_end,0:ldim))
											
		 allocate (sceq1(0:jdim,0:leqdim))	! dummy variable for equilibrium modes
		 allocate (sceq2(0:jdim,0:leqdim))	
		 allocate (sceq3(0:jdim,0:leqdim))	
		 allocate (sceq4(0:jdim,0:leqdim))	
		 allocate (sceq5(0:jdim,0:leqdim))	
		 allocate (sceq6(0:jdim,0:leqdim))	
		 allocate (sceq7(0:jdim,0:leqdim))		

		 allocate (epsi(ldim,2),ephi(ldim,2),epr(ldim,2),eprnc(ldim,2),ekenc(ldim,2),eke(ldim,2), &
			  emenc(ldim,2),eme(ldim,2),ealp(ldim,2),ealpnc(ldim,2))

	end subroutine inputlist
