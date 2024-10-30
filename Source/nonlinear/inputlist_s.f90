	subroutine inputlist_s

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
		integer :: idel,i_mn,idelp1,ip,mjst,mjnd	

!		Read the input files
		
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
		 read(11,nam_par)
		 read(11,nam_arr)
                 rewind(11)
                 close(11)
               end if

	end subroutine inputlist_s
