	subroutine energy(i)

		use param
		use processor
		use var_para
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		integer, save :: iflg = 0
		integer :: i,j,l,lp,l1,l2,mjp,mjnd
		real(IDP) :: denom,scnorm,emetot,eketot,eprtot,ealptot,emetot0,eketot0,eprtot0,ealptot0
		character(len=1) :: t
		character(len=16) :: confil
		character(len=32) :: format1='("time",1200(a1,i4,"/",i4))'
		character(len=60) :: formatt='("time",a1,"total",a1,"total (w/o 0/0)",1200(a1,i4,"/",i4))'
		character(len=32) :: formatv='(1pe13.6,1202(a1,1pe15.8))'
		integer, dimension(ldim) :: lmap
		real(IDP), dimension(ldim) :: gampsi,gamalpha,gamphi,gampr
		real(IDP), dimension(jdim) :: rint,xint,bb,cc,dd
		real(IDP), dimension(0:jdim) :: xinte,bbe,cce,dde
		
		interface
			subroutine dbydth_par(d,a,ltype,c1,c2,k)
				use param
				use var_para
				implicit none
				integer :: k,ltype
				real(IDP) :: c1,c2
				real(IDP), dimension(mj_start:,0:) :: a,d
			end subroutine dbydth_par
			subroutine dbydr_par(d,a,c1,c2,k)
				use param
				use var_para
				implicit none
				integer :: k
				real(IDP) :: c1,c2
				real(IDP), dimension(mj_start:,0:) :: a,d
			end subroutine dbydr_par
			subroutine eqtodyn_par(adyn,aeq,c1,c2)
				use param
				use var_para
				implicit none
				real(IDP) :: c1,c2
				real(IDP), dimension(mj_start:,0:) :: adyn
				real(IDP), dimension(0:,0:) :: aeq
			end subroutine eqtodyn_par
			subroutine mult(f,g,itypeg,h,itypeh,c1,c2)
				use param
				use var_para
				implicit none
				integer :: itypeg,itypeh
				real(IDP) :: c1,c2
				real(IDP), dimension(mj_start:,0:) :: f,g,h
			end subroutine mult
			subroutine trnsfr0l(f)
				use param
				use var_para
				implicit none
				real(IDP), dimension(mj_start:,0:) :: f
			end subroutine trnsfr0l
			subroutine quadq(n,x,y,x0,b,c,d,result)
				use param
				implicit none
				integer :: n
				real(IDP) :: x0,result
				real(IDP), dimension(:) :: x,y,b,c,d
			end subroutine quadq
		end interface


		epsi(:,i)=0.0_IDP 
		ephi(:,i)=0.0_IDP 
		epr(:,i)=0.0_IDP 
		eprnc(:,i)=0.0_IDP 
		ekenc(:,i)=0.0_IDP 
		eke(:,i)=0.0_IDP 
		emenc(:,i)=0.0_IDP 
		eme(:,i)=0.0_IDP 
		ealp(:,i)=0.0_IDP 
		ealpnc(:,i)=0.0_IDP 

!  store uzt-values into sc4

		do l=1,lmax
			sc4(:,l)=uzt(mj_start:mj_end,l)
		end do

!  vr up
		call dbydth_par(sc2,phi,-1,0.0_IDP,-1.0_IDP,0)
!  vth up
		call dbydr_par(sc3,phi,0.0_IDP,1.0_IDP,0)
		
		mjp=mj+1
		do l=1,lmaxn
			uzt(mj_start:mj_end,l)=jbgrr(mj_start:mj_end,leq0)*sc2(:,lln(l))**2+jbgtt(mj_start:mj_end,leq0)*sc3(:,lln(l))**2
		end do
		call trnsfr0l(uzt)
		if (myPE == 0) then
			do l=1,lmaxn
				xinte=r*denseq*uzt(:,l)
				call quadq(mjp,r,xinte,r(mj),bbe,cce,dde,ekenc(l,i))
				if (lln(l) == l0) ekenc(l,i)=2.*ekenc(l,i)
			end do
		end if

		call eqtodyn_par(sc8,jbgrr,0.0_IDP,1.0_IDP)
		call mult(sc5,sc8,1,sc2,1,0.0_IDP,1.0_IDP)
		call eqtodyn_par(sc8,jbgrt,0.0_IDP,1.0_IDP)
		call mult(sc5,sc8,-1,sc3,-1,1.0_IDP,1.0_IDP)
		call mult(sc6,sc8,-1,sc2,1,0.0_IDP,1.0_IDP)
		call eqtodyn_par(sc8,jbgtt,0.0_IDP,1.0_IDP)
		call mult(sc6,sc8,1,sc3,-1,1.0_IDP,1.0_IDP)
		do l=1,lmaxn
			uzt(mj_start:mj_end,l)=sc2(:,lln(l))*sc5(:,lln(l))+sc3(:,lln(l))*sc6(:,lln(l))
		end do
		call trnsfr0l(uzt)
		if (myPE == 0) then
			do l=1,lmaxn
				xinte=r*denseq*uzt(:,l)
				call quadq(mjp,r,xinte,r(mj),bbe,cce,dde,eke(l,i))
				if (lln(l) == l0) eke(l,i)=2.*eke(l,i)
			end do
		end if

!  br up
		call dbydth_par(sc2,psi,1,0.0_IDP,-1.0_IDP,0)
!  bth up
		call dbydr_par(sc3,psi,0.0_IDP,1.0_IDP,0)

		do l=1,lmaxn
			uzt(mj_start:mj_end,l)=grroj(mj_start:mj_end,leq0)*sc2(:,lln(l))**2+gttoj(mj_start:mj_end,leq0)*sc3(:,lln(l))**2
		end do
		call trnsfr0l(uzt)
		if (myPE == 0) then
			do l=1,lmaxn
				xinte=r*uzt(:,l)
				call quadq(mjp,r,xinte,r(mj),bbe,cce,dde,emenc(l,i))
				if (lln(l) == l0) emenc(l,i)=2.*emenc(l,i)
			end do
		end if

		call eqtodyn_par(sc8,grroj,0.0_IDP,1.0_IDP)
		call mult(sc5,sc8,1,sc2,-1,0.0_IDP,1.0_IDP)
		call eqtodyn_par(sc8,grtoj,0.0_IDP,1.0_IDP)
		call mult(sc5,sc8,-1,sc3,1,1.0_IDP,1.0_IDP)
		call mult(sc6,sc8,-1,sc2,-1,0.0_IDP,1.0_IDP)
		call eqtodyn_par(sc8,gttoj,0.0_IDP,1.0_IDP)
		call mult(sc6,sc8,1,sc3,1,1.0_IDP,1.0_IDP)
		do l=1,lmaxn
			uzt(mj_start:mj_end,l)=sc2(:,lln(l))*sc5(:,lln(l))+sc3(:,lln(l))*sc6(:,lln(l))
		end do
		call trnsfr0l(uzt)
		if (myPE == 0) then
			do l=1,lmaxn
				xinte=r*uzt(:,l)
				call quadq(mjp,r,xinte,r(mj),bbe,cce,dde,eme(l,i))
				if (lln(l) == l0) eme(l,i)=2.*eme(l,i)
			end do
		end if

		call eqtodyn_par(sc8,sqg,0.0_IDP,1.0_IDP)
		call mult(sc2,sc8,1,vprlf,-1,0.0_IDP,1.0_IDP)
		do l=1,lmaxn
			uzt(mj_start:mj_end,l)=nfeq(mj_start:mj_end)*sqg(mj_start:mj_end,leq0)*vprlf(mj_start:mj_end,lln(l))**2
		end do
		call trnsfr0l(uzt)
		if (myPE == 0) then
			do l=1,lmaxn
				xinte=r*uzt(:,l)
				call quadq(mjp,r,xinte,r(mj),bbe,cce,dde,eprnc(l,i))
				if (lln(l) == l0) eprnc(l,i)=2.*eprnc(l,i)
			end do
		end if
		do l=1,lmaxn
			uzt(mj_start:mj_end,l)=nfeq(mj_start:mj_end)*sc2(:,lln(l))*vprlf(mj_start:mj_end,lln(l))
		end do
		call trnsfr0l(uzt)
		if (myPE == 0) then
			do l=1,lmaxn
				xinte=r*uzt(:,l)
				call quadq(mjp,r,xinte,r(mj),bbe,cce,dde,epr(l,i))
				if (lln(l) == l0) epr(l,i)=2.*epr(l,i)
			end do
		end if

		if (alpha_on == 1) then
			call mult(sc2,sc8,1,vprlalp,-1,0.0_IDP,1.0_IDP)
			do l=1,lmaxn
				uzt(mj_start:mj_end,l)=nalpeq(mj_start:mj_end)*sqg(mj_start:mj_end,leq0)*vprlalp(mj_start:mj_end,lln(l))**2
			end do
			call trnsfr0l(uzt)
			if (myPE == 0) then
				do l=1,lmaxn
					xinte=r*uzt(:,l)
					call quadq(mjp,r,xinte,r(mj),bbe,cce,dde,ealpnc(l,i))
					if (lln(l) == l0) ealpnc(l,i)=2.*ealpnc(l,i)
				end do
			end if
			do l=1,lmaxn
				uzt(mj_start:mj_end,l)=nalpeq(mj_start:mj_end)*sc2(:,lln(l))*vprlalp(mj_start:mj_end,lln(l))
			end do
			call trnsfr0l(uzt)
			if (myPE == 0) then
				do l=1,lmaxn
					xinte=r*uzt(:,l)
					call quadq(mjp,r,xinte,r(mj),bbe,cce,dde,ealp(l,i))
					if (lln(l) == l0) ealp(l,i)=2.*ealp(l,i)
				end do
			end if
		end if

!  restore uzt-values

		do l=1,lmax
			uzt(mj_start:mj_end,l)=sc4(:,l)
		end do

		if (i == 1) return

		if (myPE == 0) then

			write(6,'(/"energy:numrun=",2a2,a1,",numruno=",2a2,a1,",nstep=",i10,",time=",1pe12.5,",dt=",1pe12.5)') numrun,numruno, &
					nstep,time,dt
			lp=0
			emetot=0.0
			eketot=0.0
			eprtot=0.0
			ealptot=0.0
			emetot0=0.0
			eketot0=0.0
			eprtot0=0.0
			ealptot0=0.0
			do l=1,lmaxn
				l1=lln(l)
				if (signl(l1) < 0) cycle
				if (signl(l1) > 0 .and. lo(l) > 0) then
					l2=lo(l)
					ekenc(l,1)=ekenc(l,1)+ekenc(l2,1)
					eke(l,1)=eke(l,1)+eke(l2,1)
					emenc(l,1)=emenc(l,1)+emenc(l2,1)
					eme(l,1)=eme(l,1)+eme(l2,1)
					eprnc(l,1)=eprnc(l,1)+eprnc(l2,1)
					epr(l,1)=epr(l,1)+epr(l2,1)
					ealpnc(l,1)=ealpnc(l,1)+ealpnc(l2,1)
					ealp(l,1)=ealp(l,1)+ealp(l2,1)
					ekenc(l,2)=ekenc(l,2)+ekenc(l2,2)
					eke(l,2)=eke(l,2)+eke(l2,2)
					emenc(l,2)=emenc(l,2)+emenc(l2,2)
					eme(l,2)=eme(l,2)+eme(l2,2)
					eprnc(l,2)=eprnc(l,2)+eprnc(l2,2)
					epr(l,2)=epr(l,2)+epr(l2,2)
					ealpnc(l,2)=ealpnc(l,2)+ealpnc(l2,2)
					ealp(l,2)=ealp(l,2)+ealp(l2,2)
				end if
				gampsi(l)=0.0_IDP
				denom=dt*(emenc(l,1)+emenc(l,2))
				if (denom /= 0.0_IDP) gampsi(l)=(emenc(l,2)-emenc(l,1))/denom
				gamphi(l)=0.0_IDP
				denom=dt*(ekenc(l,1)+ekenc(l,2))
				if (denom /= 0.0_IDP) gamphi(l)=(ekenc(l,2)-ekenc(l,1))/denom
				gampr(l)=0.0_IDP
				denom=dt*(eprnc(l,1)+eprnc(l,2))
				if (denom /= 0.0_IDP) gampr(l)=(eprnc(l,2)-eprnc(l,1))/denom
				emetot=emetot+eme(l,2)
				eketot=eketot+eke(l,2)
				eprtot=eprtot+epr(l,2)
				ealptot=ealptot+ealp(l,2)
				if (lln(l) /= l0) then
					emetot0=emetot0+eme(l,2)
					eketot0=eketot0+eke(l,2)
					eprtot0=eprtot0+epr(l,2)
					ealptot0=ealptot0+ealp(l,2)
				end if
				lp=lp+1
				lmap(lp)=l
			end do
			t=char(9)
			if (iflg == 0) then
				write(confil,'("eme_",2a2)') numrun(1),numrun(2)
				open(unit=83,file=confil,recl=19384)
				if (nstres_com == 0) write(83,formatt) t,t,(t,mm(lln(lmap(l))),nn(lln(lmap(l))),l=1,lp)
				write(confil,'("eke_",2a2)') numrun(1),numrun(2)
				open(unit=84,file=confil,recl=19384)
				if (nstres_com == 0) write(84,formatt) t,t,(t,mm(lln(lmap(l))),nn(lln(lmap(l))),l=1,lp)
				write(confil,'("evprlf_",2a2)') numrun(1),numrun(2)
				open(unit=85,file=confil,recl=19384)
				if (nstres_com == 0) write(85,formatt) t,t,(t,mm(lln(lmap(l))),nn(lln(lmap(l))),l=1,lp)
				write(confil,'("emenc_",2a2)') numrun(1),numrun(2)
				open(unit=93,file=confil,recl=19384)
				if (nstres_com == 0) write(93,format1) (t,mm(lln(lmap(l))),nn(lln(lmap(l))),l=1,lp)
				write(confil,'("ekenc_",2a2)') numrun(1),numrun(2)
				open(unit=94,file=confil,recl=19384)
				if (nstres_com == 0) write(94,format1) (t,mm(lln(lmap(l))),nn(lln(lmap(l))),l=1,lp)
				write(confil,'("evprlfnc_",2a2)') numrun(1),numrun(2)
				open(unit=95,file=confil,recl=19384)
				if (nstres_com == 0) write(95,format1) (t,mm(lln(lmap(l))),nn(lln(lmap(l))),l=1,lp)
			end if
			write(83,formatv) time,t,eps*eps*emetot,t,eps*eps*emetot0,(t,eps*eps*eme(lmap(l),2),l=1,lp) 
			write(84,formatv) time,t,eps*eps*eketot,t,eps*eps*eketot0,(t,eps*eps*eke(lmap(l),2),l=1,lp) 
			write(85,formatv) time,t,bet0_f*eprtot/2.0,t,bet0_f*eprtot0/2.0,(t,bet0_f*epr(lmap(l),2)/2.0,l=1,lp) 
			write(93,formatv) time,(t,eps*eps*emenc(lmap(l),2),l=1,lp) 
			write(94,formatv) time,(t,eps*eps*ekenc(lmap(l),2),l=1,lp) 
			write(95,formatv) time,(t,bet0_f*eprnc(lmap(l),2)/2.0,l=1,lp) 

			if (alpha_on == 1) then
				do l=1,lmaxn
					gamalpha(l)=0.0_IDP
					denom=dt*(ealpnc(l,1)+ealpnc(l,2))
					if (denom /= 0.0_IDP) gamalpha(l)=(ealpnc(l,2)-ealpnc(l,1))/denom
				end do
				write(6,'(/"   l   m/  n        ke        me     vprlf   vprlalp     gamke     gamme    gamvpr    gamalp")')
				do l=1,lmaxn
					l1=lln(l)
					if (signl(l1) < 0) cycle
					write(6,'(3i4,1p8e10.3)') l1,mm(l1),nn(l1),eps*eps*ekenc(l,2),eps*eps*emenc(l,2),bet0_f*eprnc(l,2)/2.0, &
								  bet0_alp*ealpnc(l,2)/2.0,gamphi(l),gampsi(l),gampr(l),gamalpha(l)
				end do
				if (iflg == 0) then
					write(confil,'("evalp_",2a2)') numrun(1),numrun(2)
					open(unit=86,file=confil,recl=19384)
					if (nstres_com == 0) write(86,formatt) t,t,(t,mm(lln(lmap(l))),nn(lln(lmap(l))),l=1,lp)
					write(confil,'("evalpnc_",2a2)') numrun(1),numrun(2)
					open(unit=96,file=confil,recl=19384)
					if (nstres_com == 0) write(96,format1) (t,mm(lln(lmap(l))),nn(lln(lmap(l))),l=1,lp)
				end if
				write(86,formatv) time,t,bet0_alp*ealptot/2.0,t,bet0_alp*ealptot0/2.0,(t,bet0_alp*ealp(lmap(l),2)/2.0,l=1,lp) 
				write(96,formatv) time,(t,bet0_alp*ealpnc(lmap(l),2)/2.0,l=1,lp) 
			else
				write(6,'(/"   l   m/  n        ke        me     vprlf     gamke     gamme    gamvpr")')
				do l=1,lmaxn
					l1=lln(l)
					if (signl(l1) < 0) cycle
					write(6,'(3i4,1p6e10.3)') l1,mm(l1),nn(l1),eps*eps*ekenc(l,2),eps*eps*emenc(l,2),bet0_f*eprnc(l,2)/2.0, &
								  gamphi(l),gampsi(l),gampr(l)
				end do
			end if

			iflg=1

			if (nonlin == 0 .and. nstep == nstep1+maxstp) then
				scnorm=max(maxval(ekenc),maxval(eke),maxval(emenc),maxval(eme))
				if (scnorm /= 0.0_IDP) then
					ekenc(1:lmaxn,2)=ekenc(1:lmaxn,2)/scnorm
					eke(1:lmaxn,2)=eke(1:lmaxn,2)/scnorm
					emenc(1:lmaxn,2)=emenc(1:lmaxn,2)/scnorm
					eme(1:lmaxn,2)=eme(1:lmaxn,2)/scnorm
					eprnc(1:lmaxn,2)=bet0_f*eprnc(1:lmaxn,2)/(2.0*eps*eps*scnorm)
					epr(1:lmaxn,2)=bet0_f*epr(1:lmaxn,2)/(2.0*eps*eps*scnorm)
					ealpnc(1:lmaxn,2)=bet0_alp*ealpnc(1:lmaxn,2)/(2.0*eps*eps*scnorm)
					ealp(1:lmaxn,2)=bet0_alp*ealp(1:lmaxn,2)/(2.0*eps*eps*scnorm)
				end if
				write(confil,'("spctr_",2a2)') numrun(1),numrun(2)
				open(unit=92,file=confil)
				if (alpha_on == 1) then
					write(92,'("l",a1,"m",a1,"n",a1,"kenc",a1,"ke",a1,"menc",a1,"me",a1,"vprlfnc",a1,"vprlf",a1,"vprlalpnc",a1,"vprlalp")') &
						 (t,l=1,10)
					do l=1,lmaxn
						if (signl(lln(l)) < 0) cycle
						write(92,'(i4,2(a1,i4),8(a1,1pe13.6))') lln(l),t,mm(lln(l)),t,nn(lln(l)),t,ekenc(l,2),t,eke(l,2),t,emenc(l,2),t,eme(l,2), &
							t,eprnc(l,2),t,epr(l,2),t,ealpnc(l,2),t,ealp(l,2)
					end do
				else
					write(92,'("l",a1,"m",a1,"n",a1,"kenc",a1,"ke",a1,"menc",a1,"me",a1,"vprlfnc",a1,"vprlf")') (t,l=1,8)
					do l=1,lmaxn
						if (signl(lln(l)) < 0) cycle
						write(92,'(i4,2(a1,i4),6(a1,1pe13.6))') lln(l),t,mm(lln(l)),t,nn(lln(l)),t,ekenc(l,2),t,eke(l,2),t,emenc(l,2),t,eme(l,2), &
							t,eprnc(l,2),t,epr(l,2)
					end do
				end if
				close(92)
			end if

		end if

	end subroutine energy