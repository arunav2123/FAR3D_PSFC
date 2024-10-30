	subroutine mmlims

	!   Calculate the minimum and maximum m values for each n value
	!   in the complex exponential representation.

		use param
		use domain
		implicit none

		integer :: n,l

		allocate (mmstart(-nmax:nmax))
		allocate (mmend(-nmax:nmax))
		allocate (mmstreq(-nmaxeq:nmaxeq))
		allocate (mmendeq(-nmaxeq:nmaxeq))

		do n = 0,nmax
			mmstart(n) = 1000
			mmend(n) = -1000
		end do

		do l = 1,lmax

			n = nn(l)
			if (n >= 0) then
				if (mm(l) < mmstart(n)) mmstart(n) = mm(l)
				if (mm(l) > mmend(n)) mmend(n) = mm(l)
			else
				if (-mm(l) < mmstart(-n)) mmstart(-n) = -mm(l)
				if (-mm(l) > mmend(-n)) mmend(-n) = -mm(l)
			end if

		end do

		!   in the complex exponential format, [fr,fi](-m,-n) = [fr,-fi](m,n),
		!   so limits for negative n rows are just negatives of corresponding
		!   positive n limits.

		do n = 0,nmax

			mmstart(-n) = -mmend(n)
			mmend(-n) = -mmstart(n)

		end do

		do n = 0,nmaxeq
			mmstreq(n) = 1000
			mmendeq(n) = -1000
		end do

		do l = 1,leqmax

			n = nneq(l)
			if (n >= 0) then
				if (mmeq(l) < mmstreq(n)) mmstreq(n) = mmeq(l)
				if (mmeq(l) > mmendeq(n)) mmendeq(n) = mmeq(l)
			else
				if (-mmeq(l) < mmstreq(-n)) mmstreq(-n) = -mmeq(l)
				if (-mmeq(l) > mmendeq(-n)) mmendeq(-n) = -mmeq(l)
			end if

		end do

		!   in the complex exponential format, [fr,fi](-m,-n) = [fr,-fi](m,n),
		!   so limits for negative n rows are just negatives of corresponding
		!   positive n limits.

		do n = 0,nmaxeq

			mmstreq(-n) = -mmendeq(n)
			mmendeq(-n) = -mmstreq(n)

		end do

	end subroutine mmlims