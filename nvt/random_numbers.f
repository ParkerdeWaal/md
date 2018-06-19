c-----------------------------------------------------------------------
c     setran calls setrn to initialize the rannum random number
c     generator and calls rannum once.  I don't know why.
c-----------------------------------------------------------------------
c
      subroutine setran
c
c
      include 'common_files.inc'
      common/rnd/pseed
c
      call setrn(pseed)
      pseed = rannum()

      return
      end 
c
c-----------------------------------------------------------------------
c     setrn initializes the rannum random number generator.
c-----------------------------------------------------------------------
c     see comments under rannum()....sjs 7/1/98
c-----------------------------------------------------------------------
c
      subroutine setrn(q)

      include 'common_files.inc'
      COMMON/CRAN/QBASE,QA1,QA2,QB1,QB2

c***  INITIALIZE WITH A CALL TO SETRN(0.D0-1.D0)

c     these don't change.  linear congruential multipliers?

      qa1=2057713.0d0
      qa2=16676923.0d0

c     24-bit integers will be generated (as doubles):

      qbase=2.d0**24

c     generate 48-bit number from seed

      qc=dint(qbase*(qbase*q))

c     first 24 bits worth:

      qb1=dint(qc/qbase)

c     next 24 bits worth:

      qb2=qc-qb1*qbase

c     correct for seed outside [0,1):

      qb1=dmod(qb1,qbase)

c     make qb2 odd:

      qb2=dint(qb2/2.0d0) * 2.0d0 + 1.0d0
c
      return
      end
c
c-----------------------------------------------------------------------
c     rannum generates a random number in [0,1). 
c-----------------------------------------------------------------------
c     it works by using 24-bit
c     integers (stored as doubles) qa1, qa2, qb1, qb2.  the most
c     significant 24 bits of qa2 * qb2 are added to the 24 least
c     significant bits of qa1 * qb2 and the 24 least significant bits
c     of qa2 * qb1 (and modded) to get a 24-bit random integer, which
c     is recycled as qb1 for next time.  qb2 gets the 24 least
c     significant bits of the qa2 * qb2 result.  
c-----------------------------------------------------------------------
c     the properties of this random number generator are better than
c     those of qran(), below, but not as good as gran().  two
c     disadvantages of this random number generator are:
c     1) it is not portable and self-contained, since it relies on
c        common variables to preserve qa1, qa2, qb1, and qb2.
c     2) there is no seed which can be used to restart a series of
c        random numbers at a particular place.  (it can be seeded, via
c        setrn, but there is no way to pick up midstream from an
c        in-progress series)
c-----------------------------------------------------------------------
c
      function rannum()

      include 'common_files.inc'
      COMMON/CRAN/QBASE,QA1,QA2,QB1,QB2

c***  FROM CLAMPS AT NRCC - FROM KALOS

      qd2=qa2*qb2
      qe2=dint(qd2/qbase)
      qc2=qd2-qbase*qe2
      qb1=dmod(qe2+dmod(qa1*qb2,qbase)+dmod(qa2*qb1,qbase),qbase)
      qb2=qc2
      rannum=qb1/qbase

      return
      end
c
c-----------------------------------------------------------------------
c     qran is a (quick) normal, linear congruential random number  
c     generator with "well-chosen" constants taken from Numerical
c     Recipes.  Call it with a positive integer seed.  It'll change that
c     seed and return a number in [0,1)...sjs 3/18/92
c-----------------------------------------------------------------------
c     Note:  Numerical Recipes doesn't SAVE the variables, which it
c     should.
c-----------------------------------------------------------------------
c     This could be better implemented using reals.
c-----------------------------------------------------------------------
c
      function qran(iseed)

      include 'common_files.inc'

      parameter (imod = 714025, imult = 1366, iadd = 150889)

      save

      real*8 qran

      iseed = mod(imult * iseed + iadd, imod)
      qran = dble(iseed) / dble(imod)
      return
      end
c
c-----------------------------------------------------------------------
c     gran is a (good) more reliable random number generator based on 3
c     linear congruential cycles, with shuffling.  It was taken straight
c     from Numerical Recipes....sjs 3/18/92
c-----------------------------------------------------------------------
c     Note:  Numerical Recipes doesn't SAVE the variables, which it 
c     should.
c-----------------------------------------------------------------------
c
      function gran(iseed)

      include 'common_files.inc'

      logical init

      parameter (isize = 97)
      parameter (imod1 = 259200, imult1 = 7141, iadd1 = 54773, 
     .     div1 = 1.d0 / imod1)
      parameter (imod2 = 134456, imult2 = 8121, iadd2 = 28411, 
     .     div2 = 1.d0 / imod2)
      parameter (imod3 = 243000, imult3 = 4561, iadd3 = 51349)

      save

      dimension store(isize)

      data init/.true./

      if (iseed .lt. 0 .or. init) then
         init = .false.
         iseed1 = mod(iadd1 - iseed, imod1)
         iseed1 = mod(imult1 * iseed1 + iadd1, imod1)
         iseed2 = mod(iseed1, imod2)
         iseed2 = mod(imult2 * iseed2 + iadd2, imod2)
         iseed3 = mod(iseed1, imod3)
         do 11 islot = 1, isize
            iseed1 = mod(imult1 * iseed1 + iadd1, imod1)
            iseed2 = mod(imult2 * iseed2 + iadd2, imod2)
            store(islot) = (dble(iseed1) + dble(iseed2) * div2) * div1
 11      continue
c         iseed = 1
      endif
      iseed1 = mod(imult1 * iseed1 + iadd1, imod1)
      iseed2 = mod(imult2 * iseed2 + iadd2, imod2)
      iseed3 = mod(imult3 * iseed3 + iadd3, imod3)
      islot = 1 + (isize * iseed3) / imod3
      gran = store(islot)
      store(islot) = (dble(iseed1) + dble(iseed2) * div2) * div1
      return
      end
c
c-----------------------------------------------------------------------
c     boxmul will return two Gaussian random variables with
c     zero mean and unit variance.  they are calculated from two
c     uniformly distributed random variables using the Box-Muller
c     algorithm:
c       given uniformly distributed and independent r1, r2,
c         x = sqrt(-2 ln(r1)) cos(2 pi r2)
c         y = sqrt(-2 ln(r1)) sin(2 pi r2)
c       are Gaussian and independent
c     ...sjs 3/17/98
c-----------------------------------------------------------------------
c     modified to use the polar form of the Box-Muller algorithm.  
c     see Numerical Recipes, Sec 7.2.... -bd
c-----------------------------------------------------------------------
c     modified to use the irand() random number generator....
c-----------------------------------------------------------------------
c     
      subroutine boxmul(x,y)

      include 'common_files.inc'

      pi2 = acos(-1.d0) * 2.d0

 100  continue
      rand1 = 2.d0 * drand() - 1.d0
      rand2 = 2.d0 * drand() - 1.d0
      radius = rand1 * rand1 + rand2 * rand2
      if (radius .gt. 1.d0 .or. radius .eq. 0.d0) then
         go to 100
      endif
      pre = sqrt(-2.d0 * log(radius) / radius)
      x = pre * rand1
      y = pre * rand2
      
      return
      end
c      
c-----------------------------------------------------------------------
C R250.F77     The R250 Pseudo-random number generator
C
C algorithm from:
C Kirkpatrick, S., and E. Stoll, 1981; A Very Fast Shift-Register
C Sequence Random Number Generator, Journal of Computational Physics,
C V. 40. p. 517
C 
C see also:
C Maier, W.L., 1991; A Fast Pseudo Random Number Generator,
C                    Dr. Dobb's Journal, May, pp. 152 - 157
C
C 
C Uses the Linear Congruential Method,
C the "minimal standard generator"
C Park & Miller, 1988, Comm of the ACM, 31(10), pp. 1192-1201
C for initialization
C
C
C For a review of BOTH of these generators, see:
C Carter, E.F, 1994; Generation and Application of Random Numbers,
C Forth Dimensions, Vol. XVI, Numbers 1,2 May/June, July/August
C
C
C $Author: stuart $
C $Workfile:   r250.f  $
C $Revision: 1.1.1.1 $
C $Date: 2000/04/05 20:56:42 $
C
C ===================================================================
c
c-----------------------------------------------------------------------
c     The minimal standard PRNG for 31 bit unsigned integers
c     designed with automatic overflow protection  
c     uses ix as the seed value if it is greater than zero
c     otherwise it is ignored
c-----------------------------------------------------------------------
c     Obtained from http://www.taygeta.com/random.html.
c-----------------------------------------------------------------------
c
      function lcmrand(ix)
c
      integer*4 ix
      integer*4 a, b, m, q, r
      integer*4 hi, lo, test
      integer*4 x
c
      save x
c
      parameter (a = 16807, b = 0, m = 2147483647)
      parameter (q = 127773, r = 2836)
c
      if ( ix .gt. 0 ) then 
         x = ix
      endif
      
      hi = x / q
      lo = mod( x, q )
      test = a * lo - r * hi
      if ( test .gt. 0 ) then
          x = test
      else
          x = test + m
      endif
      
      lcmrand = x
      return
      end
c
c-----------------------------------------------------------------------
c  R250, call R250Init with the desired initial seed BEFORE
c  the first invocation of IRAND()
c-----------------------------------------------------------------------
c     Obtained from http://www.taygeta.com/random.html.
c-----------------------------------------------------------------------
c
      subroutine R250Init(iseed)
c
      integer k, mask, msb
      integer indexf, indexb, buffer(250)
c
      common/R250COM/indexf,indexb,buffer
c
      integer ms_bit, all_bits, half_range, step
c
      parameter ( ms_bit = Z'40000000')
      parameter ( half_range = Z'20000000' )
      parameter ( all_bits = Z'7FFFFFFF' )
      parameter ( step = 7 )
c
      indexf = 1
      indexb = 104
      k = iseed
      do 10 i = 1, 250
	  buffer(i) = lcmrand( k )
	  k = -1
  10  enddo
      do 20 i = 1, 250
	 if ( lcmrand( -1 ) .gt. half_range ) then
	     buffer(i) = ior( buffer(i), ms_bit )
	 endif
 20   enddo

      msb = ms_bit
      mask = all_bits

      do 30 i = 0,30
	  k = step * i + 4
	  buffer(k) = iand( buffer(k), mask )
	  buffer(k) = ior( buffer(k), msb )
	  msb = msb / 2
          mask = mask / 2
  30  enddo
      
      return
      end
c
c-----------------------------------------------------------------------
c     R250 PRNG, run after R250_Init
c-----------------------------------------------------------------------
c     Obtained from http://www.taygeta.com/random.html.
c-----------------------------------------------------------------------
c
      function irand()
c
      integer newrand
      integer indexf, indexb, buffer(250)
c
      common/R250COM/indexf, indexb,buffer

      newrand = ieor( buffer(indexf), buffer(indexb) )
      buffer(indexf) = newrand

      indexf = indexf + 1
      if ( indexf .gt. 250 ) indexf = 1

      indexb = indexb + 1
      if ( indexb .gt. 250 ) indexb = 1
	  

      irand = newrand
      return
      end
c
c-----------------------------------------------------------------------
c     double-precision version of irand().  Uses R250 pseudo-random
c     number generator.  Modified version of irand() obtained from
c     http://www.taygeta.com/random.html.  Can be called interchangeably
c     with irand() once R250Init() has been called....sjs 7/7/99
c-----------------------------------------------------------------------
c
      real*8 function drand()
c
      integer newrand
      integer indexf, indexb, buffer(250)
c
      common/R250COM/indexf, indexb,buffer

      newrand = ieor( buffer(indexf), buffer(indexb) )
      buffer(indexf) = newrand

      indexf = indexf + 1
      if ( indexf .gt. 250 ) indexf = 1

      indexb = indexb + 1
      if ( indexb .gt. 250 ) indexb = 1
	  
c     The constant here is 2^31 - 1, which is the size of the 
c     integer*4 elements of the array buffer() which is being
c     twiddled.  This converts them to [0,1)

      drand = dble(newrand) / 2147483648.d0

      return
      end
c




