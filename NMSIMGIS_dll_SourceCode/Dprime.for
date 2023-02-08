      SUBROUTINE DPRIMF (DPSS, TRGSPC, BKGSPC)

      !DEC$ ATTRIBUTES DLLEXPORT::DPRIMF

C      ****************************************************************
C      *                                                              *
C      *  PURPOSE:  To Calculate The Auditory Dectectability Index d' *
C      *            as a Function of 1/3 Octave Band Target and       *
C      *            Background Spectra.                               *
C      *                                                              *
C      *  METHOD:   Uses basic auditory signal detection equation:    *
C      *              D' = ETA * SQRT(BW) * S/N                       *
C      *                 where: ETA = observer efficiency             *
C      *                        BW  = masking bandwidth (usually      *
C      *                              1/3 octave bandwidth)           *
C      *                        S   = target signal intensity,        *
C      *                              exp10(SPL/10)                   *
C      *                        N   = background intensity,           *
C      *                              exp10(SPL/10)                   *
C      *                                                              *
C      *  REFERENCE:                                                  *
C      *            Derived from 1983 Acoustic Range Prediction       *
C      *            Program (ARPP) software developed under USAF      *
C      *            Contract F33615-83-C-3216 by Bolt Beranek         *
C      *            and Newman, Inc (BBN) for:                        *
C      *                                                              *
C      *               FLIGHT DYNAMICS LABORATORY                     *
C      *               AIR FORCE WRIGHT AERONAUTICAL LABORATORIES     *
C      *               AIR FORCE SYSTEMS COMMAND                      *
C      *               WRIGHT-PATTERSON AFB,  OHIO  45433             *
C      *                                                              *
C      *            Documentation may by found in USAF Technical      *
C      *            Report No. AFWAL-TR-83-3115 which contains the    *
C      *            Program Description and Users' Guide.             *
C      *                                                              *
C      *  CALL SEQUENCE:                                              *
C      *            CALL DPRIME (TRGSPC, BKGSPC, IBMAX, DPMAX, DPSS)  *
C      *                                                              *
C      *  CALLING ARGUMENTS:                                          *
C      *            TRGSPC = Target 1/3 octave band source spectrum,  *
C      *                     24 bands, 50-10,000 Hz inclusive.        *
C      *            BKGSPC = Background 1/3 octave band spectrum,     *
C      *                     24 bands, 50-10,000 Hz inclusive.        *
C      *                                                              *
C      *  RETURNED VARIABLES:                                         *
C      *            IBMAX  = Band # with highest value of d'          *
C      *                     (1 = 50 Hz, 24 = 10,000 Hz)              *
C      *            DPMAX  = Highest value of d' in any of the 24     *
C      *                     1/3 octave bands.                        *
C      *            DPSS   = Cumulative value of d' across all 24     *
C      *                     1/3 octave bands (square root of the     *
C      *                     sum of the squares of the individual     *
C      *                     bands).                                  *
C      *            All these values will be zero if no target        *
C      *            spectrum band exceeds the threshold of hearing.   *
C      *                                                              *
C      *                                                              *
C      *  AUTHOR:  R. Horonjeff                                       *
C      *  CREATION DATE:  29 September 2000                           *
C      *  REVISION HISTORY:                                           *
C      *                                                              *
C      *  4 Dec 2004: KJP revised per EASN method defined in e-mail   *
C      *  from R. Horonjeff 3 Dec 2004                                *
C      *                                                              *
C      *  17 Dec 2004: KJP revised to use current hearing threshold.  *
C      *  That is, EARSPC from 2002 ISO standard.  EASN is now based  *
C      *  on d' = 1.75, instead of 1.5.                               *
C
C         4 February 2016: BJI modified to be called from Python    
C
C      *                                                              *
C      ****************************************************************
C
C
      implicit none

!     Input 
      real,intent(in) :: TRGSPC(24), BKGSPC(24)

!     Output
      real,intent(inout) :: DPSS

!     Internal
      real dpnoise, s, DPMAX, ETABW, ASL, BKGND
      real SIG, BNOISE, DP
      real EARSPC(24), easn(24)

      integer i, IBMAX
C
C     Threshold of Hearing Array
C       ISO (1961). "Normal Equal-Loudness Contours for Pure Tones
C       and Normal Threshold of Hearing Under Free Field Listening
C       Conditions," International Organization for
C       Standardization, ISO R226-1961(E), Switzerland.
c
c kjp Dec 04: Corrected #19 from -1.5 to -3.5
C
c     DATA  EARSPC / 41.7, 35.6, 29.8, 25.1, 20.8, 16.8, 13.8, 11.2,
c    1                9.0,  7.2,  6.0,  5.0,  4.4,  4.2,  3.8,  2.6,
c    2                1.0, -1.2, -3.5, -3.9, -1.1,  6.5, 15.3, 16.4 /
c
c
c kjp Dec 04: Insert new threshold values, from ISO 389-7
c
      DATA  EARSPC / 43.5, 37.5, 31.5, 26.5, 22.0, 18.0, 14.5, 11.0,
     1                8.5,  5.5,  3.5,  1.5,  1.0,  0.5,  0.0, -1.0,
     2               -1.5, -3.0, -4.5, -5.0, -3.5,  0.5,  5.3,  9.5 /
c
c Compute EASN, the effective noise associated with EARSPC.  Could do this
c within the 30 loop, and not even need for easn to be an array, but it's
c clearer like this.
c
      dpnoise = 1.75
      do i = 1,24
      s = 10.**(earspc(i)*.1)
      easn(i) = etabw(i,0.0)*s/dpnoise
      enddo

C     Initialize Variables:
C         Band number with greatest d',
C         Greatest 1/3 OB value of d',
C         Cumulative value of d'.
C
      IBMAX = 0
      DPMAX = 0.
      DPSS = 0.
C
C          Go Thru the Spectrum Band-by-Band Calculating Values of d'
C          for Those Bands Whose Target Sound Level Exceeds the ISO
C          Threshold of Hearing.
ckjp Dec 04: Do all bands, including EASN in noise instead of just skipping

      DO 30  I=1,24
      ASL   = TRGSPC(I)
      BKGND = BKGSPC(I)
C
c Old method
c     IF (ASL .GE. EARSPC(I)) THEN
c        DP = ETABW(I,BKGND) * 10.**((ASL - BKGND) / 10.)
c        DPSS = DPSS + DP**2
c        IF (DP .GT. DPMAX) THEN
c           IBMAX = I
c           DPMAX = DP
c        ENDIF
c     ENDIF
c
c New method - noise is combination of background noise and ear noise
c
      sig = 10.**(.1*asl)
      bnoise = 10**(.1*bkgnd)
c     tot = bnoise+easn(i)
c     totlev = 10.*log10(tot)

      dp = etabw(i,bkgnd)*sig/(bnoise+easn(i))
c     dp = etabw(i,totlev)*sig/(bnoise+easn(i))
c     if(bkgspc(1).lt.-90.)write(69,*)i,bkgspc(i),sig,dp
      dpss = dpss + dp**2
      if(dp.gt.dpmax) then
           ibmax = i
           dpmax = dp
        endif
C
   30 CONTINUE
C
      DPSS = SQRT(DPSS)
      RETURN
      END

c
      REAL FUNCTION ETABW(I, BKGLVL)
C
C      ****************************************************************
C      *                                                              *
C      *  PURPOSE:  To Generate the Product of Detector Efficiency    *
C      *            (ETA) and the Square Root of the Bandwidth As a   *
C      *            Function of Band Number and Background Sound      *
C      *            Level                                             *
C      *                                                              *
C      *  METHOD:   The Auditory Masking Bandwidth Can Be Greater     *
C      *            Than a 1/3 Octave Band, Especially at Low         *
C      *            Frequencies and High Background Sound Levels      *
C      *            Most of the Code Has To Do With Determining When  *
C      *            This Occurs and Making Appropriate Adjustments    *
C      *                                                              *
C      *  CALL SEQUENCE:                                              *
C      *            VAR = ETABW (I, BKGLVL)                           *
C      *                                                              *
C      *  CALLING ARGUMENTS:                                          *
C      *            I      = Band Number (1 = 50 Hz, 24 = 10 KHz)     *
C      *            BKGLVL = 1/3 Octave Band Background SPL           *
C      *                                                              *
C      ****************************************************************
C
C
      DIMENSION   NP(24), TLETA(24), TLBW(2,6,24)
      DIMENSION   TLBW1(12), TLBW2(12), TLBW3(12), TLBW4(12)
      DIMENSION   TLBW5(12), TLBW6(12), TLBW7(12), TLBW8(12)
      DIMENSION   TLBW9(12), TLBW10(12), TLBW11(12), TLBW12(12)
      DIMENSION   TLBW13(12), TLBW14(12), TLBW15(12), TLBW16(12)
      DIMENSION   TLBW17(12), TLBW18(12), TLBW19(12), TLBW20(12)
      DIMENSION   TLBW21(12), TLBW22(12), TLBW23(12), TLBW24(12)
C
      EQUIVALENCE (TLBW(1,1,1),TLBW1(1))
      EQUIVALENCE (TLBW(1,1,2),TLBW2(1))
      EQUIVALENCE (TLBW(1,1,3),TLBW3(1))
      EQUIVALENCE (TLBW(1,1,4),TLBW4(1))
      EQUIVALENCE (TLBW(1,1,5),TLBW5(1))
      EQUIVALENCE (TLBW(1,1,6),TLBW6(1))
      EQUIVALENCE (TLBW(1,1,7),TLBW7(1))
      EQUIVALENCE (TLBW(1,1,8),TLBW8(1))
      EQUIVALENCE (TLBW(1,1,9),TLBW9(1))
      EQUIVALENCE (TLBW(1,1,10),TLBW10(1))
      EQUIVALENCE (TLBW(1,1,11),TLBW11(1))
      EQUIVALENCE (TLBW(1,1,12),TLBW12(1))
      EQUIVALENCE (TLBW(1,1,13),TLBW13(1))
      EQUIVALENCE (TLBW(1,1,14),TLBW14(1))
      EQUIVALENCE (TLBW(1,1,15),TLBW15(1))
      EQUIVALENCE (TLBW(1,1,16),TLBW16(1))
      EQUIVALENCE (TLBW(1,1,17),TLBW17(1))
      EQUIVALENCE (TLBW(1,1,18),TLBW18(1))
      EQUIVALENCE (TLBW(1,1,19),TLBW19(1))
      EQUIVALENCE (TLBW(1,1,20),TLBW20(1))
      EQUIVALENCE (TLBW(1,1,21),TLBW21(1))
      EQUIVALENCE (TLBW(1,1,22),TLBW22(1))
      EQUIVALENCE (TLBW(1,1,23),TLBW23(1))
      EQUIVALENCE (TLBW(1,1,24),TLBW24(1))
C
C          1/3 Octave Band Values of 10*Log10(eta)
C             where: eta = observer efficiency (typically 0.2 to 0.4)
C
      DATA   TLETA              /                -6.96, -6.26, -5.56,
     1        -5.06, -4.66, -4.36, -4.16, -3.96, -3.76, -3.56, -3.56,
     2        -3.56, -3.56, -3.56, -3.76, -3.96, -4.16, -4.36, -4.56,
     3        -4.96, -5.36, -5.76, -6.26, -6.86 /
C
C          Masking Bandwidths as a Function of Frequency Band and
C          Background Spectrum Level (dB/Hz)
C
      DATA    NP     /  2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5,
     1                  5, 5, 3, 3, 4, 4, 4, 4, 4, 3, 3 /
C
      DATA    TLBW1  /  40., 12.1, 100., 21.1, 8*0.0 /
      DATA    TLBW2  /  40., 12.3, 100., 21.3, 8*0.0 /
      DATA    TLBW3  /  40., 12.6, 100., 21.6, 8*0.0 /
      DATA    TLBW4  /  40., 13.1,  50., 14.4, 60., 15.8, 120., 24.8,
     1                    4*0.0 /
      DATA    TLBW5  /  40., 13.5,  50., 14.6, 60., 16.1, 120., 25.1,
     1                    4*0.0 /
      DATA    TLBW6  /  40., 14.1,  50., 15.1, 60., 16.4, 120., 25.4,
     1                    4*0.0 /
      DATA    TLBW7  /  40., 14.7,  50., 15.5, 60., 16.8, 120., 25.8,
     1                    4*0.0 /
      DATA    TLBW8  /  40., 15.4,  50., 16.2, 60., 17.2, 120., 26.2,
     1                    4*0.0 /
      DATA    TLBW9  /  40., 16.1,  50., 16.8, 60., 17.8, 120., 26.8,
     1                    4*0.0 /
      DATA    TLBW10 /  40., 16.9,  50., 17.6, 60., 18.4, 120., 27.4,
     1                    4*0.0 /
      DATA    TLBW11 /  40., 17.6,  50., 18.3, 60., 19.1, 120., 28.1,
     1                    4*0.0 /
      DATA    TLBW12 /  40., 18.4,  50., 19.2, 60., 20.1, 120., 29.1,
     1                    4*0.0 /
      DATA    TLBW13 /  30., 19.1,  40., 19.3, 50., 20.2,  60., 21.2,
     1                    120., 30.2,  2*0.0 /
      DATA    TLBW14 /  30., 19.8,  40., 20.2, 50., 21.2,  60., 22.3,
     1                    120., 31.3,  2*0.0 /
      DATA    TLBW15 /  30., 20.5,  40., 21.3, 50., 22.3,  60., 23.6,
     1                    120., 32.6,  2*0.0 /
      DATA    TLBW16 /  30., 21.4,  40., 22.1, 120., 34.1,  6*0.0 /
      DATA    TLBW17 /  30., 22.1,  40., 23.5, 120., 35.5,  6*0.0 /
      DATA    TLBW18 /  20., 22.9,  30., 23.4, 40., 24.8, 120., 36.8,
     1                    4*0.0 /
      DATA    TLBW19 /  20., 24.0,  30., 24.7, 40., 26.2, 120., 38.2,
     1                    4*0.0 /
      DATA    TLBW20 /  20., 25.1,  30., 26.2, 40., 27.6, 120., 39.6,
     1                    4*0.0 /
      DATA    TLBW21 /  20., 26.2,  30., 27.5, 40., 29.0, 120., 41.0,
     1                    4*0.0 /
      DATA    TLBW22 /  20., 27.5,  30., 28.9, 40., 30.4, 120., 42.4,
     1                    4*0.0 /
      DATA    TLBW23 /  20., 28.8,  40., 31.8, 120., 43.8,  6*0.0 /
      DATA    TLBW24 /  20., 30.1,  40., 33.1, 120., 45.1,  6*0.0 /
C
C
C          Calculate 10*LOG(BANDWIDTH) for 1/3 and 1/1 Octave Bands
C             Center Frequency Multipliers are:
C               1/3OB = 10^(+0.5/10) - 10^(-0.5/10) = 0.230768
C               1/1OB = 10^(+1.5/10) - 10^(-1.5/10) = 0.704592
C
      TL13BW = FLOAT(I+16) - 6.368253
      TL11BW = FLOAT(I+16) - 1.520624
C
C          Calculate Average Spectrum Level in Band By Subtracting
C          10*LOG(BANDWIDTH) From 1/3 Octave SPL
C
      BN0 = BKGLVL - TL13BW
C
C          Look Up 10*LOG(MASKING BANDWIDTH) As A Function of
C          Spectrum Level (N0)
C
      TLMBW = CURVE (BN0, TLBW(1,1,I), NP(I))
C
C          Bound the Masking Bandwidth by the 1/3 Octave Bandwidth
C          (lower bound) and the 1/1 Octave Bandwidth (upper bound)
C          to form usable bandwidth.
C
      TLUBW = AMIN1 (AMAX1(TLMBW,TL13BW), TL11BW)
C
C          Calculate ETA * SQRT(BANDWIDTH)
C
      ETABW = 10. ** ((TLETA(I) + 0.5*TLUBW) / 10.)
      RETURN
      END
c
      REAL FUNCTION CURVE (XL, DUMY, NPTS)
C
C      ****************************************************************
C      *                                                              *
C      *  PURPOSE:  To Look Up a Y value as a Function of X along a   *
C      *            Curve Defined By a Series of X,Y Points Connected *
C      *            By Straight Line Segments                         *
C      *                                                              *
C      *  METHOD:   Simple Linear Interpolation Between Coordinate    *
C      *            Pairs Whose X values Bound the Calling Argument   *
C      *            If X is Less Than the First X-coordinate in the   *
C      *            Curve Array, Then the First Y-coordinate is       *
C      *            Returned.  If X is Greater Than the Last          *
C      *            X-coordinate, Then the Last Y-coordinate is       *
C      *            Returned.                                         *
C      *                                                              *
C      *  CALL SEQUENCE:                                              *
C      *            VAR = CURVE (XL, DUMY, NPTS)                      *
C      *                                                              *
C      *  CALLING ARGUMENTS:                                          *
C      *            XL    = Lookup Value of X                         *
C      *            DUMY  = Two-Dimensional Array Defining the        *
C      *                     Curve                                    *
C      *            NPTS  = Number of Points Defining the Curve       *
C      *                                                              *
C      *  NOTE:                                                       *
C      *            Please pardon the antiquated IF statements.  The  *
C      *            routine is over 20 years old, but it works!       *
C      *                                                              *
C      ****************************************************************
C
C
      DIMENSION  DUMY(2,10)
C
C
      J=1
      IF (NPTS - 1)  6, 5, 1
    1 IF (XL - DUMY(1,J))  5, 5, 2
C
    2 DO 3  J=2,NPTS
      IF (XL - DUMY(1,J))  4, 5, 3
    3 CONTINUE
      GO TO 5
C
    4 I = J - 1
      CURVE = (DUMY(2,J) - DUMY(2,I))*(XL - DUMY(1,I))/(DUMY(1,J) -
     1  DUMY(1,I)) + DUMY(2,I)
      RETURN
C
    5 CURVE = DUMY(2,J)
      RETURN
C
    6 CURVE = 0.
      RETURN
      END
