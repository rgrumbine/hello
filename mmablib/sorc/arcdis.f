      FUNCTION arcdis (long1, lat1, long2, lat2)
C     Function to compute the distance between two latitude, longitude
C       points in nautical miles.
C     Algorithm from Doug MacAyeal, The University of Chicago, 1985.
C     Bullet-proofing added by Robert Grumbine 6 April 1994.
C     Returns -1 in event of errors.  Diagnostic printout eliminated.
C       Robert Grumbine 15 September 2005

      IMPLICIT none

      DOUBLE PRECISION ab, ac, bc, arg
      REAL arcdis, lat1, long1, lat2, long2

      DOUBLE PRECISION pi, rdpdg, kmtonm, rearth
      PARAMETER (pi = 3.141592654)
      PARAMETER (rdpdg = pi / 180.)
      PARAMETER (kmtonm = 1.852)
CERROR found 13 September 1995      PARAMETER (rearth = 6730.949)
      PARAMETER (rearth = 6370.949)

C     Bullet-proofing: Check that lats are within +- 90, longs
C       longs are within +- 360.
      IF (ABS(lat1) .GT. 90. .OR. ABS(lat2) .GT. 90.) THEN
CD        PRINT *,'Latitudes are off the planet.'
CD        PRINT *,'No fix within arcdis.'
CD        PRINT *,'WARNING error.'
          arcdis = -1.
          RETURN
      ENDIF
      IF (ABS(long1) .GT. 360. .OR. ABS(long2) .GT. 360.) THEN
CD        PRINT *,'Longitudes are outside range.'
CD        PRINT *,'No fix within arcdis.'
CD        PRINT *,'WARNING long1, long2 = ',long1, long2
          arcdis = -1.
          RETURN
      ENDIF

C     Operational Code      
C     Special case included because trig round off can give identical
C      points a nonzero separating distance. BG 6/3/93.
      IF (long1 .EQ. long2 .AND. lat1 .EQ. lat2) THEN
        arcdis = 0.0
        RETURN
      ENDIF

      ab = (90.-lat1)*rdpdg
      ac = (90.-lat2)*rdpdg
      bc = ABS((long1-long2)*rdpdg)
C     arg introduced 6/3/93 to proof against the special case that
C       round-off errors in trig produce a value for the argument 
C       which is out of bounds for the arc-cosine function.  This
C       should be rare though not totally impossible.
      arg = COS(ab)*COS(ac)+SIN(ab)*SIN(ac)*COS(bc)
      IF (arg .GT.  1.0) THEN 
CD        PRINT *,'Out of bounds arcdis ',lat1, long1, lat2, long2
        arg =  1.0
      ENDIF
      IF (arg .LT. -1.0) THEN 
CD        PRINT *,'Out of bounds arcdis ',lat1, long1, lat2, long2
        arg = -1.0
      ENDIF
      arcdis = ACOS(arg)*rearth

C     Bullet-proofing: Verify that distance is less than circumference
C       of the earth.
      IF (arcdis .GT. 2.*pi*rearth) THEN
CD        PRINT *,'Failure in computation of distance on earth.'
CD        PRINT *,'Derived distance = ',arcdis,' nautical miles.'
CD        PRINT *,'FATAL error'
        arcdis = -1
        RETURN
      ENDIF

CD      arcdis = arcdis/kmtonm

      RETURN
      END
