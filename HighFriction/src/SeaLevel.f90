FUNCTION getSeaLevel( Model, nodenumber, t ) RESULT( SeaLevel )
        USE types
        USE CoordinateSystems
        USE SolverUtils
        USE ElementDescription
        USE DefUtils
        IMPLICIT NONE
        !header variables
        TYPE(Model_t) :: Model
        INTEGER :: nodenumber
        REAL(KIND=dp) :: SeaLevel, t
        
        !Sea level is twice increased and decreased with a negative perturbation
        !at the end

        IF (t < 2000.0) THEN
                SeaLevel = 0.0
        ELSE
                SeaLevel = -40.0 + 0.02*t
        !ELSE IF (t < 12000.0) THEN
                !SeaLevel = 160.0
        !ELSE IF (t < 12000.0) THEN
        !        SeaLevel = 240.0 - 0.02*t
        !ELSE IF (t < 14000.0) THEN
        !        SeaLevel = 0.0
        !ELSE IF (t < 18000.0) THEN
        !        SeaLevel = -280.0 + 0.02*t
        !ELSE IF (t < 20000.0) THEN
        !        SeaLevel = 80.0
        !ELSE IF (t < 26000.0) THEN
        !        SeaLevel = 480.0 - 0.02*t
        !ELSE
        !        SeaLevel = 0.0
        END IF
        ! write(*,*) 'Sea level: ', SeaLevel
END FUNCTION getSeaLevel
