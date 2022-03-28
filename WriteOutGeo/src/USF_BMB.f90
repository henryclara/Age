! Get BMB as Function of Ice Thickness, GL Position and Distance.
FUNCTION GetBMB ( Model, nodenumber, BMBMultiplier) RESULT(BMBOut)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Variable_t), POINTER :: DepthSol,DistanceSol
   INTEGER, POINTER :: DepthPerm(:),DistancePerm(:),BMBPerm(:)
   REAL(kind=dp), POINTER :: DepthVal(:),DistanceVal(:)
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber,  NMAX, i, dim
   REAL(KIND=dp) :: BMBMultiplier,   BMBOut, alpha, G,A,rho, rhoInit
   REAL(KIND=dp), ALLOCATABLE :: BMB0(:)
   LOGICAL :: FirstTime=.True., UnFoundFatal,GotIt

   SAVE FirstTime
   SAVE BMB0
   !Get Depthwithout catching any error messages if fields don't exist
   DepthSol => VariableGet( Model % Variables, 'Depth',UnFoundFatal=UnFoundFatal)
   DepthPerm => DepthSol % Perm
   DepthVal => DepthSol % Values

   !Get Distance without catching any error messages if fields don't exist
   DistanceSol => VariableGet( Model % Variables, 'Distance',UnFoundFatal=UnFoundFatal)
   DistancePerm => DistanceSol % Perm
   DistanceVal => DistanceSol % Values

  !alpha = 0.5
  !rho = 1-np.exp(-0.0001*distance) #transition GL to Ambient
  !G = 0.001 ## melting away from GL relative to H^alpha
  !A = 0.1   ##melting near GL relative to H^alpha
  !bmb = thickness**(alpha)*(rho*G+(1-rho)*A)


   alpha = GetConstReal (Model % Constants, 'Alpha', GotIt)
   IF (.NOT.GotIt) THEN
      CALL FATAL('USF BMB', 'No alpha specified in sif file.')
   END IF
   G = GetConstReal (Model % Constants, 'G', GotIt)
   IF (.NOT.GotIt) THEN
      CALL FATAL('USF BMB', 'No G specified in sif file.')
   END IF
   A = GetConstReal (Model % Constants, 'A', GotIt)
   IF (.NOT.GotIt) THEN
      CALL FATAL('USF BMB', 'No A specified in sif file.')
   END IF
   rhoInit = GetConstReal (Model % Constants, 'rho', GotIt)
   IF (.NOT.GotIt) THEN
      CALL FATAL('USF BMB', 'No rho specified in sif file.')
   END IF
   !!alpha = 0.4
   !alpha = 0.76
   !G = 0.01
   !A = 0.06
   rho = 1 - exp(-1.0/rhoInit* DistanceVal(DistancePerm(nodenumber)))

   BMBOut = DepthVal(DepthPerm(nodenumber))**alpha*(rho*G+(1-rho)*A)
   !print *, BMBOut
   !write(*,*) 'BMBOut:',Model % Nodes % x(nodenumber),Model % Nodes % y(nodenumber),BMBOut

END FUNCTION GetBMB



