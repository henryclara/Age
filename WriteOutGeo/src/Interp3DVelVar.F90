! "Solver" written to interpolate any scalar/vector field at node points from an
! old mesh to node points of a new mesh (no matter whether it's horizontal or
! vertical resolution changes)
! Interpolation method is inverse distance weighting
! Written by Clemens Schannwell
! Last modified: 22.11.2018
SUBROUTINE Interp3DVelVar( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
  USE ParticleUtils
  USE GeneralUtils
  USE DefUtils
  USE Interpolation
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation, Found

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
  CHARACTER(LEN=MAX_NAME_LEN) :: InVarName, InFileName

  TYPE(Element_t),POINTER :: Element
  TYPE(Variable_t), POINTER :: VarSol
  TYPE(ValueList_t), POINTER :: SolverParams, BC
  REAL(KIND=dp), POINTER :: Var (:)
  INTEGER, POINTER :: VarPerm(:),NodeIndexes(:)
  INTEGER,SAVE :: VarDOFs,k,nn,n,io,i,j
  REAL(kind=dp) :: x,y,z,U,V,W, DistanceHor, DistanceVert, SearchRadius
  REAL(kind=dp) :: AnisoFactor, UInterpNum, UInterpDen, DistanceExponent
  REAL(kind=dp) :: VInterpNum, VInterpDen, WInterpNum, WInterpDen,PInterpNum,PInterpDen
  REAL, DIMENSION(:), ALLOCATABLE :: xOld,yOld,zOld,UOld,VOld,WOld,POld
!------------------------------------------------------------------------------
! Real function starts  
!------------------------------------------------------------------------------

   SolverParams => GetSolverParams() ! needed to get solver variables

!------------------------------------------------------------------------------
! Get Variables abort if not found or set to default values where possible
!------------------------------------------------------------------------------

   InVarName =  GetString( SolverParams,'Input Variable Name', Found)
   IF(.NOT.Found) THEN
       CALL FATAL(SolverName,'Keyword >Input Variable Name< not found in section >Solver<')
   END IF
   VarSol => VariableGet( Solver % Mesh % Variables, TRIM(InVarName))
   IF ( ASSOCIATED( VarSol ) ) THEN
           Var => VarSol % Values
           VarPerm => VarSol % Perm
           VarDOFs=VarSol % DOFs
   ELSE
           WRITE(Message,'(A,A,A)') &
                              'No output variable ', TRIM(InVarName), ' found'
           CALL FATAL(SolverName,Message)
   END IF

   InFileName =  GetString( SolverParams,'Input File Name', Found)
   IF(.NOT.Found) THEN
       CALL FATAL(SolverName,'Keyword >Input File Name< not found in section >Solver<')
   END IF
   SearchRadius = GetConstReal( SolverParams, 'Search Radius', Found )
   IF (.NOT.Found) THEN
            WRITE(Message,'(A)') 'Constant Search Radius not found. &
                   &Setting to 10,000 m'
            CALL INFO(SolverName, Message, level=20)
            SearchRadius = 10000.0
    END IF
   AnisoFactor = GetConstReal( SolverParams, 'Anisotropy Factor', Found )
   IF (.NOT.Found) THEN
            WRITE(Message,'(A)') 'Constant Anisotropy Factor not found. &
                   &Setting to 100 m'
            CALL INFO(SolverName, Message, level=20)
            AnisoFactor = 100.0
    END IF
   DistanceExponent = GetConstReal( SolverParams, 'Distance Weighting', Found )
   IF (.NOT.Found) THEN
            WRITE(Message,'(A)') 'Constant Distance Weighting not found. &
                   &Setting to 1'
            CALL INFO(SolverName, Message, level=20)
            DistanceExponent = 1.0
    END IF
!------------------------------------------------------------------------------
! Do a dry read of the file to determine how much memory needs to be allocated
!------------------------------------------------------------------------------
  OPEN(11, FILE=InFileName, STATUS='old', ACTION='read')
  n = 0 
  DO
    READ(11,*,iostat=io)
    IF (io/=0) EXIT
    n = n + 1
  END DO
  REWIND(11)
!------------------------------------------------------------------------------
! Allocate memory and do the actual read of the variables
!------------------------------------------------------------------------------
  IF  (VarDOFs.GT.1) THEN ! for vector variable (in this case veloctiy)
                   ! might need adjustment for something like stress
    ALLOCATE(xOld(n), yOld(n), zOld(n), UOld(n), VOld(n), WOld(n), POld(n))
    DO i=1, n
     READ (11,'(e13.5,2x,e15.8,2x,e15.8,2x,e15.8,2x,e15.8,2x,e15.8,2x,e15.8)') &
             xOld(i), yOld(i), zOld(i), UOld(i), VOld(i), WOld(i), POld(i)
    END DO
  ELSE ! for scalar variable
    ALLOCATE(xOld(n), yOld(n), zOld(n), UOld(n))
    DO i=1, n
     READ (11,'(e13.5,2x,e15.8,2x,e15.8,2x,e15.8)') xOld(i), yOld(i), zOld(i), UOld(i)
    END DO
  END IF
  nn=Solver % Mesh % NumberOfNodes
  DO k=1,nn ! loop over all nodes
      x   = Solver % Mesh % Nodes % x(k)
      y   = Solver % Mesh % Nodes % y(k)
      z   = Solver % Mesh % Nodes % z(k)
      UInterpNum = 0.0
      VInterpNum = 0.0
      WInterpNum = 0.0
      PInterpNum = 0.0
      UInterpDen = 0.0
      VInterpDen = 0.0
      WInterpDen = 0.0
      PInterpDen = 0.0
      DO i=1,n ! loop over all nodes in file and calc distance
        DistanceHor = SQRT( (x-xOld(i))**2 +(y-yOld(i))**2 )
        DistanceVert = SQRT( (z-zOld(i))**2 )
! Carry out inverse distance interpolation if within Search Radius
        IF &
        (DistanceHor.LE.SearchRadius.AND.DistanceVert.LE.SearchRadius/AnisoFactor) THEN
                UInterpNum = UInterpNum + UOld(i)/(DistanceHor + DistanceVert * &
                AnisoFactor)**DistanceExponent
                UInterpDen = UInterpDen + 1.0/(DistanceHor + DistanceVert * &
                AnisoFactor)**DistanceExponent
                IF  (VarDOFs.GT.1) THEN ! for vector variable (in this case veloctiy)
                   ! might need adjustment for something like stress
                  VInterpNum = VInterpNum + VOld(i)/(DistanceHor + DistanceVert * &
                  AnisoFactor)**DistanceExponent
                  VInterpDen = VInterpDen + 1.0/(DistanceHor + DistanceVert * &
                  AnisoFactor)**DistanceExponent
                  WInterpNum = WInterpNum + WOld(i)/(DistanceHor + DistanceVert * &
                  AnisoFactor)**DistanceExponent
                  WInterpDen = WInterpDen + 1.0/(DistanceHor + DistanceVert * &
                  AnisoFactor)**DistanceExponent
                  PInterpNum = PInterpNum + POld(i)/(DistanceHor + DistanceVert * &
                  AnisoFactor)**DistanceExponent
                  PInterpDen = PInterpDen + 1.0/(DistanceHor + DistanceVert * &
                  AnisoFactor)**DistanceExponent
                END IF
        END IF
      END DO
      ! Abort if no points where found within Search Radius
      IF (UInterpDen.EQ.0.0.AND.UInterpNum.EQ.0.0) THEN
        CALL FATAL(SolverName,'No interpolation points found in Search Radius. &
        Please increase Search Radius')
      END IF
      ! Write interpolated values to new nodal points
      IF  (VarDOFs.GT.1) THEN ! for vector variable (in this case veloctiy)
                   ! might need adjustment for something like stress
        Var(VarDOFs*(VarPerm(k)-1)+1) = UInterpNum/UInterpDen
        Var(VarDOFs*(VarPerm(k)-1)+2) = VInterpNum/VInterpDen
        Var(VarDOFs*(VarPerm(k)-1)+3) = WInterpNum/WInterpDen
        Var(VarDOFs*(VarPerm(k)-1)+4) = PInterpNum/PInterpDen
      ELSE ! for scalar variables
        Var(VarPerm(k)) = UInterpNum/UInterpDen
      END IF
   END DO
      write(*,*) 'VALUE: ', nn
     DO i=1,GetNofBoundaryElements()
        Element=>GetBoundaryElement(i)
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes
        BC=>GetBC()
        IF(.NOT.ASSOCIATED(BC)) CYCLE
        IF(ListCheckPresent(BC,'Velocity 1'))  THEN
            !write(*,*) 'AgeNode: ', NodeIndexes
            !write(*,*) 'VelocityNode: ', NodeIndexes
            !write(*,*) ListGetReal( BC, 'Velocity 2',n,NodeIndexes)
            Var(VarDOFs*(VarPerm(NodeIndexes)-1)+1) = ListGetReal( BC, 'Velocity 1',n,NodeIndexes) 
        END IF
        IF(ListCheckPresent(BC,'Velocity 2'))  THEN
            !write(*,*) 'AgeNode: ', NodeIndexes
            !write(*,*) 'VelocityNode: ', NodeIndexes
            !write(*,*) ListGetReal( BC, 'Velocity 2',n,NodeIndexes)
            Var(VarDOFs*(VarPerm(NodeIndexes)-1)+2) = ListGetReal( BC, 'Velocity 2',n,NodeIndexes) 
        END IF
        IF(ListCheckPresent(BC,'Velocity 3'))  THEN
            !write(*,*) 'AgeNode: ', NodeIndexes
            !write(*,*) 'VelocityNode: ', NodeIndexes
            !write(*,*) ListGetReal( BC, 'Velocity 2',n,NodeIndexes)
            Var(VarDOFs*(VarPerm(NodeIndexes)-1)+3) = ListGetReal( BC, 'Velocity 3',n,NodeIndexes) 
        END IF
     END DO
!------------------------------------------------------------------------------
! THE END !
!------------------------------------------------------------------------------
END SUBROUTINE Interp3DVelVar
