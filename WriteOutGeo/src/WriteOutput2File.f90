! "Solver" written to output any scalar/vector field at the node points to file
! for interpolation onto a new numerical grid for instance
! Written by Clemens Schannwell
! Last modified: 22.11.2018
SUBROUTINE WriteOutputField2File( Model,Solver,dt,TransientSimulation )
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

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName,filename,scount, OutVarName
  CHARACTER(LEN=MAX_NAME_LEN) :: OutFileName

  TYPE(Element_t),POINTER :: Element
  TYPE(Variable_t), POINTER :: VarSol
  TYPE(ValueList_t), POINTER :: SolverParams
  REAL(KIND=dp), POINTER :: Var (:)
  INTEGER, POINTER :: VarPerm(:)
  INTEGER,SAVE :: VarDOFs,k,nn,IntLength
  REAL(kind=dp) :: x,y,z,U,V,W,P,Partitions
!------------------------------------------------------------------------------
! Real function starts  
!------------------------------------------------------------------------------

   SolverParams => GetSolverParams() ! needed to get solver variables

!------------------------------------------------------------------------------
! Get Variables abort if not found
!------------------------------------------------------------------------------

   OutVarName =  GetString( SolverParams,'Output Variable Name', Found)
   IF(.NOT.Found) THEN
       CALL FATAL(SolverName,'Keyword >Output Variable Name< not found in section >Solver<')
   END IF
   VarSol => VariableGet( Solver % Mesh % Variables, TRIM(OutVarName))
   IF ( ASSOCIATED( VarSol ) ) THEN
           Var => VarSol % Values
           VarPerm => VarSol % Perm
           VarDOFs=VarSol % DOFs
   ELSE
           WRITE(Message,'(A,A,A)') &
                              'No output variable ', TRIM(OutVarName), ' found'
           CALL FATAL(SolverName,Message)
   END IF

   OutFileName =  GetString( SolverParams,'Output File Name', Found)
   IF(.NOT.Found) THEN
       CALL FATAL(SolverName,'Keyword >Output File Name< not found in section >Solver<')
   END IF
!------------------------------------------------------------------------------
! Get number of partitions for file name in parallel/ in serial this should give
! 0
!------------------------------------------------------------------------------
   Partitions = ParEnv % PEs
   IntLength = CEILING(LOG10(Partitions+1.0))
!------------------------------------------------------------------------------
! Write all nodal values to file
!------------------------------------------------------------------------------
   nn=Solver % Mesh % NumberOfNodes
   DO k=1,nn !loop over all nodes
           x   = Solver % Mesh % Nodes % x(k)
           y   = Solver % Mesh % Nodes % y(k)
           z   = Solver % Mesh % Nodes % z(k)
           IF  (VarDOFs.GT.1) THEN ! for vector variable (in this case veloctiy)
                   ! might need adjustment for something like stress
           U   = Var(VarDOFs*(VarPerm(k)-1)+1)
           V   = Var(VarDOFs*(VarPerm(k)-1)+2)
           W   = Var(VarDOFs*(VarPerm(k)-1)+3)
           P   = Var(VarDOFs*(VarPerm(k)-1)+4)
       ELSE ! for scalar variables
           U   = Var(VarPerm(k))
       END IF
           !WRITE(scount,'(i<IntLength>.<IntLength>)') ParEnv % MyPe
           WRITE(scount,'(i6.6)') ParEnv % MyPe
           filename = TRIM(OutFileName)//"."//scount ! generate filename
           OPEN(12,FILE=filename)
        IF  (VarDOFs.GT.1) THEN ! writing for vector variables
           WRITE(12,'(e13.5,2x,e15.8,2x,e15.8,2x,e15.8,2x,e15.8,2x,e15.8,2x,e15.8)') x,y,z,U,V,W,P
        ELSE ! writing for scalar variables
           WRITE(12,'(e13.5,2x,e15.8,2x,e15.8,2x,e15.8)') x,y,z,U
        END IF
   END DO
         CLOSE(12) ! close file
!------------------------------------------------------------------------------
! THE END !
!------------------------------------------------------------------------------

END SUBROUTINE WriteOutputField2File
