      SUBROUTINE GroundedMaskSolver( Model,Solver,dt,TransientSimulation )
      USE DefUtils
      USE SolverUtils
      IMPLICIT NONE

      TYPE(Model_t) :: Model
      TYPE(Solver_t):: Solver
      REAL(KIND=dp) :: dt
      LOGICAL :: TransientSimulation

      LOGICAL :: GL
      Type(Element_t), POINTER :: Element
      Type(Variable_t), POINTER :: Zb,Bed
      REAL(KIND=dp),POINTER :: Mask(:)
      INTEGER, POINTER :: Perm(:)
      INTEGER :: t,n,i
      REAL(KIND=dp), parameter :: eps=0.1_dp


      Mask => Solver % Variable % Values
      Perm => Solver % Variable % Perm
      write(*,*) 'Reading Zb..'
      Zb => VariableGet( Solver % Mesh % Variables, 'Zb' )
      IF ( ASSOCIATED( Zb ) ) THEN
            write(*,*) 'Done Zb.'
      ELSE
            write(*,*) 'Zb NOT ASSOCIATED.'
      END IF

      write(*,*) 'Reading bedrock..'
      Bed => VariableGet( Solver % Mesh % Variables, 'bedrock' )
      IF ( ASSOCIATED( Bed ) ) THEN
            write(*,*) 'Done Bedrock.'
      ELSE
            write(*,*) 'BED NOT ASSOCIATED.'
      END IF

! DEFINE INTIAL MASK
! ice sheet ==  1
! ice shelf == -1
      DO t=1,Solver % NumberOfActiveElements
         Element => GetActiveElement(t)
         IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
         n = GetElementNOFNodes()
         Do i=1,n
            if (Zb%Values(Zb%Perm(Element%NodeIndexes(i))).GT.(Bed%Values(Bed%Perm(Element%NodeIndexes(i)))+eps)) then
               Mask(Perm(Element%NodeIndexes(i)))=-1.0_dp
            else
               Mask(Perm(Element%NodeIndexes(i)))=+1.0_dp
            endif
         End do
      End do
     write(*,*) 'Done Initiliazing Mask.'
! GROUNDING LINE NEXT NEIGHBOURS
! on ice sheet side  1/2
! on ice shelf side -1/2
      !DO t=1,Solver % NumberOfActiveElements
         !Element => GetActiveElement(t)
         !IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
         !n = GetElementNOFNodes()
         !GL=.False.
         !IF (ANY(Mask(Perm(Element%NodeIndexes(1:n))).EQ.-1.0_dp)) GL=.True.
         !IF (.NOT.GL) cycle
         !Do i=1,n
           !If (Mask(Perm(Element%NodeIndexes(i))).EQ.1.0_dp) &
                           !Mask(Perm(Element%NodeIndexes(i)))=1.0/2.0
         !End do
         !GL=.False.
         !IF (ANY(Mask(Perm(Element%NodeIndexes(1:n))).EQ.1.0/2.0)) GL=.True.
         !IF (.NOT.GL) cycle
         !Do i=1,n
           !If (Mask(Perm(Element%NodeIndexes(i))).EQ.-1.0_dp) &
                           !Mask(Perm(Element%NodeIndexes(i)))=-1.0/2.0
         !End do
      !End do

!! GROUNDING LINE NEXT NEIGHBOURS SECOND ROUND
!! on ice sheet side  1/3
!! on ice shelf side -1/3
!      DO t=1,Solver % NumberOfActiveElements
!         Element => GetActiveElement(t)
!         IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
!         n = GetElementNOFNodes()
!         GL=.False.
!         IF (ANY(Mask(Perm(Element%NodeIndexes(1:n))).EQ.1.0/2.0)) GL=.True.
!         IF (.NOT.GL) cycle
!         Do i=1,n
!           If (Mask(Perm(Element%NodeIndexes(i))).EQ.1.0_dp) &
!                           Mask(Perm(Element%NodeIndexes(i)))=1.0-1.0/3.0
!         End do
!         GL=.False.
!         IF (ANY(Mask(Perm(Element%NodeIndexes(1:n))).GT.-0.9999)) GL=.True.
!         IF (.NOT.GL) cycle
!         Do i=1,n
!           If (Mask(Perm(Element%NodeIndexes(i))).EQ.-1.0_dp) &
!                           Mask(Perm(Element%NodeIndexes(i)))=-1.0+1.0/3.0
!         End do
!      End do

!! GROUNDING LINE NEXT NEIGHBOURS SECOND ROUND
!! on ice sheet side  1/4
!! on ice shelf side -1/4
!      DO t=1,Solver % NumberOfActiveElements
!         Element => GetActiveElement(t)
!         IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
!         n = GetElementNOFNodes()
!         GL=.False.
!         IF (ANY(Mask(Perm(Element%NodeIndexes(1:n))).EQ.1.0-1.0/3.0)) GL=.True.
!         IF (.NOT.GL) cycle
!         Do i=1,n
!           If (Mask(Perm(Element%NodeIndexes(i))).EQ.1.0_dp) &
!                           Mask(Perm(Element%NodeIndexes(i)))=1.0-1.0/4.0
!         End do
!         GL=.False.
!         IF (ANY(Mask(Perm(Element%NodeIndexes(1:n))).GT.-0.9999)) GL=.True.
!         IF (.NOT.GL) cycle
!         Do i=1,n
!           If (Mask(Perm(Element%NodeIndexes(i))).EQ.-1.0_dp) &
!                           Mask(Perm(Element%NodeIndexes(i)))=-1.0+1.0/4.0
!         End do
!      End do

      CALL InvalidateVariable( Model % Meshes, Solver%Mesh,Solver%Variable%name )
      End


