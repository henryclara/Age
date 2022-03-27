#!/bin/bash
###############################################################################
### This function creates an Elmer sif file for the intial remeshing ##########
###############################################################################
CreateElmerSifForwardInitPS() {
				SifFileName=$1
                BumpApl=$2
				Sigmax=$3
				Sigmay=$4
				xCenter=$5
				yCenter=$6
				ZsInit=$7
				ZbInit=$8
				BedInit=$9
				RiseIntoShelf=${10}
				MeshLayers=${11}
				Rhoi=${12}
				Rhow=${13}
				RateFactor=${14}
				GlenExponent=${15}
				BasalFrictionCoeff=${16}
				SMB=${17}
				FluxInit=${18}
				OutputName=${19}
				IceTemp=${20}
				Incline=${21}
				InclineShelf=${22}
				Alpha=${23}
				G=${24}
				A=${25}
				rho=${26}
                SlidStr=${27}
                SlidExp=${28}
				rm ${SifFileName}
cat > "${SifFileName}" << EOF
!!--------------------------------------------------------!!
!  Island ice rise setup for initial step
!!--------------------------------------------------------!!

check keywords warn
!
! working units are MPa, a, m
!
\$yearinsec = 365.25*24*60*60
\$rhoi = ${Rhoi}/(1.0e6*yearinsec^2)
\$rhow = ${Rhow}/(1.0e6*yearinsec^2)
\$A = ${RateFactor}*yearinsec*1.0e18
\$n = ${GlenExponent}
\$eta = 1.0/(2.0*A)^(1.0/n)
\$gravity = -9.8*yearinsec^2
\$C = ${BasalFrictionCoeff}/(1.0e6*yearinsec^(1.0/n))

\$Bedrock=$BedInit
\$Incline=$Incline
\$InclineShelf=$InclineShelf

$ function BedTopo(xcoord) import Bedrock,InclineShelf {\\
				_BedTopo = Bedrock + 1/100.0 * xcoord * tan(InclineShelf*pi/180)\\
}

Header
  Mesh DB "." "Mesh"
End

Constants
  Water Density = Real \$rhow
  Gas Constant = Real 8.314 !Joule/mol x  K
	Bump Amplitude = Real $BumpApl
	Sigmax = Real $Sigmax
	Sigmay = Real $Sigmay
  x0 = Real $xCenter
	y0 = Real $yCenter
	ZbInit = Real $ZbInit
	ZsInit = Real $ZsInit
	RiseIntoShelf = Real $RiseIntoShelf
	Incline = Real \$Incline
	Alpha = Real $Alpha
	G = Real $G
	A = Real $A
	rho = Real $rho
  ! For SeaSpring/SeaPressure
End

!---------------------------------------------------
!---------------- SIMULATION -----------------------
!---------------------------------------------------

Simulation
  Coordinate System  = Cartesian 3D
  Simulation Type = transient
  Extruded Mesh Levels = Integer $MeshLayers

  Timestepping Method = "bdf"
  BDF Order = 1
  Timestep Intervals = 1
  Output Intervals = 1
  Timestep Sizes = 0.1

  Initialize Dirichlet Conditions = Logical False
  Steady State Max Iterations = 6
  Steady State Min Iterations = 4

	!Specify name of result file. Used for restarts!!
  Output File = "${OutputName}FORMAT.result"
  max output level = 30
End

!---------------------------------------------------
!---------------- BODIES ---------------------------
!---------------------------------------------------

! the ice
Body 1
  Name = "ice"
  Equation = 1
  Body Force = 1
  Material = 1
  Initial Condition = 1
End

! The upper surface
Body 2
  Name= "top free surface"
  Equation = 2
  Material = 1
  Body Force = 2
  Initial Condition = 2
End

! the lower surface
Body 3
  Name= "free surface sea/ice-shelf"
  Equation = 3
  Material = 1
  Body Force = 3
  Initial Condition = 3
End

!---------------------------------------------------
!---------------- INITIAL CONDITIONS ---------------
!---------------------------------------------------

!! for ice
Initial Condition 1
  FlowVar 1 = Real 0.0
  FlowVar 2 = Real 0.0
  FlowVar 3 = Real 0.0
  V 1 = Real 0.0
  V 2 = Real 0.0
  V 3 = Real 0.0
	BedInit = Variable Coordinate 1
		Real MATC "BedTopo(tx)"
	FluxInit = Real $FluxInit
	BedBump = Variable BedInit
			Real Procedure "src/BedrockBump" "BedrockBump"
  Age = Real 0.0
End

!! for top free surface
Initial Condition 2
	Zs = Variable BedBump
			Real Procedure "src/BedrockBump" "ZsAdj"
End

!! for free surface sea/ice-shelf
Initial Condition 3
	Zb = Variable BedBump
			Real Procedure "src/BedrockBump" "ZbAdj"
End

!---------------------------------------------------
!---------------- BODY FORCES ----------------------
!---------------------------------------------------

Body Force 1
  Flow BodyForce 1 = Real 0.0
  Flow BodyForce 2 = Real 0.0
  Flow BodyForce 3 = Real \$gravity
End

!! accumulation flux in m/year
Body Force 2
   Zs Accumulation Flux 1 = Real 0.0e0
   Zs Accumulation Flux 2 = Real 0.0e0 !m/a
   Zs Accumulation Flux 3 = Real $SMB
End

!! no melting/accretion under ice/shelf
Body Force 3
  !Zb Accumulation = Real $BMB
	Zb Accumulation = Variable Time
			Real Procedure "src/USF_BMB" "GetBMB"
End

!---------------------------------------------------
!---------------- MATERIALS ------------------------
!---------------------------------------------------

!! ice material properties in MPa - m - a system
Material 1
  Viscosity Model = String "power law"
  Density = Real \$rhoi
  Viscosity = Real \$eta
  Viscosity Exponent = Real \$1.0/n
  Critical Shear Rate = Real 1.0e-15

  Sea level = Real 0.0

  Glen Enhancement Factor = Real 1.0
! the temperature to switch between the
! two regimes in the flow law
  Limit Temperature = Real -10.0
! In case there is no temperature variable
  Constant Temperature = Real $IceTemp

  Min Zs = Variable "Bottom Zb"
    Real MATC "tx + 10.0"
  Max Zs = Real 1.0e6

  !! Bed condition
  Min Zb = Equals BedBump
  Max Zb = Real 1.0e6
  Cauchy = Logical True
End

!---------------------------------------------------
!---------------- SOLVERS --------------------------
!---------------------------------------------------
!! Initialisation of the Grounded Mask
Solver 1
  !Exec Solver = Never
  Exec Solver = Before All
  Equation = "MapCoordinate"
  Procedure = "StructuredMeshMapper" "StructuredMeshMapper"

  Active Coordinate = Integer 3
  Mesh Velocity Variable = String "dSdt"
  Mesh Update Variable = String "dS"
  Mesh Velocity First Zero = Logical True

  Top Surface Variable Name = String "Zs"
  Bottom Surface Variable Name = String "Zb"

  Displacement Mode = Logical False
  Correct Surface = Logical True
  Minimum Height = Real 1.0
End

Solver 2
  !Exec Solver = Never
  Equation = GroundedMaskIni
  Procedure = "ElmerIceSolvers" "GroundedSolver"
  Variable = GroundedMask
  Variable DOFs = 1

  Toler = Real 1.0e-3
  Bedrock Variable = String "BedBump"
End


Solver 3
  !Exec Solver = Never
  Equation = "NormalVector"
  Procedure = "ElmerIceSolvers" "ComputeNormalSolver"
  Variable = String "Normal Vector"
  Variable DOFs = 3

  ComputeAll = Logical False
  Optimize Bandwidth = Logical False
End

Solver 4
  !Exec Solver = Never
  Equation = Fw
  Procedure = "ElmerIceSolvers" "GetHydrostaticLoads"
  Variable = Fw[Fwater:3]
  Variable DOFs = 3
End

Solver 5
  Exec Solver = Before Simulation
  Equation = "Velocity Preconditioning"
  Procedure = "VelocityPrecond" "VelocityPrecond"
  Variable = -dofs 3 "V"
  Optimize Bandwidth = Logical True
  Linear System Row Equilibration = Logical True
  Element = "p:1 b:4"
  Bubbles in Global System = False
  Linear System Row Equilibration = Logical True
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStabL
  Linear System Max Iterations = 1000
  Linear System Convergence Tolerance = 1.0e-3
  Linear System Preconditioning = ILU0
  Skip Compute Nonlinear Change = Logical True
  Back Rotate N-T Solution = Logical False
End

Solver 6
  Exec Solver = Before Simulation
  Equation = "Pressure Preconditioning"
  Procedure = "PressurePrecond" "PressurePrecond"
  Variable = -dofs 1 "p" ! that would be default
  Element = "p:1 b:4"
  Bubbles in Global System = False
  Linear System Solver = iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 1000
  Linear System Convergence Tolerance = 1.0e-6
  Linear System Preconditioning = None
  Skip Compute Nonlinear Change = Logical True
  Back Rotate N-T Solution = Logical False
End

Solver 7
  Equation = String "Stokes"
  Procedure = "ParStokes" "StokesSolver"
  Variable = "FlowVar"
	Flow Solver Name = String "FlowVar"
  Variable Dofs = 4
  Element = "p:1 b:4"
	Convective = Logical True
	Block Diagonal A = Logical True
	Use Velocity Laplacian = Logical True
	!-----------------------------------------------
  ! Keywords related to the block preconditioning
  !-----------------------------------------------
	! Option 1
  Block Preconditioning = Logical True
  Linear System Adaptive Tolerance = Logical True
  Linear System Base Tolerance = Real 1.0e-3
  Linear System Relative Tolerance = Real 1.0e-2

  Steady State Convergence Tolerance = 1.0E-03
  Nonlinear System Convergence Tolerance = 1.0E-05
  Nonlinear System Max Iterations = 35
  !Nonlinear System Min Iterations = 10
  Nonlinear System Min Iterations = 5
  Nonlinear System Newton After Iterations = 35
  Nonlinear System Newton After Tolerance =  5.0E-04
	! Option 2
	!  Block Preconditioning = Logical True
	! Linear System GCR Restart = Integer 200
	! Linear System Max Iterations = 200
	! Linear System Convergence Tolerance = 1.0e-6
	! !Linear System Convergence Tolerance = 1.0e-9
	! Nonlinear System Max Iterations = 100
	! Nonlinear System Convergence Tolerance = 1.0e-5
	! !Nonlinear System Convergence Tolerance = 1.0e-7
	! Nonlinear System Newton After Tolerance = 1.0e-5


  Exported Variable 1 = FlowVar Loads[Stress Vector:3 CEQ Residual:1]
  Calculate Loads = Logical True
  Exported Variable 2 = -dofs 1 "dSdt"
  Exported Variable 3 = -dofs 1 "dS"
  Exported Variable 4 = -dofs 1 "BedInit"
  Exported Variable 5 = -dofs 1 "RiseIntoShelf"
  Exported Variable 6 = -dofs 1 "BedBump"
  Exported Variable 7 = -dofs 1 "FluxInit"
  Exported Variable 8 = -dofs 1 "Mesh Velocity"
  Exported Variable 9 = -dofs 1 "Mesh Change"
End

Solver 8
  Exec Solver = Never
  Equation = String "StressSolver"
  Procedure =  File "ElmerIceSolvers" "ComputeDevStress"
  ! this is just a dummy, hence no output is needed
  Variable = -nooutput "Sij"
  Variable DOFs = 1
  ! the name of the variable containing the flow solution
  !(U,V,W,Pressure)
  Flow Solver Name = String "FlowVar"
  ! no default value anymore for "Stress Variable Name"
  Stress Variable Name = String "Stress"
  !-----------------------------------------------------------------------
   Exported Variable 1 = "Stress"
   Exported Variable 1 DOFs = 6   ! 4 in 2D, 6 in 3D
   Linear System Solver = "Iterative"
   Linear System Iterative Method = "BiCGStab"
   Linear System Max Iterations = 300
   Linear System Convergence Tolerance = 1.0E-09
   Linear System Abort Not Converged = True
   Linear System Preconditioning = "ILU3"
   Linear System Residual Output = 1
End

Solver 9
  !Exec Solver = Never
  Exec Solver = Before All
  Equation = "HeightDepth"
  Procedure = "StructuredProjectToPlane" "StructuredProjectToPlane"
  Active Coordinate = Integer 3
  Dot Product Tolerance = Real 1.0e-3

  Operator 1 = Depth
  Operator 2 = Height
  Variable 3 = Zb
  Operator 3 = Bottom
End

Solver 10
  !Exec Solver = Never
   Equation = "SolveDistance"

   Procedure = "src/DistanceSolveRD" "DistanceSolver1"
   Variable = Distance

   H scale = real 2
   Distance Pseudo DT = Real 100
! Nonlinear System Relaxation Factor = 0.25

   Nonlinear System Max Iterations = 50
   Nonlinear System Convergence Tolerance = 1.0e-5

 ! Linear System Solver = Direct
 ! Linear System Direct Method = UMFPack
   Linear System Solver = "Iterative"
   Linear System Iterative Method = "BiCGStab"
   Linear System Max Iterations = 300
   Linear System Convergence Tolerance = 1.0E-09
   Linear System Abort Not Converged = False
   Linear System Preconditioning = "ILU1"
   Linear System Residual Output = 1
   Steady State Convergence Tolerance = 1.0e-4

   Dummy Distance Computation = Logical False

End

Solver 11
  !Exec Solver = Never
  Equation = "Free Surface Top"
  Procedure =  "./src/MyFreeSurfaceSolver" "FreeSurfaceSolver"
  !Procedure =  "FreeSurfaceSolver" "FreeSurfaceSolver"
  Variable = "Zs"
  Flow Solver Name = String "FlowVar"
  Flow Loads Name = String "FlowVar Loads"
  Variable DOFs =  1
  Exported Variable 1 = "Zs Residual"
  Exported Variable 1 DOFs = 1

  !Before Linsolve = "EliminateDirichlet" "EliminateDirichlet"

  Linear System Solver = Iterative
  !Linear System Direct Method = UMFPACK
  Linear System Max Iterations = 1500
  Linear System Iterative Method = BiCGStab
  Linear System Preconditioning = ILU0
  Linear System Convergence Tolerance = Real 1.0e-6
  Linear System Abort Not Converged = False
  Linear System Residual Output = 1

  Nonlinear System Max Iterations = 100
  Nonlinear System Convergence Tolerance  = 1.0e-5
  Nonlinear System Relaxation Factor = 1.00

  Steady State Convergence Tolerance = 1.0e-03

  Stabilization Method = Stabilized
  Apply Dirichlet = Logical True

  Relaxation Factor = Real 1.0
End

Solver 12
  !Exec Solver = Never
  Equation = "Free Surface Sea/Shelf"
  Procedure =  "FreeSurfaceSolver" "FreeSurfaceSolver"
  Flow Solver Name = String "FlowVar"
  Flow Loads Name = String "FlowVar Loads"
  Variable = "Zb"
  Variable DOFS =  1
  Exported Variable 1 = "Zb Residual"
  Exported Variable 1 DOFs = 1

  Nonlinear Update Exported Variables = Logical True

  Exported Variable 2 = "Zb Accumulation "
  Exported Variable 2 DOFS = 1

  !Before Linsolve = "EliminateDirichlet" "EliminateDirichlet"

  Linear System Solver = Iterative
  Linear System Direct Method = UMFPACK
  Linear System Max Iterations = 1500
  Linear System Iterative Method = BiCGStab
  Linear System Preconditioning = ILU0
  Linear System Convergence Tolerance = Real 1.0e-6
  Linear System Abort Not Converged = False
  Linear System Residual Output = 1

  Nonlinear System Max Iterations = 100
  Nonlinear System Convergence Tolerance  = 1.0e-5
  Nonlinear System Relaxation Factor = 1.00

  Steady State Convergence Tolerance = 1.0e-03

  Stabilization Method = Stabilized
  Apply Dirichlet = Logical True

  Relaxation Factor = Real 1.0
End

Solver 13
  Exec Solver = Never
  Equation = "Age Equation"

  Variable = String "Age"
  Variable DOFs =  1

  Flow Solution Name = String "Flow Solution"
  ! Linear System Solver = Iterative
  ! Linear System Max Iterations = 1000
  ! Linear System Iterative Method = Diagonal
  ! Linear System Preconditioning = NoNe
  ! Linear System Convergence Tolerance = Real 1.0e-6
  ! Linear System Abort Not Converged = False
  ! Linear System Residual Output = 0
  Linear System Solver = Direct
  Linear System Direct Method = mumps
  mumps percentage increase working space = integer 1000

  Procedure = "./src/AgeSolverRD" "AgeSolver"
  Exported Variable 1 = -dofs 1 "age"
End

Solver 14
  Exec Solver = After Saving
  Equation = "result output"
  Procedure = "ResultOutputSolve" "ResultOutputSolver"
  Save Geometry Ids = Logical True ! add this line if you want to access boundaries in Paraview
  Output File Name = File "${OutputName}FORMAT"
  Output Format = String vtu
End


!---------------------------------------------------
!---------------- EQUATIONS ------------------------
!---------------------------------------------------

Equation 1
  Active Solvers (8) = 1 3 5 6 7 8 9 13 14
  Flow Solution Name = String "FlowVar"
  Convection = String Computed
End

Equation 2
  Active Solvers(1) = 11
  Flow Solution Name = String "FlowVar"
  Convection = String Computed
End

Equation 3
  Active Solvers(4) = 2 4 10 12 
  Flow Solution Name = String "FlowVar"
  Convection = String Computed
End

!---------------------------------------------------
!---------------- BOUNDARY CONDITIONS --------------
!---------------------------------------------------

!! Back
Boundary Condition 1
  Name = "back"
  Target Boundaries = 1
	FlowVar 1 = Variable FluxInit
		Real Procedure "src/BedrockBump" "IceFluxAtBack"
	FlowVar 2 = Real 0.0
	V 1 = Variable FluxInit
		Real Procedure "src/BedrockBump" "IceFluxAtBack"
	V 2 = Real 0.0
  Age = Real 0.0
End

Boundary Condition 2
  Name = "Looking Downhill Left"
  Target Boundaries = 2

  FlowVar 2 = Real 0.0
  V 2 = Real 0.0
End

!! BC Lateral Ice-Shelf (air or sea contact)
Boundary Condition 3
  Name = "front"
  Target Boundaries = 3


  External Pressure = Variable Coordinate 3
     Real Procedure "ElmerIceUSF" "SeaPressure"

  Compute Sea Pressure = Logical True
  ComputeNormal = Logical False

End

Boundary Condition 4
  Name = "Looking Downhill Right"
  Target Boundaries = 4

  FlowVar 2 = Real 0.0
  V 2 = Real 0.0
End

Boundary Condition 5
  Name = "bottom"
  Target Boundaries = 5
  Body Id = 3

  Normal-Tangential FlowVar = Logical True
  Normal-Tangential V = Logical True
  Flow Force BC = Logical True

!
! Condition where the bed is stuck
!
  Zb = Equals BedBump
  Zb Condition = Variable GroundedMask
    Real MATC "tx + 0.5"
!
! Bedrock conditions
!
  Slip Coefficient 2 = Variable Coordinate 1
    Real Procedure "src/USF_Contact" "SlidCoef_Contact"
  Slip Coefficient 3 = Variable Coordinate 1
    Real Procedure "src/USF_Contact" "SlidCoef_Contact"
  Flow Solution Name = String "FlowVar"
  Sliding Law = String "${SlidStr}"
EOF
if [ "$SlidStr" = "Weertman" ]; then
cat >> "${SifFileName}" << EOF
  Weertman Friction Coefficient = Real \$C
  Weertman Exponent = Real \$(1.0/${SlidExp})
  Weertman Linear Velocity = Real 0.001
EOF
elif [ "$SlidStr" = "Coulomb" ]; then
cat >> "${SifFileName}" << EOF
  Friction Law Sliding Coefficient = Real \$C
  Friction Law Post-Peak Exponent = Real 1.0
  Friction Law Maximum Value = Real 0.5
  Friction Law PowerLaw Exponent = Real \$(1.0/${SlidExp})
  Friction Law Linear Velocity = Real 0.001
EOF
fi
cat >> "${SifFileName}" << EOF
  ! Options are 'Last Grounded' (default), 'First Floating' or 'Discontinuous'
    Grounding Line Definition = String "Discontinuous"
  Test Contact Tolerance = real 1.0e-3
  !Non Detachment Inland Distance = Real 5000.0 ! distance from the GL where nodes

  FlowVar 1 = Real 0.0
  FlowVar 1 Condition = Variable GroundedMask
    Real MATC "tx + 0.5"
  V 1 = Real 0.0
  V 1 Condition = Variable GroundedMask
    Real MATC "tx + 0.5"
!
! Shelf conditions
!
  External Pressure = Variable Coordinate 3
     Real Procedure "ElmerIceUSF" "SeaPressure"

  Slip Coefficient 1 = Variable Coordinate 3
     Real Procedure "ElmerIceUSF" "SeaSpring"

  ComputeNormal Condition = Variable GroundedMask
    Real MATC "tx + 0.5"

  Compute Sea Pressure = Logical True
  Compute Sea Spring = Logical True

  Distance = Real 0.0
  Distance Condition = Variable GroundedMask
    Real MATC "tx"
End

!! BC Lateral Ice-Shelf (air or sea contact)
!! BC  Free surface Top
Boundary Condition 6
  Name = "top"
  Target Boundaries = 6
  Body Id = 2
  ComputeNormal = Logical False
  Age = Real 0.0
End
EOF
}
###############################################################################
### This function creates an Elmer sif file for the forward simulation ########
###############################################################################
CreateElmerSifForwardPS() {
				SifFileName=$1
				MeshLayers=${2}
				Rhoi=${3}
				Rhow=${4}
				RateFactor=${5}
				GlenExponent=${6}
				BasalFrictionCoeff=${7}
				IceTemp=${8}
				SMB=${9}
				FluxInit=${10}
				OutputName=${11}
				RestartName=${12}
				NoOfTimeSteps=${13}
				OutputInterval=${14}
				TimeStepSize=${15}
				Alpha=${16}
				G=${17}
				A=${18}
				rho=${19}
                SlidStr=${20}
                SlidExp=${21}
				rm ${SifFileName}.bak
cat > "${SifFileName}.bak" << EOF
!!--------------------------------------------------------!!
!  Island ice rise setup for forwards simulation
!!--------------------------------------------------------!!

check keywords warn
!
! working units are MPa, a, m
!
\$yearinsec = 365.25*24*60*60
\$rhoi = ${Rhoi}/(1.0e6*yearinsec^2)
\$rhow = ${Rhow}/(1.0e6*yearinsec^2)
\$A = ${RateFactor}*yearinsec*1.0e18
\$n = ${GlenExponent}
\$eta = 1.0/(2.0*A)^(1.0/n)
\$gravity = -9.8*yearinsec^2
\$C = ${BasalFrictionCoeff}/(1.0e6*yearinsec^(1.0/n))


Header
  Mesh DB "." "Mesh"
End

Constants
  Water Density = Real \$rhow
  Gas Constant = Real 8.314 !Joule/mol x  K
	Alpha = Real $Alpha
	G = Real $G
	A = Real $A
	rho = Real $rho
  ! For SeaSpring/SeaPressure
End

!---------------------------------------------------
!---------------- SIMULATION -----------------------
!---------------------------------------------------

Simulation
  Coordinate System  = Cartesian 3D
  Simulation Type = transient
  Extruded Mesh Levels = Integer $MeshLayers

  Timestepping Method = "bdf"
  BDF Order = 1
  Timestep Intervals = 50
  Output Intervals = 50
  Timestep Sizes = 0.5

  Initialize Dirichlet Conditions = Logical False
  Steady State Max Iterations = 1
  Steady State Min Iterations = 1

	Restart File="${RestartName}START.result"
	Restart Before Initial Conditions = Logical True
	!Specify name of result file. Used for restarts!!
  Output File = "${OutputName}END.result"
  max output level = 30
End

!---------------------------------------------------
!---------------- BODIES ---------------------------
!---------------------------------------------------

! the ice
Body 1
  Name = "ice"
  Equation = 1
  Body Force = 1
  Material = 1
  Initial Condition = 1
End

! The upper surface
Body 2
  Name= "top free surface"
  Equation = 2
  Material = 1
  Body Force = 2
  Initial Condition = 2
End

! the lower surface
Body 3
  Name= "free surface sea/ice-shelf"
  Equation = 3
  Material = 1
  Body Force = 3
  Initial Condition = 3
End

!---------------------------------------------------
!---------------- INITIAL CONDITIONS ---------------
!---------------------------------------------------

!! for ice
Initial Condition 1
  !Age = Real 0.0
End

!! for top free surface
Initial Condition 2
End

!! for free surface sea/ice-shelf
Initial Condition 3
End

!---------------------------------------------------
!---------------- BODY FORCES ----------------------
!---------------------------------------------------

Body Force 1
  Flow BodyForce 1 = Real 0.0
  Flow BodyForce 2 = Real 0.0
  Flow BodyForce 3 = Real \$gravity
End

!! accumulation flux in m/year
Body Force 2
   Zs Accumulation Flux 1 = Real 0.0e0
   Zs Accumulation Flux 2 = Real 0.0e0 !m/a
   Zs Accumulation Flux 3 = Real $SMB
End

!! no melting/accretion under ice/shelf
Body Force 3
  !Zb Accumulation = Real $BMB
	Zb Accumulation = Variable Time
			Real Procedure "src/USF_BMB" "GetBMB"
End

!---------------------------------------------------
!---------------- MATERIALS ------------------------
!---------------------------------------------------

!! ice material properties in MPa - m - a system
Material 1
  Viscosity Model = String "power law"
  Density = Real \$rhoi
  Viscosity = Real \$eta
  Viscosity Exponent = Real \$1.0/n
  Critical Shear Rate = Real 1.0e-15

  Sea level = Real 0.0

  Glen Enhancement Factor = Real 1.0
! the temperature to switch between the
! two regimes in the flow law
  Limit Temperature = Real -10.0
! In case there is no temperature variable
  Constant Temperature = Real $IceTemp

  Min Zs = Variable "Bottom Zb"
    Real MATC "tx + 10.0"
  Max Zs = Real 1.0e6

  !! Bed condition
  Min Zb = Equals BedBump
  Max Zb = Real 1.0e6
  Cauchy = Logical True
End

!---------------------------------------------------
!---------------- SOLVERS --------------------------
!---------------------------------------------------
Solver 1
  !Exec Solver = Never
  Equation = "MapCoordinate"
  Procedure = "StructuredMeshMapper" "StructuredMeshMapper"

  Active Coordinate = Integer 3
  Mesh Velocity Variable = String "dSdt"
  Mesh Update Variable = String "dS"
  Mesh Velocity First Zero = Logical True

  Top Surface Variable Name = String "Zs"
  Bottom Surface Variable Name = String "Zb"

  Displacement Mode = Logical False
  Correct Surface = Logical True
  Minimum Height = Real 1.0
End

Solver 2
  !Exec Solver = Never
  Equation = "NormalVector"
  Procedure = "ElmerIceSolvers" "ComputeNormalSolver"
  Variable = String "Normal Vector"
  Variable DOFs = 3

  ComputeAll = Logical False
  Optimize Bandwidth = Logical False
End

Solver 3
  !Exec Solver = Never
  Equation = Fw
  Procedure = "ElmerIceSolvers" "GetHydrostaticLoads"
  Variable = Fw[Fwater:3]
  Variable DOFs = 3
End

Solver 4
  Exec Solver = Before Simulation
  Equation = "Velocity Preconditioning"
  Procedure = "VelocityPrecond" "VelocityPrecond"
  Variable = -dofs 3 "V"
  Optimize Bandwidth = Logical True
  Linear System Row Equilibration = Logical True
  Element = "p:1 b:4"
  Bubbles in Global System = False
  Linear System Row Equilibration = Logical True
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStabL
  Linear System Max Iterations = 1000
  Linear System Convergence Tolerance = 1.0e-3
  Linear System Preconditioning = ILU0
  Skip Compute Nonlinear Change = Logical True
  Back Rotate N-T Solution = Logical False
End

Solver 5
  Exec Solver = Before Simulation
  Equation = "Pressure Preconditioning"
  Procedure = "PressurePrecond" "PressurePrecond"
  Variable = -dofs 1 "p" ! that would be default
  Element = "p:1 b:4"
  Bubbles in Global System = False
  Linear System Solver = iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 1000
  Linear System Convergence Tolerance = 1.0e-6
  Linear System Preconditioning = None
  Skip Compute Nonlinear Change = Logical True
  Back Rotate N-T Solution = Logical False
End

Solver 6
  Equation = String "Stokes"
  Procedure = "ParStokes" "StokesSolver"
  Variable = "FlowVar"
	Flow Solver Name = String "FlowVar"
  Variable Dofs = 4
  Element = "p:1 b:4"
	Convective = Logical False
	Block Diagonal A = Logical True
	Use Velocity Laplacian = Logical True
	!-----------------------------------------------
  ! Keywords related to the block preconditioning
  !-----------------------------------------------
	! Option 1
  Block Preconditioning = Logical True
	Linear System Convergence Tolerance = 1e-7
  Linear System Adaptive Tolerance = Logical True
  Linear System Base Tolerance = Real 1.0e-3
  Linear System Relative Tolerance = Real 1.0e-1

  Steady State Convergence Tolerance = 1.0E-03
  Nonlinear System Convergence Tolerance = 1.0E-05
  Nonlinear System Max Iterations = 35
  !Nonlinear System Min Iterations = 10
  Nonlinear System Min Iterations = 5
  Nonlinear System Newton After Iterations = 35
  Nonlinear System Newton After Tolerance =  1.0E-05
	! Option 2
	!  Block Preconditioning = Logical True
	! Linear System GCR Restart = Integer 200
	! Linear System Max Iterations = 200
	! Linear System Convergence Tolerance = 1.0e-6
	! !Linear System Convergence Tolerance = 1.0e-9
	! Nonlinear System Max Iterations = 35
	! Nonlinear System Convergence Tolerance = 1.0e-5
	! !Nonlinear System Convergence Tolerance = 1.0e-7
	! Nonlinear System Newton After Tolerance = 1.0e-5


  Exported Variable 1 = FlowVar Loads[Stress Vector:3 CEQ Residual:1]
  Calculate Loads = Logical True
  Exported Variable 2 = -dofs 1 "dSdt"
  Exported Variable 3 = -dofs 1 "dS"
  Exported Variable 4 = -dofs 1 "BedInit"
  Exported Variable 5 = -dofs 1 "RiseIntoShelf"
  Exported Variable 6 = -dofs 1 "BedBump"
  Exported Variable 7 = -dofs 1 "FluxInit"
  Exported Variable 8 = -dofs 1 "Mesh Velocity"
  Exported Variable 9 = -dofs 1 "Mesh Change"
End

Solver 7
  Equation = String "StressSolver"
  Procedure =  File "ElmerIceSolvers" "ComputeDevStress"
  ! this is just a dummy, hence no output is needed
  Variable = -nooutput "Sij"
  Variable DOFs = 1
  ! the name of the variable containing the flow solution
  !(U,V,W,Pressure)
  Flow Solver Name = String "FlowVar"
  ! no default value anymore for "Stress Variable Name"
  Stress Variable Name = String "Stress"
  !-----------------------------------------------------------------------
   Exported Variable 1 = "Stress"
   Exported Variable 1 DOFs = 6   ! 4 in 2D, 6 in 3D
   Linear System Solver = "Iterative"
   Linear System Iterative Method = "BiCGStab"
   Linear System Max Iterations = 300
   Linear System Convergence Tolerance = 1.0E-09
   Linear System Abort Not Converged = True
   Linear System Preconditioning = "ILU3"
   Linear System Residual Output = 1
End

Solver 8
  !Exec Solver = Never
  Equation = "HeightDepth"
  Procedure = "StructuredProjectToPlane" "StructuredProjectToPlane"
  Active Coordinate = Integer 3
  Dot Product Tolerance = Real 1.0e-3

  Operator 1 = Depth
  Operator 2 = Height
  Variable 3 = Zb
  Operator 3 = Bottom
End

Solver 9
  !Exec Solver = Never
   Equation = "SolveDistance"

   Procedure = "src/DistanceSolveRD" "DistanceSolver1"
   Variable = Distance

   H scale = real 2
   Distance Pseudo DT = Real 100
! Nonlinear System Relaxation Factor = 0.25

   Nonlinear System Max Iterations = 50
   Nonlinear System Convergence Tolerance = 1.0e-5

 ! Linear System Solver = Direct
 ! Linear System Direct Method = UMFPack
   Linear System Solver = "Iterative"
   Linear System Iterative Method = "BiCGStab"
   Linear System Max Iterations = 300
   Linear System Convergence Tolerance = 1.0E-09
   Linear System Abort Not Converged = False
   Linear System Preconditioning = "ILU1"
   Linear System Residual Output = 1
   Steady State Convergence Tolerance = 1.0e-4

   Dummy Distance Computation = Logical False

End

Solver 10
  !Exec Solver = Never
  Equation = "Free Surface Top"
  Procedure =  "./src/MyFreeSurfaceSolver" "FreeSurfaceSolver"
  !Procedure =  "FreeSurfaceSolver" "FreeSurfaceSolver"
  Flow Solver Name = String "FlowVar"
  Flow Loads Name = String "FlowVar Loads"
  Variable = "Zs"
  Variable DOFs =  1
  Exported Variable 1 = "Zs Residual"
  Exported Variable 1 DOFs = 1

  !Before Linsolve = "EliminateDirichlet" "EliminateDirichlet"

  Linear System Solver = Iterative
  !Linear System Direct Method = UMFPACK
  Linear System Max Iterations = 1500
  Linear System Iterative Method = BiCGStab
  Linear System Preconditioning = ILU0
  Linear System Convergence Tolerance = Real 1.0e-6
  Linear System Abort Not Converged = False
  Linear System Residual Output = 1

  Nonlinear System Max Iterations = 100
  Nonlinear System Convergence Tolerance  = 1.0e-5
  Nonlinear System Relaxation Factor = 1.00

  Steady State Convergence Tolerance = 1.0e-03

  Stabilization Method = Stabilized
  Apply Dirichlet = Logical True

  Relaxation Factor = Real 1.0
End

Solver 11
  !Exec Solver = Never
  Equation = "Free Surface Sea/Shelf"
  Procedure =  "FreeSurfaceSolver" "FreeSurfaceSolver"
  Flow Solver Name = String "FlowVar"
  Flow Loads Name = String "FlowVar Loads"
  Variable = "Zb"
  Variable DOFS =  1
  Exported Variable 1 = "Zb Residual"
  Exported Variable 1 DOFs = 1

  Nonlinear Update Exported Variables = Logical True

  Exported Variable 2 = "Zb Accumulation "
  Exported Variable 2 DOFS = 1

  !Before Linsolve = "EliminateDirichlet" "EliminateDirichlet"

  Linear System Solver = Iterative
  Linear System Direct Method = UMFPACK
  Linear System Max Iterations = 1500
  Linear System Iterative Method = BiCGStab
  Linear System Preconditioning = ILU0
  Linear System Convergence Tolerance = Real 1.0e-6
  Linear System Abort Not Converged = False
  Linear System Residual Output = 1

  Nonlinear System Max Iterations = 100
  Nonlinear System Convergence Tolerance  = 1.0e-5
  Nonlinear System Relaxation Factor = 1.00

  Steady State Convergence Tolerance = 1.0e-03

  Stabilization Method = Stabilized
  Apply Dirichlet = Logical True

  Relaxation Factor = Real 1.0
End

Solver 12
  !Exec Solver = Never
  Equation = GroundedMask
  Procedure = "ElmerIceSolvers" "GroundedSolver"
  Variable = GroundedMask
  Variable DOFs = 1

  Toler = Real 1.0e-3
  Bedrock Variable = String "BedBump"
End

Solver 13
  Exec Solver = Never
  Equation = "Age Equation"

  Variable = String "Age"
  Variable DOFs =  1

  Flow Solution Name = String "Flow Solution"
  ! Linear System Solver = Iterative
  ! Linear System Max Iterations = 1000
  ! Linear System Iterative Method = Diagonal
  ! Linear System Preconditioning = NoNe
  ! Linear System Convergence Tolerance = Real 1.0e-6
  ! Linear System Abort Not Converged = False
  ! Linear System Residual Output = 0
  Linear System Solver = Direct
  Linear System Direct Method = mumps
  mumps percentage increase working space = integer 1000

  Procedure = "./src/AgeSolverRD" "AgeSolver"
  Exported Variable 1 = -dofs 1 "age"
End

Solver 14
  Exec Solver = After Saving
  Equation = "result output"
  Procedure = "ResultOutputSolve" "ResultOutputSolver"
  Save Geometry Ids = Logical True ! add this line if you want to access boundaries in Paraview
  Output File Name = File "${OutputName}END"
  Output Format = String vtu
End


!---------------------------------------------------
!---------------- EQUATIONS ------------------------
!---------------------------------------------------

Equation 1
  Active Solvers (9) = 1 2 4 5 6 7 8 13 14
End

Equation 2
  Active Solvers(1) = 10
  Flow Solution Name = String "FlowVar"
  Convection = String Computed
End

Equation 3
  Active Solvers(4) = 3 9 11 12 
  Flow Solution Name = String "FlowVar"
  Convection = String Computed
End
!---------------------------------------------------
!---------------- BOUNDARY CONDITIONS --------------
!---------------------------------------------------

!! Back
Boundary Condition 1
  Name = "back"
  Target Boundaries = 1
  Age = Real 0.0
	FlowVar 1 = Variable FluxInit
		Real Procedure "src/BedrockBump" "IceFluxAtBack"
	FlowVar 2 = Real 0.0
	V 1 = Variable FluxInit
		Real Procedure "src/BedrockBump" "IceFluxAtBack"
	V 2 = Real 0.0
End

Boundary Condition 2
  Name = "Looking Downhill Left"
  Target Boundaries = 2

  FlowVar 2 = Real 0.0
  V 2 = Real 0.0
End

!! BC Lateral Ice-Shelf (air or sea contact)
Boundary Condition 3
  Name = "front"
  Target Boundaries = 3


  External Pressure = Variable Coordinate 3
     Real Procedure "ElmerIceUSF" "SeaPressure"

  Compute Sea Pressure = Logical True
  ComputeNormal = Logical False

End

Boundary Condition 4
  Name = "Looking Downhill Right"
  Target Boundaries = 4

  FlowVar 2 = Real 0.0
  V 2 = Real 0.0
End

Boundary Condition 5
  Name = "bottom"
  Target Boundaries = 5
  Body Id = 3

  Normal-Tangential FlowVar = Logical True
  Normal-Tangential V = Logical True
  Flow Force BC = Logical True

!
! Condition where the bed is stuck
!
  Zb = Equals BedBump
  Zb Condition = Variable GroundedMask
    Real MATC "tx + 0.5"
!
! Bedrock conditions
!
  Slip Coefficient 2 = Variable Coordinate 1
    Real Procedure "src/USF_Contact" "SlidCoef_Contact"
  Slip Coefficient 3 = Variable Coordinate 1
    Real Procedure "src/USF_Contact" "SlidCoef_Contact"
  Flow Solution Name = String "FlowVar"
  Sliding Law = String "${SlidStr}"
EOF
if [ "$SlidStr" = "Weertman" ]; then
cat >> "${SifFileName}.bak" << EOF
  Weertman Friction Coefficient = Real \$C
  Weertman Exponent = Real \$(1.0/${SlidExp})
  Weertman Linear Velocity = Real 0.001
EOF
elif [ "$SlidStr" = "Coulomb" ]; then
cat >> "${SifFileName}.bak" << EOF
  Friction Law Sliding Coefficient = Real \$C
  Friction Law Post-Peak Exponent = Real 1.0
  Friction Law Maximum Value = Real 0.5
  Friction Law PowerLaw Exponent = Real \$(1.0/${SlidExp})
  Friction Law Linear Velocity = Real 0.001
EOF
fi
cat >> "${SifFileName}.bak" << EOF
  ! Options are 'Last Grounded' (default), 'First Floating' or 'Discontinuous'
    Grounding Line Definition = String "Discontinuous"
  Test Contact Tolerance = real 1.0e-3
  !Non Detachment Inland Distance = Real 5000.0 ! distance from the GL where nodes

  FlowVar 1 = Real 0.0
  FlowVar 1 Condition = Variable GroundedMask
    Real MATC "tx + 0.5"
  V 1 = Real 0.0
  V 1 Condition = Variable GroundedMask
    Real MATC "tx + 0.5"
!
! Shelf conditions
!
  External Pressure = Variable Coordinate 3
     Real Procedure "ElmerIceUSF" "SeaPressure"

  Slip Coefficient 1 = Variable Coordinate 3
     Real Procedure "ElmerIceUSF" "SeaSpring"

  ComputeNormal Condition = Variable GroundedMask
    Real MATC "tx + 0.5"

  Compute Sea Pressure = Logical True
  Compute Sea Spring = Logical True

  Distance = Real 0.0
  Distance Condition = Variable GroundedMask
    Real MATC "tx"
End

!! BC Lateral Ice-Shelf (air or sea contact)
!! BC  Free surface Top
Boundary Condition 6
  Name = "top"
  Target Boundaries = 6
  Body Id = 2
  ComputeNormal = Logical False
End
EOF
}

