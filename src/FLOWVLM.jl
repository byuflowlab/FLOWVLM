# @Author: Eduardo Alvarez <user>
# @Date:   2017-07-20T12:37:14-06:00
# @Email:  Edo.AlvarezR@gmail.com
# @Last modified by:   user
# @Last modified time: 2017-07-20T12:37:14-06:00
# @Comments: Vortex lattice solver

module FLOWVLM

using PyPlot # Comment PyPlot out if using ProfileView
include("FLOWVLM_solver.jl")
include("FLOWVLM_wing.jl")
include("FLOWVLM_wingsystem.jl")
include("FLOWVLM_rotor.jl")
include("FLOWVLM_tools.jl")
include("FLOWVLM_postprocessing.jl")
include("utilities.jl")

const pm = 3/4 # Default chord-position of the control point
const pn = 1/4 # Default chord-position of the bound vortex


# Fields that can be calculated (implementation exists)
# FIELDS[i] = [[Dependent fields], field type (scalar/vector)]
# Each new field must be implemented into `calculate_field()`
const FIELDS = Dict(
    "Gamma" =>    [[], "scalar"],         # Vortex strength
    ################## LIFT AND SIDEWASH ####################
    "Ftot"  =>    [["Gamma"], "vector"],     # Aerodynamic force (D+L+S)
    "D"     =>    [["Gamma"], "vector"],     # Drag
    "L"     =>    [["Gamma"], "vector"],     # Lift
    "S"     =>    [["Gamma"], "vector"],     # Sideslip force
    "CFtot" =>    [["Gamma"], "vector"],     # COEFFICIENTS PER PANEL
    "CD"    =>    [["Gamma"], "vector"],     #
    "CL"    =>    [["Gamma"], "vector"],     #
    "CS"    =>    [["Gamma"], "vector"],     #
    "Cftot/CFtot" =>    [["CFtot"], "scalar"], # NORMALIZED UNIT-SPAN
    "Cd/CD" =>    [["CFtot"], "scalar"],     #      COEFFICIENTS PER PANEL
    "Cl/CL" =>    [["CFtot"], "scalar"],     #
    "Cs/CS" =>    [["CFtot"], "scalar"],     #
    # ################# INDUCED DRAG ###########################
    # "Dind"  =>    [["Gamma"], "vector"],  # Induced drag
    # "CDind" =>    [["Gamma"], "vector"],  # Induced drag coefficient
    # ################## MOMENTS ###############################
    "A"     =>    [[], "scalar"],         # Area of each panel
    "Mtot"  =>    [["Ftot"], "vector"],  # Total moment (-M_L + M_M - M_N)
    "M_L"   =>    [["Ftot"], "vector"],  # Rolling moment
    "M_M"   =>    [["Ftot"], "vector"],  # Pitching moment
    "M_N"   =>    [["Ftot"], "vector"],  # Yawing moment
    "CMtot" =>    [["Mtot"], "vector"],  # COEFFICIENTS PER PANEL
    "CM_L"  =>    [["Mtot"], "vector"],  #
    "CM_M"  =>    [["Mtot"], "vector"],  #
    "CM_N"  =>    [["Mtot"], "vector"],  #
  )




################################################################################
# WING AND WINGSYSTEM COMMON FUNCTIONS
################################################################################
"Solves the VLM of the Wing or WingSystem"
function solve(wing, Vinf; t::Float64=0.0,
                vortexsheet=nothing, extraVinf=nothing, extraVinfArgs...)
  setVinf(wing, Vinf)
  HSs = getHorseshoes(wing; t=t, extraVinf=extraVinf, extraVinfArgs...)
  Gammas = VLMSolver.solve(HSs, Vinf; t=t, vortexsheet=vortexsheet,
                            extraVinf=extraVinf, extraVinfArgs...)
  _addsolution(wing, "Gamma", Gammas; t=t)
end

"Returns all the horseshoes of the Wing or WingSystem"
function getHorseshoes(wing; t::Float64=0.0, extraVinf...)
  m = get_m(wing)
  HSs = Array{Any,1}[]
  for i in 1:m
    push!(HSs, getHorseshoe(wing, i; t=t, extraVinf...))
  end
  return HSs
end

"Returns the velocity induced at point X"
function Vind(wing, X; t::Float64=0.0, ign_col::Bool=false)
  V = zeros(3)
  # Adds the velocity induced by each horseshoe
  for i in 1:get_m(wing)
    HS = getHorseshoe(wing, i; t=t)
    V += VLMSolver.V(HS, X; ign_col=ign_col)
  end
  return V
end
##### END OF COMMONS ###########################################################


end # END OF MODULE
