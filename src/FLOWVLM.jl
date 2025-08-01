"""
  Vortex lattice solver.

  # AUTHORSHIP
    * Author    : Eduardo J Alvarez
    * Email     : Edo.AlvarezR@gmail.com
    * Created   : Jul 2017
    * License   : GNU AFFERO GENERAL PUBLIC LICENSE
"""
module FLOWVLM

# ------------ GENERIC MODULES -------------------------------------------------
import Dates
using PyPlot # Comment PyPlot out if using ProfileView
using Dierckx
using LinearAlgebra: dot, norm, cross


# ------------ FLOW LAB MODULES ------------------------------------------------
import AirfoilPrep
ap = AirfoilPrep

import GeometricTools
gt = GeometricTools

import CCBlade
ccb = CCBlade

# ------------ HEADERS ---------------------------------------------------------
for header_name in ["dt", "solver", "wing", "wingsystem",
                            "postprocessing", "rotor_ccb", "rotor", "tools"]
  include("FLOWVLM_"*header_name*".jl")
end


# ------------ GLOBAL VARIABLES ------------------------------------------------
const pm = 3/4 # Default chord-position of the control point
const pn = 1/4 # Default chord-position of the bound vortex

# Fields that can be calculated (implementation exists)
# FIELDS[i] = [[Dependent fields], field type (scalar/vector)]
# Each new field must be implemented into `calculate_field()`
const FIELDS = Dict(
    "Gamma" =>    [[], "scalar"],         # Vortex strength
    "Vinf"  =>    [[], "vector"],         # Velocity at each CP used for Gamma
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
    ################## EXTRA FIELDS ####################
    "mu"    =>    [[], "scalar"],        # Dipole strength of dragging line
    "Vind"  =>    [[], "vector"],        # Any induced velocity field
    "Vvpm"  =>    [[], "vector"],        # Velocity induced by VPM
    "Vvpm_ApA"  =>    [[], "vector"],        # Velocity induced by VPM
    "Vvpm_AB"   =>    [[], "vector"],        # Velocity induced by VPM
    "Vvpm_BBp"  =>    [[], "vector"],        # Velocity induced by VPM
    "Vkin"  =>    [[], "vector"],        # Kinematic velocity
    "ftot"  =>    [[], "vector"],        # Aerodynamic force (D+L+S) per unit span
    "default-vector"  =>    [[], "vector"],  # Place holder for a vector field
    "default-scalar"  =>    [[], "scalar"],  # Place holder for a scalar field
  )







################################################################################
# WING AND WINGSYSTEM COMMON FUNCTIONS
################################################################################
"Solves the VLM of the Wing or WingSystem"
function solve(wing::AbstractWing, Vinf; t::FWrap=0.0,
                vortexsheet=nothing, extraVinf=nothing, keep_sol=false,
                extraVinfArgs...)

  # Sets Vinf (this forces to recalculate horseshoes)
  setVinf(wing, Vinf; keep_sol=keep_sol)

  # Obtain horseshoes
  HSs = getHorseshoes(wing; t=t, extraVinf=extraVinf, extraVinfArgs...)
  Vinfs = getVinfs(wing; t=t, extraVinf=extraVinf, extraVinfArgs...)
  
  # Calls the solver
  Gammas = VLMSolver.solve(HSs, Vinfs; t=t, vortexsheet=vortexsheet)
                            # extraVinf=extraVinf, extraVinfArgs...)

  _addsolution(wing, "Gamma", Gammas; t=t)
  _addsolution(wing, "Vinf", Vinfs; t=t)
end

"Returns all the horseshoes of the Wing or WingSystem"
function getHorseshoes(wing; t::FWrap=0.0, extraVinf...)
  m = get_m(wing)
  HSs = Array{Any,1}[]
  for i in 1:m
    push!(HSs, getHorseshoe(wing, i; t=t, extraVinf...))
  end
  return HSs
end

"Returns the velocity induced at point X"
function Vind(wing::AbstractWing{TF_design,TF_trajectory}, X; t::FWrap=0.0, ign_col::Bool=false,
                        ign_infvortex::Bool=false, only_infvortex::Bool=false) where {TF_design,TF_trajectory}
  TF_promoted = promote_type(TF_design,TF_trajectory,typeof(t))
  V = zeros(TF_promoted, 3)
  # Adds the velocity induced by each horseshoe
  for i in 1:get_m(wing)
    HS = getHorseshoe(wing, i; t=t)
    V += VLMSolver.V(HS, X; ign_col=ign_col, ign_infvortex=ign_infvortex,
                                                only_infvortex=only_infvortex)
  end
  return V
end

function get_hash(var::String)
  return VLMSolver.HS_hash[var]
end
##### END OF COMMONS ###########################################################


end # END OF MODULE
