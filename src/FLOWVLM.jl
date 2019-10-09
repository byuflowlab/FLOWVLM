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
using PyPlot # Comment PyPlot out if using ProfileView
using Dierckx


# ------------ FLOW LAB MODULES ------------------------------------------------
# The following modules are under development, hence imports here are hardcoded

# Airfoil processing https://github.com/EdoAlvarezR/airfoil
airfoil_path = "/home/edoalvar/Dropbox/FLOWResearch/FLOWCodes/airfoil/"
# push!(LOAD_PATH, joinpath(airfoil_path,"src/"))
# using airfoilprep
include(airfoil_path*"src/airfoilprep.jl")
ap = airfoilprep

# VTKtools https://github.com/byuflowlab/VTKtools.jl
# vtktools_path = "/home/edoalvar/Dropbox/FLOWResearch/FLOWCodes/GeometricTools/"
# push!(LOAD_PATH, joinpath(vtktools_path,"src/"))
# using VTKtools
# include(vtktools_path*"src/GeometricTools.jl")
import GeometricTools
vtk = GeometricTools


# ------------ HEADERS ---------------------------------------------------------
for header_name in ["dt", "solver", "wing", "wingsystem",
                                    "tools", "postprocessing"]
  include("FLOWVLM_"*header_name*".jl")
end
include("utilities.jl")

# try # Rotor module is under developments
  include("FLOWVLM_rotor.jl")
# catch e
#   warn("FLOWVLM_rotor.jl module failed to load: $e")
# end


# ------------ GLOBAL VARIABLES ------------------------------------------------
const pm = 3/4 # Default chord-position of the control point
const pn = 1/4 # Default chord-position of the bound vortex

const def_airfoil = ap.data_path*"oval00.txt"  # Default airfoil shape

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
    "Vvpm"  =>    [[], "vector"],        # Velocity induced by VPM
    "Vkin"  =>    [[], "vector"],        # Kinematic velocity
    "ftot"  =>    [[], "vector"],        # Aerodynamic force (D+L+S) per unit span
  )







################################################################################
# WING AND WINGSYSTEM COMMON FUNCTIONS
################################################################################
"Solves the VLM of the Wing or WingSystem"
function solve(wing, Vinf; t::FWrap=0.0,
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
function Vind(wing, X; t::FWrap=0.0, ign_col::Bool=false, ign_infvortex::Bool=false)
  V = zeros(3)
  # Adds the velocity induced by each horseshoe
  for i in 1:get_m(wing)
    HS = getHorseshoe(wing, i; t=t)
    V += VLMSolver.V(HS, X; ign_col=ign_col, ign_infvortex=ign_infvortex)
  end
  return V
end

function get_hash(var::String)
  return VLMSolver.HS_hash[var]
end
##### END OF COMMONS ###########################################################


end # END OF MODULE
