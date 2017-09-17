# @Author: Eduardo Alvarez <user>
# @Date:   2017-07-20T12:37:14-06:00
# @Email:  Edo.AlvarezR@gmail.com
# @Last modified by:   user
# @Last modified time: 2017-07-20T12:37:14-06:00
# @Comments: Vortex lattice solver

module FLOWVLM

using PyPlot # Comment PyPlot out if using ProfileView
include("FLOWVLM_solver.jl")
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
# WING CLASS
################################################################################
"""
  `Wing(leftxl, lefty, leftz, leftchord, leftchordtwist)`

Wing constructor that automatically discretizes the surface into lattices. The
wing must we built from left to right. Chords are parallel to the zx-plane,
leading edge in the direction of the -xaxis and trailing in the direction of the
+xaxis (hence, the span is aligned with the y-axis).

<!-- The i-th control point is located at (_xm[i], _ym[i]), with left bound -->
<!-- vortex (_xn[i], _yn[i]) and right bound vortex (_xn[i+1], _yn[i+1]) -->

  # Arguments
  *  `leftxl`         : x-position of leading edge of the left tip.
  *  `lefty`          : y-position of left tip.
  *  `leftzl          : z-position of leading edge of the left tip.
  *  `leftchord`      : Chord of the left tip.
  *  `leftchordtwist` : Twist left tip's chord in degrees.

  # Example
  `julia julia> wing = Wing(0.0, 0.0, 0.0, 10.0, 3.0);`
"""
type Wing

  # Initialization variables (USER INPUT)
  leftxl::Float64       # x-position of leading edge of the left tip
  lefty::Float64        # y-position of left tip
  leftzl::Float64       # z-position of leading edge of the left tip
  leftchord::Float64    # Chord of the left tip
  leftchordtwist::Float64    # Angle of the chord of the left tip in degrees

  # Properties
  m::Int64                      # Number of lattices
  O::Array{Float64,1}           # Origin of local reference frame
  Oaxis::Array{Float64,2}       # Unit vectors of the local reference frame
  invOaxis::Array{Float64,2}    # Inverse unit vectors
  Vinf::Any                     # Vinf function used in current solution

  # Data storage
  ## Solved fields
  sol::Dict{String,Any}         # Dictionary storing every solved field
  ## Discretized wing geometry
  _xlwingdcr::typeof(Float64[]) # x-position of leading edge
  _xtwingdcr::typeof(Float64[]) # x-position of trailing edge
  _ywingdcr::typeof(Float64[])  # y-position of the chord
  _zlwingdcr::typeof(Float64[]) # z-position of leading edge
  _ztwingdcr::typeof(Float64[]) # z-position of trailing edge
  ## VLM domain
  _xm::typeof(Float64[]) # x-position of the control point
  _ym::typeof(Float64[]) # y-position of the control point
  _zm::typeof(Float64[]) # z-position of the control point
  _xn::typeof(Float64[]) # x-position of the bound vortex
  _yn::typeof(Float64[]) # y-position of the bound vortex
  _zn::typeof(Float64[]) # z-position of the bound vortex
  ## Calculation data
  _HSs::Any              # Horseshoes

  Wing(leftxl, lefty, leftzl, leftchord, leftchordtwist,
                  m=0,
                    O=[0.0,0.0,0.0],
                    Oaxis=[1.0 0 0; 0 1 0; 0 0 1],
                    invOaxis=[1.0 0 0; 0 1 0; 0 0 1],
                    Vinf=nothing,
                  sol=Dict(),
                  _xlwingdcr=[leftxl],
                    _xtwingdcr=[leftxl+leftchord*cos(leftchordtwist*pi/180)],
                    _ywingdcr=[lefty],
                    _zlwingdcr=[leftzl],
                    _ztwingdcr=[leftzl-leftchord*sin(leftchordtwist*pi/180)],
                  _xm=Float64[], _ym=Float64[], _zm=Float64[],
                  _xn=[leftxl+pn*leftchord*cos(leftchordtwist*pi/180)],
                    _yn=[lefty],
                    _zn=[leftzl-pn*leftchord*sin(leftchordtwist*pi/180)],
                  _HSs=nothing
              ) = new(leftxl, lefty, leftzl, leftchord, leftchordtwist,
                      m, O, Oaxis, invOaxis, Vinf,
                      sol,
                      _xlwingdcr, _xtwingdcr, _ywingdcr, _zlwingdcr, _ztwingdcr,
                      _xm, _ym, _zm, _xn, _yn, _zn,
                      _HSs
                      )
end



"""
  `addchord(wing, x, y, z, c, twist, n, r=1.0)`

Adds a new chord to the wing and creates n lattices in the new section. Wing
must be build from left to right.

  # Arguments
  *     x       : x-position leading edge of the chord.
  *     y       : y-position of the chord.
  *     z       : z-position leading edge of the chord.
  *     c       : Chord length.
  *     twist   : Twist of the chord in degrees.
  *     n::Int64: Number of lattices in the new section.

  # Optional arguments
  *     r       : Ratio between lengths of first and last lattices.
  *     central : Give it true to take the length ratio between the lattice
                  midway and first and last. Give it a number between 0 and 1
                  to define the position of the reference midway.
  *     refinement : Use this option for more complex refinements. It
                  receives an array `[sec1, sec2, ...]` with sections of
                  refinements in the format `sec=[c, n, r]`, with `c` the length
                  of this section (sum of all c = 1), `n` the ratio of lattices
                  in this section, and `r` the increment ratio. If this option
                  is used, it will ignore arguments `r` and `central`.

  # Examples
    `julia> wing = Wing(0.0, 0.0, 0.0, 10.0, 3.0);`
    `julia> addchord(wing, 2.5, 10.0, 5.0, 5.0, 0.0, 10);`
"""
function addchord(self::Wing,
                  x::Float64, y::Float64, z::Float64,
                  c::Float64, twist::Float64,
                  n::Int64; r::Float64=1.0, central=false, refinement=[])
  # ERROR CASES
  if c <= 0
    error("Invalid chord length (c <= 0)")
  elseif n <= 0
    error("Invalid number of lattices (n <= 0)")
  elseif r <= 0
    error("Invalid expansion ratio (r <= 0)")
  end

  # RESETS IF EXISTING SOLUTION
  if "Gamma" in keys(self.sol)
    _reset(self)
  end

  _twist = twist*pi/180

  # Precaclulations of complex refinement
  if size(refinement)[1]!=0
    nsecs = size(refinement)[1]
    ntot = sum([refinement[i][2] for i in 1:nsecs])
    ctot = sum([refinement[i][1] for i in 1:nsecs])
    ns = [] # Number of lattices in each section
    for i in 1:nsecs
      if i==nsecs
        push!(ns, n-sum(ns))
      else
        push!(ns, floor(n*refinement[i][2]/ntot))
      end
    end
    # println("$nsecs\t$ntot\t$ctot\t$ns")
  end

  # Left boundary of the new section
  Ll = [self._xlwingdcr[end], self._ywingdcr[end], self._zlwingdcr[end]]
  Lt = [self._xtwingdcr[end], self._ywingdcr[end], self._ztwingdcr[end]]
  # Right boundary of the section
  Rl = [x, y, z]
  Rt = [x + c*cos(_twist), y, z - c*sin(_twist)]

  # Discretizes the section in n lattices
  l = sqrt(dot( Rl-Ll , Rl-Ll )) # Lenght of the new section's leading edge
  cumlen = 0 # Cumulative length of the leading edge already walked
  # println("New section - n=$n\tr=$r\tcentral=$central")
  for i in 1:n

    # Wing discretization
    if size(refinement)[1]!=0 # Complex refinement
      sec_i = 1 # Current section
      for j in 1:nsecs
        if i>sum([ns[k] for k in 1:j])
          sec_i+=1
        else
          break
        end
      end
      _prev_ns = sec_i==1 ? 0 : sum(ns[1:sec_i-1])
      _n = ns[sec_i]
      _r = refinement[sec_i][3]
      _l = l*(refinement[sec_i][1]/ctot)
      _i = i-_prev_ns
      p = _l/( (_n*(_n-1)/2)*(_r+1)/(_r-1) )
      d1 = p*(_n-1)/(_r-1)
      len = d1 + p*(_i-1)

    elseif r==1.0 # i.e., uniform discretization
      len = l/n # This lattice's leading edge length


    else # Linear increment (see notebook entry 20170519)
      # Case of no central expansion
      if central==false
        p = l/( (n*(n-1)/2)*(r+1)/(r-1) )
        d1 = p*(n-1)/(r-1)
        len = d1 + p*(i-1)
      # Case of central expansion
      else
        _central = central==true ? 0.5 : central
        # Left of the center
        # if cumlen/l < _central
        if i<=floor(n*_central)
          _l = l*_central
          _n = floor(n*_central)
          _r = r
          _i = i
        # Right of the center
        else
          _l = l*(1-_central)
          _n = n-floor(n*_central)
          _r = 1/r
          _i = i-floor(n*_central)
        end
        p = _l/( (_n*(_n-1)/2)*(_r+1)/(_r-1) )
        d1 = p*(_n-1)/(_r-1)
        len = d1 + p*(_i-1)
        # println("$i\t$_r\t$cumlen\t_n=$_n\t_i=$_i")
      end

    end
    cumlen += len
    THISl = Ll + (cumlen/l)*(Rl-Ll)
    THISt = Lt + (cumlen/l)*(Rt-Lt)

    push!(self._xlwingdcr, THISl[1])
    push!(self._ywingdcr, THISl[2])
    push!(self._zlwingdcr, THISl[3])
    push!(self._xtwingdcr, THISt[1])
    push!(self._ztwingdcr, THISt[3])

    # Bound vortex
    Cn = THISt - THISl
    N = THISl + pn*Cn
    push!(self._xn, N[1])
    push!(self._yn, N[2])
    push!(self._zn, N[3])

    # Control point
    Ml = Ll + ((cumlen-len/2)/l)*(Rl-Ll)
    Mt = Lt + ((cumlen-len/2)/l)*(Rt-Lt)
    Cm = Mt - Ml
    M = Ml + pm*Cm
    push!(self._xm, M[1])
    push!(self._ym, M[2])
    push!(self._zm, M[3])
  end

  # Verifies correct discretization
  if abs((cumlen-l)/l)>0.0000000001
    error("Critical logic error! cumlen!=l ($cumlen!=$l)")
  end

  # Updates the constructor
  self.m += n
  _reset(self; keep_Vinf=true)
end

"""
  `setcoordsystem(self, O, Oaxis; check=true)`

Redefines the local coordinate system of the wing, where `O` is the new origin
and `Oaxis` is the matrix [i; j; k] of unit vectors
"""
function setcoordsystem(self::Wing, O::Array{Float64,1},
                            Oaxis::Array{Float64,2};
                            check=true)

  if check; check_coord_sys(Oaxis); end;

  self.O = O
  self.Oaxis = Oaxis
  self.invOaxis = inv(Oaxis)
  _reset(self; keep_Vinf=true)
end


function setcoordsystem(self::Wing, O::Array{Float64,1},
                            Oaxis::Array{Array{Float64,1},1};
                            check=true)
  dims = 3
  M = zeros(dims, dims)
  for i in 1:dims
    M[i, :] = Oaxis[i]
  end
  setcoordsystem(self, O, M; check=check)
end

"Sets Vinf(X,t) as the freestream of this wing"
function setVinf(self::Wing, Vinf)
  _reset(self)
  self.Vinf = Vinf
end

"Returns the m-th control point"
function getControlPoint(self::Wing, m::Int64)
  # Local coordinate system
  CP = [self._xm[m], self._ym[m], self._zm[m]]
  # Global coordinate system
  CP = countertransform(CP, self.invOaxis, self.O)
  return CP
end

"Returns the m-th horseshoe in the global coordinate system"
function getHorseshoe(self::Wing, m::Int64; t::Float64=0.0)
  # ERROR CASES
  if m>self.m || m<=0
    error("Invalid m (m>self.m or m<=0)")
  elseif false==("Gamma" in keys(self.sol)) && self.Vinf==nothing
    error("Freestream hasn't been define yet, please call function set_Vinf()")
  end

  # Calculates horseshoes if not available
  if self._HSs==nothing
    _calculateHSs(self; t=t)
  end

  return self._HSs[m]
end

"Returns leading-edge coordinates of the n-th chord"
function getLE(self::Wing, n::Int64)
  # Local coordinate system
  LE = [self._xlwingdcr[n], self._ywingdcr[n], self._zlwingdcr[n]]
  # Global coordinate system
  LE = countertransform(LE, self.invOaxis, self.O)
  return LE
end

"Returns trailing-edge coordinates of the n-th chord"
function getTE(self::Wing, n::Int64)
  # Local coordinate system
  TE = [self._xtwingdcr[n], self._ywingdcr[n], self._ztwingdcr[n]]
  # Global coordinate system
  TE = countertransform(TE, self.invOaxis, self.O)
  return TE
end

"Returns total number of lattices in the wing"
function get_m(self::Wing)
  return self.m
end


##### INTERNAL FUNCTIONS #######################################################
function _reset(self::Wing; verbose=false, keep_Vinf=false)
  if verbose; println("Resetting wing"); end;
  self.sol = Dict()
  if keep_Vinf==false; self.Vinf = nothing; end;
  self._HSs = nothing
end

function _addsolution(self::Wing, field_name::String, sol_field; t::Float64=0.0)
  self.sol[field_name] = sol_field
  if field_name=="Gamma"
    _calculateHSs(self; t=t)
  end
end

function _calculateHSs(self::Wing; t::Float64=0.0)
  HSs = Array{Any,1}[]
  for i in 1:get_m(self)
    # Horseshoe geometry
    ## Points in the local coordinate system
    Ap = [self._xtwingdcr[i], self._ywingdcr[i], self._ztwingdcr[i]]
    A = [self._xn[i], self._yn[i], self._zn[i]]
    B = [self._xn[i+1], self._yn[i+1], self._zn[i+1]]
    Bp = [self._xtwingdcr[i+1], self._ywingdcr[i+1], self._ztwingdcr[i+1]]
    ## Points in the global coordinate system
    Ap, A, B, Bp = countertransform([Ap, A, B, Bp], self.invOaxis, self.O)

    # Control point
    CP = getControlPoint(self, i)

    # Direction of semi-infinite vortices
    infDA = self.Vinf(Ap,t); infDA = infDA/norm(infDA)
    infDB = self.Vinf(Bp,t); infDB = infDB/norm(infDB)

    # Circulation
    Gamma = "Gamma" in keys(self.sol) ? self.sol["Gamma"][i] : nothing

    HS = [Ap, A, B, Bp, CP, infDA, infDB, Gamma]
    push!(HSs, HS)
  end
  self._HSs = HSs
end
##### END OF CLASS #############################################################



################################################################################
# WING AND WINGSYSTEM COMMON FUNCTIONS
################################################################################
"Solves the VLM of the Wing or WingSystem"
function solve(wing, Vinf; t::Float64=0.0)
  setVinf(wing, Vinf)
  HSs = getHorseshoes(wing; t=t)
  Gammas = VLMSolver.solve(HSs, Vinf; t=t)
  _addsolution(wing, "Gamma", Gammas; t=t)
end

"Returns all the horseshoes of the Wing or WingSystem"
function getHorseshoes(wing; t::Float64=0.0)
  m = get_m(wing)
  HSs = Array{Any,1}[]
  for i in 1:m
    push!(HSs, getHorseshoe(wing, i; t=t))
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

include("FLOWVLM_wingsystem.jl")
include("FLOWVLM_tools.jl")
include("FLOWVLM_postprocessing.jl")

end # END OF MODULE
