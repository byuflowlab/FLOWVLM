
################################################################################
# WING CLASS
################################################################################
abstract type AbstractWing{TF_design,TF_trajectory} end

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
mutable struct Wing{TF_design<:FWrap,TF_trajectory} <: AbstractWing{TF_design,TF_trajectory}

  # Initialization variables (USER INPUT)
  leftxl::TF_design                 # x-position of leading edge of the left tip
  lefty::TF_design                  # y-position of left tip
  leftzl::TF_design                 # z-position of leading edge of the left tip
  leftchord::TF_design              # Chord of the left tip
  leftchordtwist::TF_design         # Angle of the chord of the left tip in degrees

  # Properties
  m::IWrap                      # Number of lattices
  O::Vector{TF_trajectory}                   # Origin of local reference frame
  Oaxis::Matrix{TF_trajectory}                 # Unit vectors of the local reference frame
  invOaxis::Matrix{TF_trajectory}              # Inverse unit vectors
  Vinf::Union{Nothing,Function}                     # Vinf function used in current solution

  # Data storage
  ## Solved fields
  sol::Dict{String,Any}                  # Dictionary storing every solved field
  ## Discretized wing geometry
  _xlwingdcr::Vector{TF_design}          # x-position of leading edge
  _xtwingdcr::Vector{TF_design}          # x-position of trailing edge
  _ywingdcr::Vector{TF_design}           # y-position of the chord
  _zlwingdcr::Vector{TF_design}          # z-position of leading edge
  _ztwingdcr::Vector{TF_design}          # z-position of trailing edge
  ## VLM domain
  _xm::Vector{TF_design}                 # x-position of the control point
  _ym::Vector{TF_design}                 # y-position of the control point
  _zm::Vector{TF_design}                 # z-position of the control point
  _xn::Vector{TF_design}                 # x-position of the bound vortex
  _yn::Vector{TF_design}                 # y-position of the bound vortex
  _zn::Vector{TF_design}                 # z-position of the bound vortex
  ## Calculation data
  _HSs::Union{Nothing,Vector{Vector{Any}}}              # Horseshoes
end

function Wing(leftxl::TF_design, lefty::TF_design, leftzl::TF_design, leftchord::TF_design, leftchordtwist::TF_design, TF_trajectory=Float64; 
  m=0,
  O=[0.0,0.0,0.0],
  Oaxis=[1.0 0 0; 0 1 0; 0 0 1],
  invOaxis=[1.0 0 0; 0 1 0; 0 0 1],
  Vinf=nothing,
  sol=Dict(),
  _xlwingdcr::Vector{TF_design}=[leftxl],
  _xtwingdcr::Vector{TF_design}=[leftxl+leftchord*cos(leftchordtwist*pi/180)],
  _ywingdcr::Vector{TF_design}=[lefty],
  _zlwingdcr::Vector{TF_design}=[leftzl],
  _ztwingdcr::Vector{TF_design}=[leftzl-leftchord*sin(leftchordtwist*pi/180)],
  _xm=TF_design[], _ym=TF_design[], _zm=TF_design[],
  _xn::Vector{TF_design}=[leftxl+pn*leftchord*cos(leftchordtwist*pi/180)],
  _yn::Vector{TF_design}=[lefty],
  _zn::Vector{TF_design}=[leftzl-pn*leftchord*sin(leftchordtwist*pi/180)],
  _HSs=nothing
) where {TF_design}
  O=convert(Vector{TF_trajectory}, O)
  Oaxis=convert(Matrix{TF_trajectory}, Oaxis)
  invOaxis=convert(Matrix{TF_trajectory}, invOaxis)
  return Wing{TF_design, TF_trajectory}(leftxl, lefty, leftzl, leftchord, leftchordtwist,
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
                  x::FWrap, y::FWrap, z::FWrap,
                  c::FWrap, twist::FWrap,
                  n::IWrap; r::FWrap=1.0, central=false, refinement=[])
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
      if i==nsecs && nsecs!=1
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
    setcoordsystem(wing::Wing, O::Vector, Oaxis::Matrix)

Redefines the local coordinate system of the wing, where `O` is the new origin
and `Oaxis` is the matrix of unit vectors
"""
function setcoordsystem(self::Wing, O::Vector{<:FWrap},
                            Oaxis::Matrix{<:FWrap};
                            check=true)

  if check; check_coord_sys(Oaxis); end;

  self.O .= O
  self.Oaxis .= Oaxis
  self.invOaxis .= inv(Oaxis)
  _reset(self; keep_Vinf=true)
end


function setcoordsystem(self::Wing, O::Vector{<:FWrap},
                            Oaxis::Array{T,1} where {T<:AbstractArray};
                            check=true)
  dims = 3
  M = zeros(eltype(Oaxis), dims, dims)
  for i in 1:dims
    M[i, :] = Oaxis[i]
  end
  setcoordsystem(self, O, M; check=check)
end

"Sets Vinf(X,t) as the freestream of this wing"
function setVinf(self::Wing, Vinf; keep_sol=false)
  _reset(self; keep_sol=keep_sol)
  self.Vinf = Vinf
end

"Returns the m-th control point"
function getControlPoint(self::Wing, m::IWrap)
  # Local coordinate system
  CP = [self._xm[m], self._ym[m], self._zm[m]]
  # Global coordinate system
  CP = countertransform(CP, self.invOaxis, self.O)
  return CP
end

"Returns the undisturbed freestream at each control point, or at the horseshoe
point indicated as `target`."
function getVinfs(self::Wing{TF_design,TF_trajectory}; t::FWrap=0.0, target="CP",
                              extraVinf=nothing, extraVinfArgs...) where {TF_design,TF_trajectory}
  if !(target in keys(VLMSolver.HS_hash))
    error("Logic error! Invalid target $target.")
  end
  t_i = VLMSolver.HS_hash[target]

  # Calculates Vinf at each control point
  TF_promoted = promote_type(TF_design,TF_trajectory)
  Vinfs = Vector{Vector{TF_promoted}}(undef,get_m(self))
  for i in 1:get_m(self)
    if target=="CP"
      T = getControlPoint(self, i)      # Targeted point
    else
      T = getHorseshoe(self, i; t=t, extraVinf=extraVinf, extraVinfArgs...)[t_i]
    end

    this_Vinf = self.Vinf(T, t)
    if extraVinf!=nothing; this_Vinf += extraVinf(i, t; extraVinfArgs..., wing=self); end;

    Vinfs[i] = this_Vinf
  end

  return Vinfs
end

"Returns the m-th horseshoe in the global coordinate system"
function getHorseshoe(self::Wing, m::IWrap; t::FWrap=0.0, extraVinf...)
  # ERROR CASES
  if m>self.m || m<=0
    error("Invalid m (m>self.m or m<=0)")
  elseif false==("Gamma" in keys(self.sol)) && self.Vinf==nothing
    error("Freestream hasn't been define yet, please call function setVinf()")
  end

  # Calculates horseshoes if not available
  if self._HSs==nothing
    _calculateHSs(self; t=t, extraVinf...)
  end

  return self._HSs[m]
end

"Returns leading-edge coordinates of the n-th chord"
function getLE(self::Wing, n::IWrap)
  # Local coordinate system
  LE = [self._xlwingdcr[n], self._ywingdcr[n], self._zlwingdcr[n]]
  # Global coordinate system
  LE = countertransform(LE, self.invOaxis, self.O)
  return LE
end

"Returns trailing-edge coordinates of the n-th chord"
function getTE(self::Wing, n::IWrap)
  # Local coordinate system
  TE = [self._xtwingdcr[n], self._ywingdcr[n], self._ztwingdcr[n]]
  # Global coordinate system
  TE = countertransform(TE, self.invOaxis, self.O)
  return TE
end

"""
    get_m(wing::Wing)

Returns the number of horseshoes in the wing
"""
function get_m(self::Wing)
  return self.m
end

"Returns a deep copy of this wing"
function copy(self::Wing)
  return deepcopy(self)
end

##### INTERNAL FUNCTIONS #######################################################
function _reset(self::Wing; verbose=false, keep_Vinf=false, keep_sol=false)
  if verbose; println("Resetting wing"); end;
  if keep_sol==false
      self.sol = Dict()
  else
      self.sol = Dict([entry for entry in self.sol if entry[1]!="Gamma"])
  end
  if keep_Vinf==false; self.Vinf = nothing; end;
  self._HSs = nothing
end

function _addsolution(self::Wing, field_name::String, sol_field; t::FWrap=0.0)
  self.sol[field_name] = sol_field
  if field_name=="Gamma"
    # _calculateHSs(self; t=t)
    for (i, gamma) in enumerate(sol_field)
      self._HSs[i][8] = gamma
    end
  end
end

function _calculateHSs(self::Wing; t::FWrap=0.0, extraVinf=nothing, extraVinfArgs...)
  HSs = Vector{Any}[]
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
    infDA = self.Vinf(Ap,t)
    infDB = self.Vinf(Bp,t)
    # Extra freestream
    if extraVinf!=nothing
      this_extraVinf = extraVinf(i, t; extraVinfArgs..., wing=self)
      infDA += this_extraVinf
      infDB += this_extraVinf
    end
    infDA = infDA/norm(infDA)
    infDB = infDB/norm(infDB)

    # Circulation
    Gamma = "Gamma" in keys(self.sol) ? self.sol["Gamma"][i] : nothing

    HS = [Ap, A, B, Bp, CP, infDA, infDB, Gamma]
    push!(HSs, HS)
  end
  self._HSs = HSs
end

function Base.deepcopy_internal(x::Wing, stackdict::IdDict)
    if haskey(stackdict, x)
        return stackdict[x]
    end

    y = Wing(Base.deepcopy_internal(x.leftxl, stackdict),
             Base.deepcopy_internal(x.lefty, stackdict),
             Base.deepcopy_internal(x.leftzl, stackdict),
             Base.deepcopy_internal(x.leftchord, stackdict),
             Base.deepcopy_internal(x.leftchordtwist, stackdict),
             Base.deepcopy_internal(x.m, stackdict),
             Base.deepcopy_internal(x.O, stackdict),
             Base.deepcopy_internal(x.Oaxis, stackdict),
             Base.deepcopy_internal(x.invOaxis, stackdict),
             x.Vinf,
             Base.deepcopy_internal(x.sol, stackdict),
             Base.deepcopy_internal(x._xlwingdcr, stackdict),
             Base.deepcopy_internal(x._xtwingdcr, stackdict),
             Base.deepcopy_internal(x._ywingdcr, stackdict),
             Base.deepcopy_internal(x._zlwingdcr, stackdict),
             Base.deepcopy_internal(x._ztwingdcr, stackdict),
             Base.deepcopy_internal(x._xm, stackdict),
             Base.deepcopy_internal(x._ym, stackdict),
             Base.deepcopy_internal(x._zm, stackdict),
             Base.deepcopy_internal(x._xn, stackdict),
             Base.deepcopy_internal(x._yn, stackdict),
             Base.deepcopy_internal(x._zn, stackdict),
             Base.deepcopy_internal(x._HSs, stackdict))

    stackdict[x] = y
    return y
end
##### END OF CLASS #############################################################
