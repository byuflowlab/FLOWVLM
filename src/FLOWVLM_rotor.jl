# FLOWVLM module for the modeling of rotors/propellers/wind turbines. This
# module is under development, hence the following dependecies are being
# hardcoded:

# CCBlade https://github.com/byuflowlab/ccblade
ccblade_path = "/home/user/Dropbox/FLOWResearch/FLOWCodes/CCBlade/"
include(ccblade_path*"src/CCBlade.jl")
ccb = CCBlade



################################################################################
# ROTOR CLASS
################################################################################
"""
  `Rotor(CW, r, chord, theta, LE_x, LE_z, B, airfoil)`

Object defining the geometry of a rotor/propeller/wind turbine.

  # Arguments
  * CW::Bool                   : True for clockwise rotation, false for CCW.
  * r::Array{Float64,1}        : Radius position for the following variables.
  * chord::Array{Float64,1}    : Chord length.
  * theta::Array{Float64,1}    : Angle of attack (deg) from the rotor's plane
                                  of rotation.
  * LE_x::Array{Float64,1}     : x-position of leading edge.
  * LE_z::Array{Float64,1}     : z-position of leading edge (height from plane
                                  of rotation).
  * B::Int64                   : Number of blades.

  # Optional Arguments
  * airfoils::Array{Tuple{Float64, airfoilprep.Polar},1} : 2D airfoil properties
                                  along blade in the form [ (r_i, Polar_i) ]
                                  with Polar_i describes the airfoil at i-th
                                  radial position r_i (both the airfoil geometry
                                  in Polar_i and r_i must be normalized). At
                                  least root (r=0) and tip (r=1) must be given
                                  so all positions in between can be
                                  extrapolated. This properties are only used
                                  when calling CCBlade and for generating good
                                  loking visuals; ignore if only solving the VLM.

NOTE: r here is the radial position after precone is included in the geometry,
hence the need of explicitely declaring LE_z.

  # PROPERTIES
  * `sol` : Contains solution fields specific for Rotor types. They are formated
            as sol[field_name] = Dict(
                            "field_name" => output_field_name,
                            "field_type" => "scalar" or "vector",
                            "field_data" => data
                            )
            where `data` is an array data[i] = [val1, val2, ...] containing
            this field values (scalar or vector) of all control points in the
            i-th blade.

<!-- NOTE TO SELF: r is the y-direction on a wing, hence, remember to build the
               blade from root in the direction of positive y. -->
"""
type Rotor

  # Initialization variables (USER INPUT)
  CW::Bool                      # True for clockwise rotation
  r::Array{Float64,1}           # Radius position for the following variables
  chord::Array{Float64,1}       # Chord length
  theta::Array{Float64,1}       # Angle of attack (deg) from the rotor's axis
  LE_x::Array{Float64,1}        # x-position of leading edge
  LE_z::Array{Float64,1}        # z-position of leading edge (Height from plane
                                #                                  of rotation)
  B::Int64                      # Number of blades
  # Optional inputs
  airfoils::Array{Tuple{Float64,ap.Polar},1} # 2D airfoil properties along blade

  # Properties
  hubR::Float64                 # Hub radius
  rotorR::Float64               # Rotor radius
  m::Int64                      # Number of control points (per blade)
  sol::Dict{String,Any}         # Solution fields for CCBlade (not FLOWVLM)

  # Data storage
  _wingsystem::WingSystem       # Rotor assembly
  _r::Array{Float64,1}          # Radius of each control point (on one blade)
  _chord::Array{Float64,1}      # Chord length at each control point
  _theta::Array{Float64,1}      # Angle of attack (deg) at each control point
  _LE_x::Array{Float64,1}
  _LE_z::Array{Float64,1}
  _polars::Array{ap.Polar, 1}   # Polar object at each control point (with x,y
                                #  containing the exact geometric airfoil)
  _polarroot::ap.Polar          # Polar at the root
  _polartip::ap.Polar           # Polar at the tip

  Rotor(
          CW, r, chord, theta, LE_x, LE_z, B,
          airfoils=Tuple{Float64, ap.Polar}[],
          hubR=r[1], rotorR=r[end],
            m=0, sol=Dict(),
          _wingsystem=WingSystem(),
            _r=Float64[], _chord=Float64[], _theta=Float64[],
            _LE_x=Float64[], _LE_z=Float64[],
            _polars=ap.Polar[],
              _polarroot=ap.dummy_polar(), _polartip=ap.dummy_polar()
        ) = new(
          CW, r, chord, theta, LE_x, LE_z, B,
          airfoils,
          hubR, rotorR,
            m, sol,
          _wingsystem,
            _r, _chord, _theta,
            _LE_x, _LE_z,
            _polars, _polarroot, _polartip
        )
end

"Initializes the geometry of the rotor, discretizing each blade into n lattices"
function initialize(self::Rotor, n::Int64; r_lat::Float64=1.0,
                          central=false, refinement=[])
  # Checks for arguments consistency
  _check(self)

  # Flag for calculating airfoils
  rfl_flag = size(self.airfoils)[1]!=0

  # Generates blade
  blade, r, chord, theta, LE_x, LE_z = _generate_blade(self, n; r=r_lat,
                                        central=central, refinement=refinement)
  self._r, self._chord, self._theta = r, chord, theta
  self._LE_x, self._LE_z = LE_x, LE_z
  self.m = get_m(blade)

  # Generates airfoil properties at all control points of this blade
  if rfl_flag; _calc_airfoils(self, n, r_lat, central, refinement); end;

  # Generates full rotor
  init_angle = pi/2
  d_angle = 2*pi/self.B
  blades_Oaxis = self.CW ? [0 0 -1; 0 1 0; 1 0 0] : [0 0 1; 0 1 0; -1 0 0]
  for i in 1:self.B
    this_blade = i==1 ? blade : copy(blade)
    this_angle = init_angle + (i-1)*d_angle

    # Sets the blade in the rotor coordinate system, and rotates it
    # this_Oaxis = [0 0 -1;
    #               -1 0 0;
    #               0 1 0]*[
    #               1 0 0;
    #               0 cos(this_angle) sin(this_angle);
    #               0 -sin(this_angle) cos(this_angle)]
    this_Oaxis = blades_Oaxis*[
                  1 0 0;
                  0 cos(this_angle) sin(this_angle);
                  0 -sin(this_angle) cos(this_angle)]
    setcoordsystem(this_blade, [0.0,0,0], this_Oaxis)

    # Adds it to the rotor
    addwing(self._wingsystem, "Blade$(i)", this_blade)
  end
end

"Sets Vinf(X,t) as the incoming freestream of this rotor"
function setVinf(self::Rotor, Vinf)
  _reset(self)
  setVinf(self._wingsystem, Vinf)
end

"Saves the rotor in VTK legacy format"
function save(self::Rotor, filename::String; args...)
  save(self._wingsystem, filename; args...)

  if size(self.airfoils)[1]!=0
    save_loft(self, filename; args...)
  end
end

function setcoordsystem(self::Rotor, args...)
  setcoordsystem(self._wingsystem, args...)
end

"Returns total number of lattices on each blade"
function get_mBlade(self::Rotor)
  return self.m
end

"Returns the requested blade"
function get_blade(self::Rotor, blade_i::Int64)
  return get_wing(self._wingsystem, blade_i)
end

"Returns total number of lattices in the rotor"
function get_m(self::Rotor)
  return get_m(self._wingsystem)
end

"Returns the m-th control point of the system"
function getControlPoint(self::Rotor, m::Int64)
  return getControlPoint(self._wingsystem, m)
end

"Returns the m-th horseshoe of the system in the global coordinate system"
function getHorseshoe(self::Rotor, m::Int64; t::Float64=0.0, extraVinf...)
  getHorseshoe(self._wingsystem, m; t=t, extraVinf...)
end

"Saves the lofted surface of the blade"
function save_loft(self::Rotor, filename::String; path="", num=nothing, args...)
  suf = "loft"

  CP_index = []       # Stores the CP index in order to hash the points
  lines = []          # Contour lines of cross sections

  # Iterates over each airfoil creating cross sections
  for (i,polar) in enumerate(self._polars)
    theta = pi/180*self._theta[i]   # Angle of attack

    # Actual airfoil contour
    x, y = self._chord[i]*polar.x, self._chord[i]*polar.y
    # Flips it if clockwise rotating rotor
    if self.CW; y = -y; end;
    # Reformats x,y into point
    points = [ [x[i], y[i], 0.0] for i in 1:size(x)[1] ]
    # Rotates the contour in the right angle of attack
    Oaxis = [ cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1]
    points = vtk.countertransform(points, inv(Oaxis), zeros(3))

    # Position of leading edge in FLOVLM blade's c.s.
    # Airfoil's x-axis = FLOWVLM blade's x-axis
    # Airfoil's y-axis = FLOWVLM blade's z-axis
    O = [self._LE_x[i], self._r[i], self._LE_z[i]]
    # Reformats the contour into FLOWVLM blade's c.s.
    points = [ O+[p[1], p[3], p[2]] for p in points]

    # Adds this airfoil
    push!(lines, points)

    # Adds the CP index as point data for all points in this airfoil
    push!(CP_index, [i for p in points])
  end

  # Generates vtk cells from cross sections
  sections = [ [(1.0, 1, 1.0, false)] for i in 1:size(lines)[1]-1]
  points, vtk_cells, CP_index = vtk.multilines2vtkmulticells(lines, sections;
                                                      point_datas=CP_index)

  # Generates each blade
  for i in 1:self.B # Iterates over blades
    this_blade = self._wingsystem.wings[i]

    # Transforms points from FLOWVLM blade's c.s. to global c.s.
    this_points = Array{Float64,1}[ vtk.countertransform(p, this_blade.invOaxis,
                                    this_blade.O) for p in points]

    # Formats the point data for generateVTK
    data = []
    # Control point indices
    push!(data, Dict(
                "field_name" => "ControlPoint_Index",
                "field_type" => "scalar",
                "field_data" => CP_index
                )
          )
    # Stored fields
    for (field_name, field) in self.sol # Iterates over fields
      data_points = [] # Field data associated to each geometric point

      for (j,p) in enumerate(this_points) # Iterates over geometric points
        CP_i = Int(CP_index[j])  # Hashes the control point index of this geom point
        push!(data_points, field["field_data"][i][CP_i])
      end

      push!(data, Dict(
                  "field_name" => field["field_name"],
                  "field_type" => field["field_type"],
                  "field_data" => data_points
                  )
            )
    end

    # Generates the vtk file
    this_name = filename*"_"*self._wingsystem.wing_names[i]*"_"*suf
    vtk.generateVTK(this_name, this_points; cells=vtk_cells, point_data=data,
                                path=path, num=num)
  end
end



"""
Returns a CCBlade's rotor object corresponding to the i-th blade. The r position
used to generate the object are the control points, `Rhub` is the minimum radius
given when initializing the FLOWVLM Rotor object, and `Rtip` is the maximum. The
precone given to CCBlade is 0 since precone it is expected to be already
accounted for in the geometry given by the user when initializing the FLOWVLM
Rotor.

It applies 3D corrections to the airfoil polars given for initializing the
FLOWVLM Rotor, hence the RPMs for those corrections must be specified.

WARNING: This function will use the Inflow field, hence make sure they are
updated.

NOTE: In the current implementation it assumes a constant Reynolds number on
each airfoil throughout operation as captured in the polar used for initializing
the FLOWVLM Rotor. The implementation of a varying Reynolds will be left for
future development as needed.
"""
function FLOWVLM2CCBlade(self::Rotor, RPM, blade_i::Int64)
  # ERROR CASES
  if size(self.airfoils)[1]<2
    error("Airfoil data not found when generating CCBlade Rotor.")
  elseif size(self._polars)[1]==0
    error("Control point polars haven't been calculated yet."*
              " Run `_calc_airfoils()` before calling this function.")
  elseif !("CCBInflow" in keys(self.sol))
    error("CCBInflow field not found. "*
              "Call `calc_inflow()` before calling this function")
  end

  Rhub = self.hubR
  Rtip = self.rotorR
  precone = 0.0
  inflows = self.sol["CCBInflow"]["field_data"][blade_i]

  # Prepares airfoil polars
  af = ccb.AirfoilData[]
  for (i,polar) in enumerate(self._polars)

    # 3D corrections
    r_over_R = self._r[i] / Rtip
    c_over_r = self._chord[i] / self._r[i]
    #   NOTE: Here I'm taking the freestream to be the absolute value CCBlade's
    #   x-component of inflow. This may cause problems when the flow is reversed
    tsr = (2*pi*RPM/60 * Rtip) / abs(inflows[i][1])
    this_polar = ap.correction3D(polar, r_over_R, c_over_r, tsr)

    # 360 extrapolation
    CDmax = 1.3
    this_polar = ap.extrapolate(this_polar, CDmax)

    # Makes sure the polar is injective for easing the spline
    this_polar = ap.injective(this_polar)

    # Converts to CCBlade's AirfoilData object
    alpha, cl = ap.get_cl(this_polar)
    _, cd = ap.get_cd(this_polar)
    ccb_polar = ccb.af_from_data(alpha*pi/180, cl, cd)

    push!(af, ccb_polar)
  end

  rotor = ccb.Rotor(self._r, self._chord,
                    self._theta*pi/180, af,
                    Rhub, Rtip, self.B, precone)

  return rotor
end


##### CALCULATION OF SOULTION FIELDS############################################
"Receives the freestream velocity function V(x,t) and the current RPM of the
rotor, and it calculates the inflow velocity field that each control point
sees in the global coordinate system"
function calc_inflow(self::Rotor, Vinf, RPM; t::Float64=0.0)
  omega = 2*pi*RPM/60

  data_Vtots = Array{Array{Float64,1}}[]     # Inflow in the global c.s.
  data_Vccbs = Array{Array{Float64,1}}[]     # Inflow in CCBlade's c.s.

  for (i,blade) in enumerate(self._wingsystem.wings) # Iterates over blades
    Vtots = Array{Float64,1}[]
    Vccbs = Array{Float64,1}[]

    for j in 1:get_m(blade) # Iterates over control points
      CP = getControlPoint(blade, j)

      # Freestream velocity in global c.s.
      this_Vinf = Vinf(CP, t)

      # Velocity due to rotation in FLOWVLM blade's c.s.
      this_Vrot = [omega*self._r[j], 0.0, 0.0]
      # Velocity due to rotation in global c.s.
      this_Vrot = countertransform(this_Vrot, blade.invOaxis, zeros(3))

      this_Vtot = this_Vinf + this_Vrot
      push!(Vtots, this_Vtot)
      push!(Vccbs, _global2ccblade(blade, this_Vtot, self.CW))
    end

    push!(data_Vtots, Vtots)
    push!(data_Vccbs, Vccbs)
  end

  # Adds solution fields
  field_tots = Dict(
              "field_name" => "GlobInflow",
              "field_type" => "vector",
              "field_data" => data_Vtots
              )
  self.sol[field_tots["field_name"]] = field_tots

  field_ccbs = Dict(
              "field_name" => "CCBInflow",
              "field_type" => "vector",
              "field_data" => data_Vccbs
              )
  self.sol[field_ccbs["field_name"]] = field_ccbs
end

"Calculates the distributed loads from CCBlade"
function calc_distributedloads(self::Rotor, Vinf, RPM, rho::Float64;
                                    t::Float64=0.0)
  data = Array{Array{Float64,1}}[]

  # Calculates inflows
  calc_inflow(self, Vinf, RPM; t=t)

  # Distributed load on each blade
  for blade_i in 1:self.B

    # Formats the inflow as CCBlade's Inflow
    inflow = self.sol["CCBInflow"]["field_data"][blade_i]
    inflow_x = [V[1] for V in inflow]
    inflow_y = [V[2] for V in inflow]
    ccbinflow = ccb.Inflow(inflow_x, inflow_y, rho)

    # Generates CCBlade Rotor object
    ccbrotor = FLOWVLM2CCBlade(self, RPM, blade_i)

    # Calls CCBlade
    # NOTE TO SELF: Forces normal and tangential to the plane of rotation
    Np, Tp, uvec, vvec = ccb.distributedloads(ccbrotor, ccbinflow, false)

    # Convert forces from CCBlade's c.s. to global c.s.
    ccb_Fs = [ [Np[i], Tp[i], 0.0] for i in 1:get_mBlade(self) ]
    Fs = [ _ccblade2global( get_blade(self, blade_i), ccb_F, self.CW;
                          translate=false) for ccb_F in ccb_Fs ]
    # Stores the field
    push!(data, Fs)
  end

  # Adds solution fields
  field = Dict(
              "field_name" => "DistributedLoad",
              "field_type" => "vector",
              "field_data" => data
              )
  self.sol[field["field_name"]] = field
end

##### INTERNAL FUNCTIONS #######################################################
"Checks for consistency in internal variables"
function _check(self::Rotor)
  # Matching dimensions
  nr = size(self.r)[1]
  for (label, arr) in [("chord",self.chord), ("theta",self.theta),
                        ("LE_x",self.LE_x), ("LE_z",self.LE_z)]
    if size(arr)[1] != nr
      error("Invalid dimensions in array $label."*
                                    " Expected $nr, found $(size(arr)[1]).")
    end
  end
  # r from root to tip
  prev_ri = 0
  for ri in self.r
    if ri<prev_ri
      error("Array r must be in increasing order and non-negative.")
    end
    prev_ri = ri
  end
  # Airfoil data
  if size(self.airfoils)[1]!=0
    # Enough polars
    if size(self.airfoils)[1]==1
      error("Two or more airfoils expected. "*
                  "Received $(size(self.airfoils)[1])")
    elseif self.airfoils[1][1]!=0.0
      error("Missing root airfoil.")
    elseif self.airfoils[end][1]!=1.0
      error("Missing tip airfoil.")
    end
  end
end

function _reset(self::Rotor; verbose=false, keep_Vinf=false)
  if verbose; println("Resetting Rotor"); end;
  _reset(self._wingsystem; verbose=verbose, keep_Vinf=keep_Vinf)
end

"Generates the blade and discretizes it into lattices"
function _generate_blade(self::Rotor, n::Int64; r::Float64=1.0,
                          central=false, refinement=[])

  # Splines
  spline_k = min(size(self.r)[1]-1, 3)
  spline_bc = "error"
  spline_s = 0.001
  _spl_chord = Dierckx.Spline1D(self.r, self.chord; k=spline_k,s=spline_s)
  _spl_theta = Dierckx.Spline1D(self.r, -self.theta; k=spline_k,s=spline_s)
  _spl_LE_x = Dierckx.Spline1D(self.r, self.LE_x; k=spline_k,s=spline_s)
  _spl_LE_z = Dierckx.Spline1D(self.r, self.LE_z; k=spline_k,s=spline_s)
  spl_chord(x) = Dierckx.evaluate(_spl_chord, x)
  spl_theta(x) = (-1)^(self.CW==false)*Dierckx.evaluate(_spl_theta, x)
  spl_LE_x(x) = Dierckx.evaluate(_spl_LE_x, x)
  spl_LE_z(x) = Dierckx.evaluate(_spl_LE_z, x)

  # Outputs
  out_r = Float64[]         # Radial position of each control point
  out_chord = Float64[]     # Chord at each control point
  out_theta = Float64[]     # Angle of attach at each control point
  out_LE_x = Float64[]
  out_LE_z = Float64[]


  # Precalulations of complex refinement
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
  end


  # Initializes the blade
  blade = Wing(spl_LE_x(self.r[1]), self.r[1], spl_LE_z(self.r[1]),
                spl_chord(self.r[1]), spl_theta(self.r[1]))

  # Discretizes the leading edge in n lattices
  l = self.rotorR - self.hubR # Lenght of the leading edge from root to tip
  cumlen = 0 # Cumulative length of the leading edge already walked
  prev_r = self.hubR
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
      end

    end
    cumlen += len
    this_r = self.hubR + (cumlen/l)*(self.rotorR-self.hubR)
    # println("i=$i\tr=$this_r\tchord=$(spl_chord(this_r))")

    # Adds this lattice
    addchord(blade, spl_LE_x(this_r), this_r, spl_LE_z(this_r),
              spl_chord(this_r), spl_theta(this_r), 1)

    # Properties at control point
    CP_r = (this_r+prev_r)/2
    push!(out_r, CP_r)
    push!(out_chord, spl_chord(CP_r))
    push!(out_theta, spl_theta(CP_r))
    push!(out_LE_x, spl_LE_x(CP_r))
    push!(out_LE_z, spl_LE_z(CP_r))

    prev_r = this_r
  end

  # Verifies correct discretization
  if abs((cumlen-l)/l)>0.0000000001
    error("Critical logic error! cumlen!=l ($cumlen!=$l)")
  end

  return blade, out_r, out_chord, out_theta, out_LE_x, out_LE_z
end

"Calculates the airfoils at each control point"
function _calc_airfoils(self::Rotor, n::Int64, r::Float64,
                        central, refinement; rediscretize::Bool=true,
                        rfl_n_lower::Int64=15, rfl_n_upper::Int64=15,
                        rfl_r::Float64=14.0, rfl_central::Bool=true)

  # Erases previous polars
  self._polars = ap.Polar[]

  # Normalized position of each interval boundary
  aux_f(x) = x
  if size(refinement)[1]!=0
    # Precalulations of complex refinement
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
    points = vtk.multidiscretize(aux_f, 0, 1,
          [(sec[1], ns[i], sec[3], false) for (i,sec) in enumerate(refinement)])
  else
    points = vtk.discretize(aux_f, 0, 1, n, r; central=central)
  end

  # Normalized position of control points
  norm_CPs = [ (points[i]+points[i-1])/2 for i in 2:size(points)[1] ]

  # ------ Iterates over CPs calculating the polar of the airfoil at that CP ---

  prev_contour, next_contour = nothing, nothing # Airfoil contours
  prev_i, next_i = nothing, 1 # Indices of airfoils previous and next to cur CP

  # Iterates over each control point
  for this_CP in norm_CPs

    # Finds the next airfoil closest to it
    while this_CP > self.airfoils[next_i][1]
      prev_i, prev_contour = next_i, next_contour
      next_i += 1
      next_contour = nothing
    end

    # Rediscretizes any contours as needed
    if prev_contour==nothing
      if rediscretize
        prev_contour = _rediscretize_airfoil(self.airfoils[prev_i][2].x,
                                              self.airfoils[prev_i][2].y,
                                              rfl_n_lower, rfl_n_upper, rfl_r,
                                              rfl_central)
      else
        prev_contour = (self.airfoils[prev_i][2].x, self.airfoils[prev_i][2].y)
      end
    end
    if next_contour==nothing
      if rediscretize
        next_contour = _rediscretize_airfoil(self.airfoils[next_i][2].x,
                                              self.airfoils[next_i][2].y,
                                              rfl_n_lower, rfl_n_upper, rfl_r,
                                              rfl_central)
      else
        next_contour = (self.airfoils[next_i][2].x, self.airfoils[next_i][2].y)
      end
    end

    # Weight of the next airfoil vs the previous
    weight = (this_CP-self.airfoils[prev_i][1]
                ) / (self.airfoils[next_i][1]-self.airfoils[prev_i][1])

    # Blends Cd and Cl curves
    prev_pypolar = self.airfoils[prev_i][2].pyPolar
    next_pypolar = self.airfoils[next_i][2].pyPolar
    blended_pypolar = prev_pypolar[:blend](next_pypolar, weight)

    # Blends airfoil geometry
    blended_x = weight*next_contour[1] + (1-weight)*prev_contour[1]
    blended_y = weight*next_contour[2] + (1-weight)*prev_contour[2]

    # Blended Polar object
    blended_polar = ap._pyPolar2Polar(blended_pypolar, blended_x, blended_y)

    push!(self._polars, blended_polar) # Stores CP airfoils in self._polars
  end

  # Stores root and tip airfoils
  if rediscretize
    root_x, root_y = _rediscretize_airfoil(self.airfoils[1][2].x,
                                            self.airfoils[1][2].y,
                                            rfl_n_lower, rfl_n_upper, rfl_r,
                                            rfl_central)
    tip_x, tip_y = _rediscretize_airfoil(self.airfoils[end][2].x,
                                            self.airfoils[end][2].y,
                                            rfl_n_lower, rfl_n_upper, rfl_r,
                                            rfl_central)
    self._polarroot = ap._pyPolar2Polar(self.airfoils[1][2].pyPolar,
                                            root_x, root_y)
    self._polartip = ap._pyPolar2Polar(self.airfoils[end][2].pyPolar,
                                            tip_x, tip_y)
  else
    self._polarroot = self.airfoils[1][2]
    self._polartip = self.airfoils[end][2]
  end

end

function _rediscretize_airfoil(x, y, n_lower::Int64, n_upper::Int64, r::Float64,
                                central::Bool)
  # Separate upper and lower sides to make the contour injective in x
  upper, lower = ap.splitcontour(x, y)

  # Parameterize both sides independently
  fun_upper = vtk.parameterize(upper[1], upper[2], zeros(upper[1]); inj_var=1)
  fun_lower = vtk.parameterize(lower[1], lower[2], zeros(lower[1]); inj_var=1)

  # New discretization for both surfaces
  upper_points = vtk.discretize(fun_upper, 0, 1, n_upper, r[1]; central=central)
  lower_points = vtk.discretize(fun_lower, 0, 1, n_lower, r[1]; central=central)

  # Put both surfaces back together from TE over the top and from LE over the bottom.
  reverse!(upper_points)                        # Trailing edge over the top
  new_x = [point[1] for point in upper_points]
  new_y = [point[2] for point in upper_points]  # Leading edge over the bottom
  new_x = vcat(new_x, [point[1] for point in lower_points])
  new_y = vcat(new_y, [point[2] for point in lower_points])

  return new_x, new_y
end

"""Receives a vector in the global coordinate system and transforms it into
CCBlade's coordinate system relative to `blade`. NOTE: This function only
rotates `V` into the new axis without translating it unless otherwise indicated.

NOTE TO SELF:
CCblade's blade x-axis = FLOWVLM Rotor's blade z-axis
CCblade's blade y-axis = FLOWVLM Rotor's blade x-axis
CCblade's blade z-axis = FLOWVLM Rotor's blade y-axis"""
function _global2ccblade(blade::Wing, V::Array{Float64,1}, CW::Bool;
                                                        translate::Bool=false)
  # V in FLOWVLM Rotor's blade c.s.
  V_vlm = transform(V, blade.Oaxis, translate ? blade.O : zeros(3))
  # V in CCBlade's c.s.
  V_ccb = typeof(V)([ (-1)^(CW)*V_vlm[3],  V_vlm[1], V_vlm[2]])

  return V_ccb
end


"""Receives a vector in CCBlade's coordinate system relative to `blade` and
transforms it into the global coordinate system. NOTE: This function only
rotates `V` into the new axis without translating it unless otherwise indicated.

NOTE TO SELF:
FLOWVLM Rotor's blade x-axis = CCblade's blade y-axis
FLOWVLM Rotor's blade y-axis = CCblade's blade z-axis
FLOWVLM Rotor's blade z-axis = CCblade's blade x-axis"""
function _ccblade2global(blade::Wing, V::Array{Float64,1}, CW::Bool;
                                                        translate::Bool=false)
  # V in FLOWVLM Rotor's blade c.s.
  V_vlm = typeof(V)([ V[2],  V[3], (-1)^(CW)*V[1] ])
  # V in global c.s.
  V_glob = countertransform(V_vlm, blade.invOaxis, translate ? blade.O : zeros(3))

  return V_glob
end

##### END OF ROTOR CLASS #######################################################












#
