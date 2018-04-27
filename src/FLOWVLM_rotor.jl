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

Object defining the geometry of a rotor/propeller/wind turbine. This class
behaves as an extension of the WingSystem class, hence all functions of
WingSystem can be applied to a Rotor object.

  # Arguments
  * CW::Bool                   : True for clockwise rotation, false for CCW.
  * r::Array{Float64,1}        : Radius position for the following variables.
  * chord::Array{Float64,1}    : Chord length.
  * theta::Array{Float64,1}    : Angle of attack (deg) from the rotor's plane
                                  of rotation.
  * LE_x::Array{Float64,1}     : x-position of leading edge (positive is ahead
                                  of radial axis relative to rotation).
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
  r::FArrWrap                   # Radius position for the following variables
  chord::FArrWrap               # Chord length
  theta::FArrWrap               # Angle of attack (deg) from the rotor's axis
  LE_x::FArrWrap                # x-position of leading edge
  LE_z::FArrWrap                # z-position of leading edge (Height from plane
                                #                                  of rotation)
  B::IWrap                      # Number of blades
  # Optional inputs
  airfoils::Array{Tuple{FWrap,ap.Polar},1} # 2D airfoil properties along blade

  # Properties
  RPM::Any                      # Current revs per minute
  hubR::FWrap                   # Hub radius
  rotorR::FWrap                 # Rotor radius
  m::IWrap                      # Number of control points (per blade)
  sol::Dict{String,Any}         # Solution fields for CCBlade (not FLOWVLM)

  # Data storage
  _wingsystem::WingSystem       # Rotor assembly
  _r::FArrWrap                  # Radius of each control point (on one blade)
  _chord::FArrWrap              # Chord length at each control point
  _theta::FArrWrap              # Angle of attack (deg) at each control point
  _LE_x::FArrWrap
  _LE_z::FArrWrap
  _polars::Array{ap.Polar, 1}   # Polar object at each control point (with x,y
                                #  containing the exact geometric airfoil)
  _polarroot::ap.Polar          # Polar at the root
  _polartip::ap.Polar           # Polar at the tip

  Rotor(
          CW, r, chord, theta, LE_x, LE_z, B,
          airfoils=Tuple{FWrap, ap.Polar}[],
          RPM=nothing,
            hubR=r[1], rotorR=r[end],
            m=0, sol=Dict(),
          _wingsystem=WingSystem(),
            _r=FWrap[], _chord=FWrap[], _theta=FWrap[],
            _LE_x=FWrap[], _LE_z=FWrap[],
            _polars=ap.Polar[],
              _polarroot=ap.dummy_polar(), _polartip=ap.dummy_polar()
        ) = new(
          CW, r, chord, theta, LE_x, LE_z, B,
          airfoils,
          RPM,
            hubR, rotorR,
            m, sol,
          _wingsystem,
            _r, _chord, _theta,
            _LE_x, _LE_z,
            _polars, _polarroot, _polartip
        )
end

"Initializes the geometry of the rotor, discretizing each blade into n lattices"
function initialize(self::Rotor, n::IWrap; r_lat::FWrap=1.0,
                          central=false, refinement=[], verif=false, rfl_args...)
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

  # Verifies correct lattice and blade elements discretization
  if verif; _verif_discr(self, blade, r, chord, theta, LE_x, LE_z); end;

  # Generates airfoil properties at all control points of this blade
  if rfl_flag; _calc_airfoils(self, n, r_lat, central, refinement; rfl_args...); end;

  # ------------ Generates full rotor -----------------
  # Default blade c.s. relative to rotor c.s.
  blades_Oaxis = self.CW ? [0 -1 0; 0 0 1; -1.0 0 0] : [0 1 0; 0 0 1; 1.0 0 0]
  init_angle = 0.0
  d_angle = 2*pi/self.B
  for i in 1:self.B
    this_blade = i==1 ? blade : copy(blade)
    this_angle = init_angle + (i-1)*d_angle # Azumithal angle of this blade

    # Sets the blade in the rotor coordinate system, and rotates it
    this_Oaxis = [ cos(this_angle) sin(this_angle) 0;
                  -sin(this_angle) cos(this_angle) 0;
                   0 0 1]*blades_Oaxis
    setcoordsystem(this_blade, [0.0,0,0], this_Oaxis)

    # Adds it to the rotor
    addwing(self, "Blade$(i)", this_blade; force=true)
  end

  # Sets the rotor in the global coordinate system
  rotor_Oaxis = [-1 0 0; 0 -1 0; 0 0 1.0]
  setcoordsystem(self._wingsystem, [0.0,0,0], rotor_Oaxis)
end

"Solves for the Gamma field (circulation) on each blade using CCBlade. It also
includes the fields Ftot, L, D, and S.

If include_comps==true it stores CCBlade-calculated normal and tangential forces
in the Rotor."
function solvefromCCBlade(self::Rotor, Vinf, RPM, rho::FWrap; t::FWrap=0.0,
                            include_comps::Bool=false, return_performance::Bool=false,
                            Vref=nothing, sound_spd=nothing)

  setVinf(self, Vinf)
  setRPM(self, RPM)

  # (Calls a HS to make sure they have been calculated)
  _ = getHorseshoe(self, 1)

  if sound_spd==nothing
    warn("No sound speed has been provided. No Mach corrections will be applied.")
  end

  # Calculates distributed load from CCBlade
  prfrmnc = calc_distributedloads(self, Vinf, RPM, rho; t=t,
                                        include_comps=include_comps,
                                        return_performance=return_performance,
                                        Vref=Vref, sound_spd=sound_spd)

  # Decomposes load into aerodynamic forces and calculates circulation
  gamma = calc_aerodynamicforces(self, rho)

  new_gamma = FWrap[]
  new_Ftot = FArrWrap[]
  new_L = FArrWrap[]
  new_D = FArrWrap[]
  new_S = FArrWrap[]

  # Formats solution fields as a FLOWVLM solution
  for i in 1:self.B  # Iterates over blades

    for j in 1:get_mBlade(self) # Iterates over lattices on blade
      push!(new_gamma, gamma[i][j])
      push!(new_Ftot, self.sol["DistributedLoad"]["field_data"][i][j])
      push!(new_L, self.sol["Lift"]["field_data"][i][j])
      push!(new_D, self.sol["Drag"]["field_data"][i][j])
      push!(new_S, self.sol["RadialForce"]["field_data"][i][j])
    end

  end

  # Adds the fields as FLOWVLM solutions
  _addsolution(self._wingsystem, "Gamma", new_gamma; t=t)
  _addsolution(self._wingsystem, "Ftot", new_Ftot; t=t)
  _addsolution(self._wingsystem, "L", new_L; t=t)
  _addsolution(self._wingsystem, "D", new_D; t=t)
  _addsolution(self._wingsystem, "S", new_S; t=t)

  return prfrmnc
end

"Sets Vinf(X,t) as the incoming freestream of this rotor"
function setVinf(self::Rotor, Vinf)
  _reset(self)
  setVinf(self._wingsystem, Vinf)
end

"Sets `RPM` as the revolutions per minutes of this rotor"
function setRPM(self::Rotor, RPM)
  _reset(self; keep_Vinf=true)
  _resetRotor(self)
  self.RPM = RPM
end

"Saves the rotor in VTK legacy format"
function save(self::Rotor, filename::String; addtiproot=true, args...)
  _ = getHorseshoe(self, 1) # Makes sure the wake is calculated right

  save(self._wingsystem, filename; args...)

  if size(self.airfoils)[1]!=0
    save_loft(self, filename; addtiproot=addtiproot, args...)
  end
end

"Sets a coordinate system for the rotor. If the user is calling this function,
give `user=true`, otherwise it won't do the automatic translation to blade c.s."
function setcoordsystem(self::Rotor, O::FArrWrap,
                            Oaxis::FMWrap; user=false, args...)
  if user
    setcoordsystem(self._wingsystem, O, Oaxis*[-1 0 0; 0 -1 0; 0 0 1.0],args...)
  else
    setcoordsystem(self._wingsystem, O, Oaxis ,args...)
  end

  _resetRotor(self; keep_RPM=true)
end

"Rotates the rotor `degs` degrees in the direction of rotation"
function rotate(self::Rotor, degs::FWrap)
  rotOaxis = vtk.rotation_matrix(0.0, 0.0, (-1)^!self.CW*degs)
  newOaxis = rotOaxis*self._wingsystem.Oaxis
  # setcoordsystem(self._wingsystem, self._wingsystem.O, newOaxis)
  setcoordsystem(self, self._wingsystem.O, newOaxis)
end


"Returns the undisturbed freestream at each control point"
function getVinfs(self::Rotor; t::FWrap=0.0, target="CP",
                              extraVinf=nothing, extraVinfArgs...)
  if !(target in keys(VLMSolver.HS_hash))
    error("Logic error! Invalid target $target.")
  end

  _ = getHorseshoe(self, 1; t=t, extraVinf=extraVinf, extraVinfArgs...)

  Vinfs = FArrWrap[]
  for i in 1:self.B
    blade_Vinfs = _calc_inflow(get_blade(self, i), get_RPM(self), t;
                                                                target=target)
    for V in blade_Vinfs
      push!(Vinfs, V)
    end
  end

  # Adds any extra terms
  if extraVinf!=nothing
    for i in 1:self.B
      blade = get_blade(self, i)
      for j in 1:get_m(blade)
        Vinfs[i][j] += extraVinf(j, t; extraVinfArgs..., wing=blade)
      end
    end
  end

  return Vinfs
end

function get_RPM(self::Rotor)
  if self.RPM==nothing
    error("RPM not defined yet."*
          " Call function `setRPM()` before calling this function.")
  end
  return self.RPM
end

"Returns total number of lattices on each blade"
function get_mBlade(self::Rotor)
  return self.m
end

"Returns the requested blade"
function get_blade(self::Rotor, blade_i::IWrap)
  return get_wing(self, blade_i)
end

"Returns total number of lattices in the rotor"
function get_m(self::Rotor)
  return get_m(self._wingsystem)
end

"Returns the m-th control point of the system"
function getControlPoint(self::Rotor, m::IWrap)
  return getControlPoint(self._wingsystem, m)
end

"Returns the m-th horseshoe of the system in the global coordinate system"
function getHorseshoe(self::Rotor, m::IWrap; t::FWrap=0.0, extraVinf...)
  # Checks if horseshoes will be recalculated
  flag = true in [get_blade(self, i)._HSs==nothing for i in 1:self.B]

  # ERROR CASE IF NEEDS TO CALCULATE HORSESHOES
  if flag
    if self.RPM==nothing
      error("RPM hasn't been defined yet."*
            " Call function `setRPM()` before calling this function.")
    elseif self._wingsystem.Vinf==nothing
        error("Freestream hasn't been define yet, please call function set_Vinf()")
    end
  end

  # If horseshoes will be calculated, it forces to do it now and replaces
  # the regular infinite vortex with vortex in the direction of inflow at the
  # control point
  if flag
    for i in 1:self.B # Iterates over blades
      blade = get_blade(self, i)

      # Case horseshoes haven't been calculated yet
      if blade._HSs==nothing
        O = blade.O # Center of rotation

        # Forces to calculate horseshoes now
        _calculateHSs(blade; t=t, extraVinf...)

        # Calculates the inflow at each side Ap and Bp of each HS
        VAp = _calc_inflow(blade, get_RPM(self), t; target="Ap")
        VBp = _calc_inflow(blade, get_RPM(self), t; target="Bp")


        # Corrects each infinite vortex (infDA and infDB)
        for j in 1:size(blade._HSs)[1] # Iterates over horseshoes
          blade._HSs[j][6] = VAp[j]/norm(VAp[j])
          blade._HSs[j][7] = VBp[j]/norm(VBp[j])
        end
      end

    end
  end

  return getHorseshoe(self._wingsystem, m; t=t, extraVinf...)
end

"Saves the lofted surface of the blade"
function save_loft(self::Rotor, filename::String; addtiproot=false, path="", num=nothing, args...)
  # ERROR CASES
  if size(self.airfoils)[1]<2
    error("Requested lofted surface, but no airfoil geometry was given.")
  elseif size(self._polars)[1]<2
    error("Polars not found. Run `initialize()` and try again")
  end


  suf = "loft"

  CP_index = []       # Stores the CP index in order to hash the points
  lines = []          # Contour lines of cross sections

  # Iterates over each airfoil creating cross sections
  for (i,polar) in enumerate(self._polars)

    theta = pi/180*self._theta[i]   # Angle of attack

    # Actual airfoil contour
    x, y = self._chord[i]*polar.x, self._chord[i]*polar.y
    # Reformats x,y into point
    points = [ [x[i], y[i], 0.0] for i in 1:size(x)[1] ]
    # Rotates the contour in the right angle of attack
    # and orients the airfoil for CCW or CW rotation
    Oaxis = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1]
    Oaxis = [1 0 0; 0 (-1.0)^self.CW 0; 0 0 (-1.0)^self.CW]*Oaxis
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


    # Case of root or tip
    if addtiproot && (i==1 || i==size(self._polars,1))
      root_flg = i==1
      ind = root_flg ? 1 : size(self.r,1)
      alt_polar = root_flg ? self._polarroot : self._polartip

      theta = (-1)^(self.CW)*pi/180*self.theta[ind]   # Angle of attack

      # Actual airfoil contour
      x, y = self.chord[ind]*alt_polar.x, self.chord[ind]*alt_polar.y
      # Reformats x,y into point
      points = [ [x[i], y[i], 0.0] for i in 1:size(x)[1] ]
      # Rotates the contour in the right angle of attack
      # and orients the airfoil for CCW or CW rotation
      Oaxis = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1]
      Oaxis = [1 0 0; 0 (-1.0)^self.CW 0; 0 0 (-1.0)^self.CW]*Oaxis
      points = vtk.countertransform(points, inv(Oaxis), zeros(3))

      # Position of leading edge in FLOVLM blade's c.s.
      # Airfoil's x-axis = FLOWVLM blade's x-axis
      # Airfoil's y-axis = FLOWVLM blade's z-axis
      O = [(-1)*self.LE_x[ind], self.r[ind], (-1)^(self.CW)*self.LE_z[ind]]

      # Reformats the contour into FLOWVLM blade's c.s.
      points = [ O+[p[1], p[3], p[2]] for p in points]

      if root_flg
        lines = vcat([points], lines)
        CP_index = vcat([[i for p in points]], CP_index)
      else
        push!(lines, points)
        push!(CP_index, [i for p in points])
      end
    end
  end

  # Generates vtk cells from cross sections
  sections = [ [(1.0, 1, 1.0, false)] for i in 1:size(lines)[1]-1]
  points, vtk_cells, CP_index = vtk.multilines2vtkmulticells(lines, sections;
                                                      point_datas=CP_index)

  # Generates each blade
  for i in 1:self.B # Iterates over blades
    this_blade = self._wingsystem.wings[i]

    # Transforms points from FLOWVLM blade's c.s. to global c.s.
    this_points = FArrWrap[ vtk.countertransform(p, this_blade.invOaxis,
                                    this_blade.O) for p in points]

    # Formats the point data for generateVTK
    data = []
    # # Control point indices
    # push!(data, Dict(
    #             "field_name" => "ControlPoint_Index",
    #             "field_type" => "scalar",
    #             "field_data" => CP_index
    #             )
    #       )
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
function FLOWVLM2CCBlade(self::Rotor, RPM, blade_i::IWrap, turbine_flag::Bool;
                          sound_spd=nothing)
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

    r_over_R = self._r[i] / Rtip
    c_over_r = self._chord[i] / self._r[i]
    #   NOTE: Here I'm taking the freestream to be the absolute value CCBlade's
    #   x-component of inflow. This may cause problems when the flow is reversed
    this_Vinf = abs(inflows[i][1])
    tsr = this_Vinf < 1e-4 ? nothing : (2*pi*RPM/60 * Rtip) / this_Vinf

    # Mach correction
    if sound_spd!=nothing
      Ma = norm(inflows[i])/sound_spd
      alpha, cl = ap.get_cl(polar)
      geom = ap.get_geometry(polar)
      this_polar = ap.Polar(ap.get_Re(polar), alpha, cl/sqrt(1-Ma^2),
                              ap.get_cd(polar)[2], ap.get_cm(polar)[2],
                              geom[1], geom[2])
    else
      this_polar = polar
    end

    # 3D corrections
    this_polar = ap.correction3D(this_polar, r_over_R, c_over_r, tsr)

    # 360 extrapolation
    CDmax = 1.3
    this_polar = ap.extrapolate(this_polar, CDmax)

    # Makes sure the polar is injective for easing the spline
    this_polar = ap.injective(this_polar)

    # Converts to CCBlade's AirfoilData object
    alpha, cl = ap.get_cl(this_polar)
    _, cd = ap.get_cd(this_polar)
    ccb_polar = ccb.af_from_data(alpha, cl, cd; spl_k=5)

    push!(af, ccb_polar)
  end

  rotor = ccb.Rotor(self._r, self._chord,
                    (-1)^(turbine_flag)*(-1)^self.CW*self._theta*pi/180, af,
                    Rhub, Rtip, self.B, precone)

  return rotor
end


##### CALCULATION OF SOLUTION FIELDS ###########################################
"Receives the freestream velocity function V(x,t) and the current RPM of the
rotor, and it calculates the inflow velocity field that each control point
sees in the global coordinate system"
function calc_inflow(self::Rotor, Vinf, RPM; t::FWrap=0.0)
  omega = 2*pi*RPM/60

  data_Vtots = Array{FArrWrap}[]     # Inflow in the global c.s.
  data_Vccbs = Array{FArrWrap}[]     # Inflow in CCBlade's c.s.

  for (i,blade) in enumerate(self._wingsystem.wings) # Iterates over blades
    Vtots = FArrWrap[]
    Vccbs = FArrWrap[]

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

"Calculates the distributed loads from CCBlade. It also stores normal (Np) and
tangential (Tp) components relative to the plane of rotation as given by
CCBlade, if `include_comps==true`.

If `return_performance==true`, it returns propulsive efficiency `eta`,
thrust coefficient `CT`, and torque coefficient `CQ` of each blade.

NOTE: These loads are per unit length of span"
function calc_distributedloads(self::Rotor, Vinf, RPM, rho::FWrap;
                                t::FWrap=0.0, include_comps=false,
                                return_performance=false, Vref=nothing,
                                sound_spd=nothing)
  data = Array{FArrWrap}[]
  if include_comps
    data_Np = FArrWrap[]
    data_Tp = FArrWrap[]
  end

  turbine_flag = false  # This is a flag for ccblade to swap signs

  # Calculates inflows
  calc_inflow(self, Vinf, RPM; t=t)

  if return_performance; coeffs = []; end;

  # Distributed load on each blade
  for blade_i in 1:self.B

    # Formats the inflow as CCBlade's Inflow
    inflow = self.sol["CCBInflow"]["field_data"][blade_i]
    inflow_x = [V[1] for V in inflow]
    inflow_y = [V[2] for V in inflow]
    inflow_x = (-1)^(!turbine_flag)*inflow_x # The negative is needed to counteract the
    ccbinflow = ccb.Inflow(inflow_x, inflow_y, rho) # propeller swapping sign in CCBlade

    # Generates CCBlade Rotor object
    ccbrotor = FLOWVLM2CCBlade(self, RPM, blade_i, turbine_flag;
                                                            sound_spd=sound_spd)
    # Calls CCBlade
    # NOTE TO SELF: Forces normal and tangential to the plane of rotation
    Np, Tp, uvec, vvec = ccb.distributedloads(ccbrotor, ccbinflow, turbine_flag)

    # Convert forces from CCBlade's c.s. to global c.s.
    ccb_Fs = [ [Np[i], Tp[i], 0.0] for i in 1:get_mBlade(self) ]
    Fs = [ _ccblade2global( get_blade(self, blade_i), ccb_F, self.CW;
                          translate=false) for ccb_F in ccb_Fs ]
    # Stores the field
    push!(data, Fs)
    if include_comps
      push!(data_Np, Np)
      push!(data_Tp, Tp)
    end

    if return_performance
      if Vref==nothing
        error("Performance coefficients requested, but no reference freestream"*
              " has been provided.")
      end
      T, Q = ccb.thrusttorque(ccbrotor, [ccbinflow], turbine_flag)
      eta, CT, CQ = ccb.nondim(T, Q, Vref, 2*pi*RPM/60, rho,
                                ccbrotor.Rtip, ccbrotor.precone, turbine_flag)
      push!(coeffs, [eta,CT,CQ])
    end
  end

  # Adds solution fields
  field = Dict(
              "field_name" => "DistributedLoad",
              "field_type" => "vector",
              "field_data" => data
              )
  self.sol[field["field_name"]] = field
  if include_comps
    field = Dict(
                "field_name" => "Np",
                "field_type" => "scalar",
                "field_data" => data_Np
                )
    self.sol[field["field_name"]] = field
    field = Dict(
                "field_name" => "Tp",
                "field_type" => "scalar",
                "field_data" => data_Tp
                )
    self.sol[field["field_name"]] = field
  end

  if return_performance
    return coeffs
  else
    return nothing
  end
end

"Calculates sectional aerodynamic forces in a rotor where the field
DistributedLoad has already been solved for. It also calculates the bound
circulation Gamma"
function calc_aerodynamicforces(self::Rotor, rho::FWrap)
  if !("DistributedLoad" in keys(self.sol))
    error("Field `DistributedLoad` not found."*
          " Call `calc_distributedloads()` before calling this function.")
  end

  data_L = Array{FArrWrap}[]
  data_D = Array{FArrWrap}[]
  data_R = Array{FArrWrap}[]
  data_gamma = FArrWrap[]

  for blade_i in 1:self.B
    V = self.sol["GlobInflow"]["field_data"][blade_i]   # Inflow at each element
    Vmag = [norm(Vi) for Vi in V]
    F = self.sol["DistributedLoad"]["field_data"][blade_i]# Total force
    n_elem = size(F)[1]                                   # Number of elements

    # Unit directions
    dunit = V./Vmag                                       # Drag direction
    runit = [cross(F[i],V[i]) for i in 1:n_elem]          # Radial direction
    runit = runit./[norm(r) for r in runit]
    lunit = [cross(dunit[i], runit[i]) for i in 1:n_elem] # Lift direction
    lunit = lunit./[norm(l) for l in lunit]

    # Decomposition of total force
    Lmag = [dot(F[i], lunit[i]) for i in 1:n_elem]
    L = Lmag.*lunit                                       # Lift
    D = [dot(F[i], dunit[i]) for i in 1:n_elem].*dunit    # Drag
    R = [dot(F[i], runit[i]) for i in 1:n_elem].*runit    # Radial

    # Calculates circulation
    gamma = (-1)^(self.CW) * Lmag./(rho*Vmag)

    push!(data_L, L)
    push!(data_D, D)
    push!(data_R, R)
    push!(data_gamma, gamma)
  end

  # Adds solution fields
  field = Dict(
              "field_name" => "Lift",
              "field_type" => "vector",
              "field_data" => data_L
              )
  self.sol[field["field_name"]] = field
  field = Dict(
              "field_name" => "Drag",
              "field_type" => "vector",
              "field_data" => data_D
              )
  self.sol[field["field_name"]] = field
  field = Dict(
              "field_name" => "RadialForce",
              "field_type" => "vector",
              "field_data" => data_R
              )
  self.sol[field["field_name"]] = field
  field = Dict(
              "field_name" => "Gamma",
              "field_type" => "scalar",
              "field_data" => data_gamma
              )
  self.sol[field["field_name"]] = field

  return data_gamma
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

function _reset(self::Rotor;
                verbose=false, keep_Vinf=false, keep_RPM=true, keep_sol=true)
  _reset(self._wingsystem; verbose=verbose, keep_Vinf=keep_Vinf)
  _resetRotor(self; verbose=verbose, keep_RPM=keep_RPM)
end

function _resetRotor(self::Rotor; verbose=false, keep_RPM=false)
  if verbose; println("Resetting Rotor"); end;
  self.sol = Dict()
  if !keep_RPM; self.RPM = nothing; end;
end

"Generates the blade and discretizes it into lattices"
function _generate_blade(self::Rotor, n::IWrap; r::FWrap=1.0,
                          central=false, refinement=[], spl_k=5)

  # Splines
  spline_k = min(size(self.r)[1]-1, spl_k)
  spline_bc = "error"
  # spline_s = 0.001
  spline_s = 0.01
  _spl_chord = Dierckx.Spline1D(self.r, self.chord; k=spline_k,s=spline_s)
  _spl_theta = Dierckx.Spline1D(self.r, self.theta; k=spline_k,s=spline_s)
  _spl_LE_x = Dierckx.Spline1D(self.r, self.LE_x; k=spline_k,s=spline_s)
  _spl_LE_z = Dierckx.Spline1D(self.r, self.LE_z; k=spline_k,s=spline_s)
  spl_chord(x) = Dierckx.evaluate(_spl_chord, x)
  spl_theta(x) = (-1)^(self.CW)*Dierckx.evaluate(_spl_theta, x)
  spl_LE_x(x) =(-1)*Dierckx.evaluate(_spl_LE_x, x)
  spl_LE_z(x) = (-1)^(self.CW)*Dierckx.evaluate(_spl_LE_z, x)

  # Outputs
  out_r = FWrap[]         # Radial position of each control point
  out_chord = FWrap[]     # Chord at each control point
  out_theta = FWrap[]     # Angle of attach at each control point
  out_LE_x = FWrap[]
  out_LE_z = FWrap[]


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

"Verifies correct splining for lattice and element discretization"
function _verif_discr(self, blade, elem_r, elem_chord, elem_theta,
                                      elem_LE_x, elem_LE_z)

  # Original data
  r, chord, theta = self.r, self.chord, self.theta
  LE_x, LE_z = self.LE_x, self.LE_z
  Rtip = self.r[end]

  # Lattice discretization
  vlm_r, vlm_chord, vlm_theta, vlm_LE_x, vlm_LE_z = [[] for i in 1:5]
  for i in 1:get_m(blade)+1
    lex = blade._xlwingdcr[i]
    ley = blade._ywingdcr[i]
    lez = blade._zlwingdcr[i]
    tex = blade._xtwingdcr[i]
    tey = ley
    tez = blade._ztwingdcr[i]
    le = [lex,ley,lez]
    te = [tex,tey,tez]

    c = norm(le-te)
    tht = 180/pi*atan2((le-te)[3],-(le-te)[1])

    push!(vlm_r, ley)
    push!(vlm_chord, c)
    push!(vlm_theta, tht)
    push!(vlm_LE_x, lex)
    push!(vlm_LE_z, lez)
  end

  # Plots
  fig = figure("discret_verif", figsize=(7*2,5*2))
  for (i,(lbl, cr, cchord, ctheta, cLE_x, cLE_z)) in enumerate([
              ["Element", elem_r, elem_chord, elem_theta, elem_LE_x, elem_LE_z],
              ["Lattice", vlm_r, vlm_chord, vlm_theta, vlm_LE_x, vlm_LE_z]])
    subplot(220+2*(i-1)+1)

    title("Discretization Verification - $lbl")
    plot(r/Rtip, chord/Rtip, "ok", label="Chord data")
    plot(cr/Rtip, cchord/Rtip, "--or", label="Chord Spline")
    plot(r/Rtip, LE_x/Rtip, "^k", label="LE-x data")
    plot(cr/Rtip, -cLE_x/Rtip, "--^g", label="LE-x Spline")
    plot(r/Rtip, LE_z/Rtip, "*k", label="LE-z data")
    plot(cr/Rtip, cLE_z/Rtip, "--*b", label="LE-z Spline")
    xlabel(L"$r/R$")
    ylabel(L"$c/R$, $x/R$, $z/R$")
    legend(loc="best")
    grid(true, color="0.8", linestyle="--")

    subplot(220+2*(i-1)+2)
    plot(r/Rtip, theta, "ok", label="Twist data")
    plot(cr/Rtip, ctheta, "--^r", label="Twist Spline")
    xlabel(L"$r/R$")
    ylabel(L"Twist $\theta$ ($^\circ$)")
    legend(loc="best")
    grid(true, color="0.8", linestyle="--")
  end

end

"Calculates the airfoils at each control point"
function _calc_airfoils(self::Rotor, n::IWrap, r::FWrap,
                        central, refinement; rediscretize::Bool=true,
                        rfl_n_lower::IWrap=15, rfl_n_upper::IWrap=15,
                        rfl_r::FWrap=14.0, rfl_central::Bool=true)

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

function _rediscretize_airfoil(x, y, n_lower::IWrap, n_upper::IWrap, r::FWrap,
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

"""
Returns the inflow velocity at the requested target point on all horseshoes,
where the inflow is calculated as freestream + rotational velocity.
"""
function _calc_inflow(blade::Wing, RPM, t::FWrap; target="CP")
  if !(target in keys(VLMSolver.HS_hash))
    error("Logic error! Invalid target $target.")
  elseif blade._HSs==nothing
    error("Horseshoes haven't been initiated yet!")
  end
  t_i = VLMSolver.HS_hash[target]

  omega = 2*pi*RPM/60
  out = FArrWrap[]
  O = blade.O                   # Center of rotation
  runit = blade.Oaxis[2,:]      # Radial direction
  tunit = blade.Oaxis[1,:]      # Tangential direction

  for HS in blade._HSs # Iterates over horseshoes
    T = HS[t_i]                    # Targeted point

    # Calculates the rotational velocity at each side of the HS
    ## A side
    Tr = abs(dot(T-O, runit))      # Radius of rotation of target
    rotT = omega*Tr*tunit          # Rotational velocity

    # Total velocity
    VT = blade.Vinf(T, t) + rotT

    push!(out, VT)
  end

  return out
end

"""Receives a vector in the global coordinate system and transforms it into
CCBlade's coordinate system relative to `blade`. NOTE: This function only
rotates `V` into the new axis without translating it unless otherwise indicated.
(for definition of axes see notebook entry 20171202)
"""
function _global2ccblade(blade::Wing, V::FArrWrap, CW::Bool;
                                                        translate::Bool=false)
  # V in FLOWVLM Rotor's blade c.s.
  V_vlm = transform(V, blade.Oaxis, translate ? blade.O : zeros(3))

  # CCBlade c.s. transformation matrix
  ccb_Oaxis = _ccbladeOaxis(blade, CW)

  # V in CCBlade's c.s.
  V_ccb = transform(V_vlm, ccb_Oaxis, zeros(3))

  return V_ccb
end


"""Receives a vector in CCBlade's coordinate system relative to `blade` and
transforms it into the global coordinate system. NOTE: This function only
rotates `V` into the new axis without translating it unless otherwise indicated.
"""
function _ccblade2global(blade::Wing, V::FArrWrap, CW::Bool;
                                                        translate::Bool=false)
  # CCBlade c.s. transformation matrix
  ccb_Oaxis = _ccbladeOaxis(blade, CW)

  # V in FLOWVLM Rotor's blade c.s.
  V_vlm = countertransform(V, inv(ccb_Oaxis), zeros(3))

  # V in global c.s.
  V_glob = countertransform(V_vlm, blade.invOaxis, translate ? blade.O : zeros(3))

  return V_glob
end

"Returns the CCBlade's transformation matrix relative to the blade's c.s."
function _ccbladeOaxis(blade::Wing, CW::Bool)
  # CCBlade c.s. matrix
  ccb_Oaxis = [ 0 0 (-1)^CW;  # CC x-dir = Blade z-dir
                1.0 0 0;      # CC y-dir = Blade x-dir
                0 (-1)^CW 0]  # CC z-dir = Blade y-dir

  return ccb_Oaxis
end

"Extension of WingSystem's `addwing()` function"
function addwing(self::Rotor, wing_name::String, wing; force=false, args...)
  if !force
    error("This function is preserved for development."*
          " If doing so give `force=true`.")
  end
  addwing(self._wingsystem, wing_name, wing; args...)
end

"Extension of WingSystem's `get_wing()` function"
function get_wing(self::Rotor, args...; optargs...)
  return get_wing(self._wingsystem, args...; optargs...)
end

"Extension of WingSystem's `_addsolution()` function"
function _addsolution(self::Rotor, args...; optargs...)
  _addsolution(self._wingsystem, args...; optargs...)
end

"Extension of WingSystem's `_fetch_wing()` function"
function _fetch_wing(self::Rotor, args...)
  return _fetch_wing(self._wingsystem, args...)
end

"Extension of WingSystem's `_get_O()` function"
function _get_O(rotor::Rotor)
  return _get_O(rotor._wingsystem)
end

"Extension of WingSystem's `_get_O()` function"
function _get_Oaxis(rotor::Rotor)
  return _get_Oaxis(rotor._wingsystem)
end
##### END OF ROTOR CLASS #######################################################












#
