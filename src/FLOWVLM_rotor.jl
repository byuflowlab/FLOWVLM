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
mutable struct Rotor

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
  turbine_flag::Bool            # Whether this is a wind turbine or a propeller

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
          turbine_flag=false,
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
          turbine_flag,
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
                          central=false, refinement=[], verif=false,
                          figsize_factor=2/3, genblade_args=[], rfl_args...)
  # Checks for arguments consistency
  _check(self)

  # Flag for calculating airfoils
  rfl_flag = size(self.airfoils)[1]!=0

  # Generates blade
  blade, r, chord, theta, LE_x, LE_z = _generate_blade(self, n; r=r_lat,
                                        central=central, refinement=refinement,
                                        genblade_args...)
  self._r, self._chord, self._theta = r, chord, theta
  self._LE_x, self._LE_z = LE_x, LE_z
  self.m = get_m(blade)

  # Verifies correct lattice and blade elements discretization
  if verif; _verif_discr(self, blade, r, chord, theta, LE_x, LE_z; figsize_factor=figsize_factor); end;

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


"""Given the velocity induced at each control point (Vind = Vwake, no lifting
surface), solves for the Gamma field (circulation) on each blade by looking at
the airfoil polar at the effective angle of attack of every section. It also
includes the fields Ftot, L, D, and S. (WARNING: These Ftot, L, D, and S are
forces per unit length!)

This method solves iteratively until the circulation distribution converges.

NOTE: Vind is expected to be in the global coordinate system.
NOTE: Vind is expected to be formated as Vind[i][j] being the velocity vector
of the j-th control point in the i-th blade.
"""
function solvefromVite(self::Rotor, Vind::Array{Array{T, 1}, 1}, args...;
                          maxite::Int64=100, tol::Real=0.01, rlx=0.0, optargs...
                          ) where{T<:FArrWrap}

                          println(rlx)

  out = nothing
  if "Gamma" in keys(get_blade(self, 1).sol)
      prev_sol = [get_blade(self, j).sol["Gamma"] for j in 1:self.B]
  else
      prev_sol = nothing
  end
  ite = 0
  err = nothing

  for i in 1:maxite

    if i==1
      surfVind = Vind
    else
      # Adds V induced by lifting surface
      surfVind = [[ Vind[j][k] + Vind(self, getControlPoint(get_blade(self, j), k); ign_infvortex=true)
                            for k in 1:size(Vind[j],1)] for j in 1:size(Vind,1)]
    end

    out = solvefromV(self, surfVind, args...; optargs...)
    this_sol = [get_blade(self, j).sol["Gamma"] for j in 1:self.B]

    if prev_sol != nothing
      # Checking convergence: Average variation
      err = mean( [mean( abs.(prev_sol[j]-this_sol[j])./abs.(prev_sol[j]) ) for j in 1:self.B] )
      if err < tol
        break
      end

      # Relaxation
      for j in 1:rotor.B
        blade = vlm.get_blade(rotor, j)
        blade.sol["Gamma"][:] = rlx*prev_sol[j] .+ (1-rlx)*this_sol[j]
      end
    end

    prev_sol = this_sol
    ite += 1
  end

  if ite==maxite
    @warn "Iterative Rotor solvefromV reached max iterations without converging."*
            " maxite:$maxite\t error:$err"
  end

  return out
end

"""Given the velocity induced at each control point (Vind = Vliftsurface+Vwake),
solves for the Gamma field (circulation) on each blade by looking at the airfoil
polar at the effective angle of attack of every section. It also includes the
fields Ftot, L, D, and S. (WARNING: These Ftot, L, D, and S are
forces per unit length!)

THIS METHOD IS UNSTABLE.

NOTE: Vind is expected to be in the global coordinate system.
NOTE: Vind is expected to be formated as Vind[i][j] being the velocity vector
of the j-th control point in the i-th blade.
"""
function solvefromV(self::Rotor, Vind::Array{Array{T, 1}, 1}, args...;
                                                  optargs...) where{T<:FArrWrap}
  # ERROR CASES
  if size(Vind, 1)!=self.B
    error("Expected $(self.B) Vind entries; got $(size(Vind, 1)).")
  else
    for bi in 1:self.B
      if size(Vind[bi],1)!=get_mBlade(self)
        error("Expected $(get_mBlade(self)) Vind[$bi] entries;"*
              " got $(size(Vind[i],1)).")
      end
    end
  end

  return solvefromCCBlade(self, args...; _lookuptable=true, _Vinds=Vind,
                                                                    optargs...)
end

"Solves for the Gamma field (circulation) on each blade using CCBlade. It also
includes the fields Ftot, L, D, and S. (WARNING: These Ftot, L, D, and S are
forces per unit length!)

If include_comps==true it stores CCBlade-calculated normal and tangential forces
in the Rotor."
function solvefromCCBlade(self::Rotor, Vinf, RPM, rho::FWrap; t::FWrap=0.0,
                            include_comps::Bool=true, return_performance::Bool=false,
                            Vref=nothing, sound_spd=nothing, Uinds=nothing,
                            _lookuptable::Bool=false, _Vinds=nothing,
                            tiploss_correction=false, AR_to_360extrap=true)

  setVinf(self, Vinf)
  setRPM(self, RPM)

  # (Calls a HS to make sure they have been calculated)
  _ = getHorseshoe(self, 1)

  if sound_spd==nothing
    @warn "No sound speed has been provided. No Mach corrections will be applied."
  end

  # Calculates distributed load from CCBlade
  prfrmnc, gammas = calc_distributedloads(self, Vinf, RPM, rho; t=t,
                                        include_comps=include_comps,
                                        return_performance=return_performance,
                                        Vref=Vref, sound_spd=sound_spd,
                                        Uinds=Uinds,
                                        _lookuptable=_lookuptable, _Vinds=_Vinds,
                                        tiploss_correction=tiploss_correction,
                                        AR_to_360extrap=AR_to_360extrap)

  # Decomposes load into aerodynamic forces and calculates circulation
  gamma = calc_aerodynamicforces(self, rho; overwritegammas=gammas)

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
  # (WARNING: These Ftot, L, D, and S are forces per unit length!)
  _addsolution(self._wingsystem, "Ftot", new_Ftot; t=t)
  _addsolution(self._wingsystem, "L", new_L; t=t)
  _addsolution(self._wingsystem, "D", new_D; t=t)
  _addsolution(self._wingsystem, "S", new_S; t=t)

  return prfrmnc
end

"Sets Vinf(X,t) as the incoming freestream of this rotor"
function setVinf(self::Rotor, Vinf; keep_sol=false)
  _reset(self; keep_sol=keep_sol)
  setVinf(self._wingsystem, Vinf; keep_sol=keep_sol)
end

"Sets `RPM` as the revolutions per minutes of this rotor"
function setRPM(self::Rotor, RPM)
  _reset(self; keep_Vinf=true)
  _resetRotor(self)
  self.RPM = RPM
end

"Saves the rotor in VTK legacy format"
function save(self::Rotor, filename::String; addtiproot=true, airfoils=false,
                                     wopwop=false, wopbin=true, wopext="wop",
                                     wopv=1.0, save_horseshoes=true,
                                     args...)
  if save_horseshoes
      _ = getHorseshoe(self, 1) # Makes sure the wake is calculated right
  end

  strn = save(self._wingsystem, filename; save_horseshoes=save_horseshoes, args...)

  if size(self.airfoils)[1]!=0
    strn *= save_loft(self, filename; addtiproot=addtiproot, airfoils=airfoils,
                                wopwop=wopwop, wopbin=wopbin, wopext=wopext,
                                wopv=wopv, args...)
  end

  return strn
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
  rotOaxis = gt.rotation_matrix(0.0, 0.0, (-1)^!self.CW*degs)
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
function save_loft(self::Rotor, filename::String; addtiproot=false, path="",
                      num=nothing, airfoils=false,
                      wopwop=false, wopext="wop", wopbin=true, wopv=1.0,
                      args...)

  if wopwop; addtiproot=true; end;

  # ERROR CASES
  if size(self.airfoils)[1]<2
    error("Requested lofted surface, but no airfoil geometry was given.")
  elseif size(self._polars)[1]<2
    error("Polars not found. Run `initialize()` and try again")
  end

  strn = ""

  suf = "loft"
  rfl_suf = "rfl"

  CP_index = []       # Stores the CP index in order to hash the points
  lines = []          # Contour lines of cross sections

  # Iterates over each airfoil creating cross sections
  for (i,polar) in enumerate(self._polars)

    theta = pi/180*self._theta[i]   # Angle of attack

    # Actual airfoil contour
    x, y = self._chord[i]*polar.x, self._chord[i]*polar.y*(-1)^self.turbine_flag
    # Reformats x,y into point
    points = [ [x[i], y[i], 0.0] for i in 1:size(x)[1] ]
    # Rotates the contour in the right angle of attack
    # and orients the airfoil for CCW or CW rotation
    Oaxis = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1]
    Oaxis = [1 0 0; 0 (-1.0)^self.CW 0; 0 0 (-1.0)^self.CW]*Oaxis
    points = gt.countertransform(points, inv(Oaxis), fill(0.0, 3))

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
      points = gt.countertransform(points, inv(Oaxis), fill(0.0, 3))

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
  points, vtk_cells, CP_index = gt.multilines2vtkmulticells(lines, sections;
                                                      point_datas=CP_index)

  # Flips the cells in clockwise rotor to have normals pointing out
  if self.CW
      vtk_cells = reverse.(vtk_cells)
  end

  if airfoils || wopwop
    line_points, vtk_lines, _ = gt.lines2vtk(lines)
  end

  # Generates each blade
  for i in 1:self.B # Iterates over blades
    this_blade = self._wingsystem.wings[i]

    # Transforms points from FLOWVLM blade's c.s. to global c.s.
    this_points = FArrWrap[ gt.countertransform(p, this_blade.invOaxis,
                                    this_blade.O) for p in points]
    if airfoils || wopwop
      this_line_points = FArrWrap[ gt.countertransform(p, this_blade.invOaxis,
                                    this_blade.O) for p in line_points]
    end

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
      if field["field_type"] != "not-vtk"
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
    end

    # Generates the vtk file
    this_name = filename*"_"*self._wingsystem.wing_names[i]*"_"*suf
    strn *= gt.generateVTK(this_name, this_points; cells=vtk_cells,
                                            point_data=data, path=path, num=num)
    if airfoils
      this_linename = filename*"_"*self._wingsystem.wing_names[i]*"_"*rfl_suf
      strn *= gt.generateVTK(this_linename, this_line_points; cells=vtk_lines,
                                path=path, num=num)
    end

    # Generates PLOT3D-like geometry for PSU-WOPWOP
    if wopwop

        # PRECALCULATIONS
        nc = size(vtk_cells, 1)             # Number of cells
        CPs = fill(0.0, 3, nc)                  # Control point of every cell
        Ns = fill(0.0, 3, nc)                   # Normal of every cell
        for k in 1:nc

            p = [this_points[j+1] for j in vtk_cells[k]]    # Points of cell

            # NOTE: Here we assume that cells are quadrilaterals
            # NOTE: Here the negative sign is necessary because it turns out
            # `vtk_cells` are going clockwise
            crss1 = -cross(p[2]-p[1], p[3]-p[1])
            crss2 = -cross(p[4]-p[3], p[1]-p[3])

            # # Area
            # A1 = 0.5*norm(crss1)                            # Area of triangle 1
            # A2 = 0.5*norm(crss2)                            # Area of triangle 2
            # A = A1+A2                                       # Area quadrilateral

            # # Normal
            # N1 = crss1/(2*A1)
            # N2 = crss2/(2*A2)
            # N = (A1*N1 + A2*N2)/A

            # Normal scaled by area
            # NA = N*A
            NA = crss1/2 + crss2/2
            Ns[:, k] = NA

            # Centroid (control point)
            for pnt in p
                CPs[:, k] .+= pnt
            end
            CPs[:, k] /= size(p, 1)
        end

        nHS = get_m(this_blade)                     # Number of horseshoes
        Cs = fill(0.0, 3, nHS)                          # Midpoints of lifting line
        NCs = fill(0.0, 3, nHS)                         # Ficticious normals to line
        lift_points = []                            # Lifting line points
        lift_vtk_cells = []                         # VTK lifting line
        for k in 1:nHS
            Ap, A, B, Bp, _, _, _, _ = getHorseshoe(this_blade, k)
            Cs[:, k] = (A+B)/2

            # Estimate the area of every triangle of the horseshoe extending the
            # bound vortex to the leading edge
            crss1 = cross((Ap-A)/(1-pn), B-A)
            crss2 = cross(A-B, (Bp-B)/(1-pn))
            NCs[:, k] = crss1/2 + crss2/2

            # Unit vector normals
            NCs[:, k] /= norm(NCs[:, k])

            # Normals scaled by length
            NCs[:, k] *= norm(B-A)

            push!(lift_vtk_cells, (size(lift_points, 1), size(lift_points, 1)+1))
            push!(lift_points, A)
            # NOTE: Here I'm assuming the horseshoes are contiguous
            if k==nHS
                push!(lift_points, B)
            end
        end

        # ----- Root and tip end caps ---------
        points_root = [this_line_points[pi+1] for pi in vtk_lines[1]]
        points_tip = [this_line_points[pi+1] for pi in vtk_lines[end]]

        # Forces end caps to be even number of points
        if size(points_root, 1)%2 != 0
            points_root = points_root[1:end-1]
        end
        if size(points_tip, 1)%2 != 0
            points_tip = points_tip[1:end-1]
        end
        npr = size(points_root, 1)              # Number of points in root
        npt = size(points_tip, 1)               # Number of points in tip
        ncr = Int(size(points_root, 1)/2)-1     # Number of cells in root
        nct = Int(size(points_tip, 1)/2)-1      # Number of cells in tip

        # 1-indexed cells
        cells_root = [[ci, ci+1, (npr+1)-(ci+1), (npr+1)-ci] for ci in 1:ncr]
        cells_tip =  [[ci, ci+1, (npt+1)-(ci+1), (npt+1)-ci] for ci in 1:nct]

        # Normals scaled by area
        # NOTE: Direction of normal vs clockwise nodes?
        Nroot = fill(0.0, 3, ncr)
        for ci in 1:ncr
            p1, p2, p3, p4 = (points_root[pi] for pi in cells_root[ci])
            crss1 = cross(p2-p1, p3-p1)
            crss2 = cross(p4-p3, p1-p3)
            Nroot[:, ci] = crss1/2 + crss2/2
        end
        Ntip = fill(0.0, 3, nct)
        for ci in 1:nct
            p1, p2, p3, p4 = (points_tip[pi] for pi in cells_tip[ci])
            crss1 = -cross(p2-p1, p3-p1)
            crss2 = -cross(p4-p3, p1-p3)
            Ntip[:, ci] = crss1/2 + crss2/2
        end


        # ------------------------------------------------------------------
        # ----------------- LOFTED BLADE FOR THICKNESS ---------------------
        # ------------------------------------------------------------------
        # Create file
        f = open(joinpath(path,
                    this_name*"."*(num!=nothing ? "$(num)." : "")*wopext), "w")

        # Binary / ASCII printing
        prnt(x) = wopbin ? write(f, x) : print(f, x)
        prntln(x) = wopbin ? write(f, x) : print(f, x, "\n")

        # Convertion to 4-bytes numbers
        # NOTE: 4 bytes = 4*8 bites = 32 bites
        fl(x) = Float32(x)
        nt(x) = Int32(x)
        # Convertion to n-bytes string
        st(x::String, n) = x * " "^(n-length(x))

        if  wopv==0.0
            # Patch declaration line
            prntln("Patch"*" "^27)
            # Number of zones
            prntln(nt(1))

            # ----------------- FIRST PATCH: LOFT ------------------------------
            # imax jmax
            imax = self.m-1 + 2*addtiproot
            jmax = Int(nc/imax)
            prnt(nt(imax))
            prntln(nt(jmax))

            # imax × jmax floating point x coordinates
            # imax × jmax floating point y coordinates
            # imax × jmax floating point z coordinates
            for k in 1:3
                for j in 1:nc
                    prntln(fl(CPs[k, j]))
                end
            end

            # imax × jmax floating point normal vector x coordinates
            # imax × jmax floating point normal vector y coordinates
            # imax × jmax floating point normal vector z coordinates
            for k in 1:3
                for j in 1:nc
                    prntln(fl(Ns[k, j]))
                end
            end

        elseif wopv==1.0
            # Magic number
            prntln(nt(42))
            # Version number
            prnt(nt(1))
            prntln(nt(0))
            # Units for Tecplot
            prntln(st("N/m^2", 32))
            # Comments
            prntln(st("Geometry input file for PSU-WOPWOP (Format v1.0)\n"*
                      "------------------------------------------------\n"*
                      "Created by FLOWVLM (written by Eduardo Alvarez)\n"*
                      "https://github.com/byuflowlab/FLOWVLM\n"*
                      "Creation date: $(Dates.now())\n"*
                      "Units: SI\n"*
                      "Format: Unstructured grid, face-centered", 1024))

            # Format string
            prntln(nt(1))               # Geometry file flag
            prntln(nt(3))               # Number of zones
            prntln(nt(2))               # 1==structured, 2==unstructured
            prntln(nt(1))               # Geometry 1==constant, 2==periodic, 3==aperiodic
            prntln(nt(2))               # Normal vectors 1==node, 2==face
            prntln(nt(1))               # Floating point 1==single, 2==double
            prntln(nt(0))               # iblank values 1==included, 0==not
            prntln(nt(0))               # WOPWOP secret conspiracy

            # ----------------- FIRST PATCH: LOFT ------------------------------
            # Name
            prntln(st("loft", 32))
            # nbNodes
            prntln(nt( size(this_points, 1) ))
            # nbFaces
            prntln(nt( size(vtk_cells, 1) ))
            # Connectivity
            for cell in vtk_cells
                prnt(nt( size(cell, 1) ))         # Number of nodes in this cell
                # for pi in reverse(cell)         # Clockwise node ordering
                for pi in cell                    # NOTE: Turns out that `vtk_cells` are already clockwise
                    prnt(nt( pi+1 ))              # 1-indexed node index
                end
                if !wopbin; prntln(""); end;
            end

            # ----------------- SECOND PATCH: ROOT CAP -------------------------
            # Name
            prntln(st("root", 32))
            # nbNodes
            prntln(nt( size(points_root, 1) ))
            # nbFaces
            prntln(nt( size(cells_root, 1) ))
            # Connectivity
            for cell in cells_root
                prnt(nt( size(cell, 1) ))         # Number of nodes in this cell
                for pi in reverse(cell)           # NOTE: Clockwise node ordering
                # for pi in cell
                    prnt(nt( pi ))                # 1-indexed node index
                end
                if !wopbin; prntln(""); end;
            end

            # ----------------- THIRD PATCH: TIP CAP ---------------------------
            # Name
            prntln(st("tip", 32))
            # nbNodes
            prntln(nt( size(points_tip, 1) ))
            # nbFaces
            prntln(nt( size(cells_tip, 1) ))
            # Connectivity
            for cell in cells_tip
                prnt(nt( size(cell, 1) ))         # Number of nodes in this cell
                # for pi in reverse(cell)           # NOTE: Clockwise node ordering
                for pi in cell
                    prnt(nt( pi ))                # 1-indexed node index
                end
                if !wopbin; prntln(""); end;
            end

            # ----------------- DATA FIRST PATCH -------------------------------
            # nbNodes floating point x coordinates
            # nbNodes floating point y coordinates
            # nbNodes floating point z coordinates
            for k in 1:3
                for p in this_points
                    prntln(fl(p[k]))
                end
            end

            # nbFaces floating point normal vector x coordinates
            # nbFaces floating point normal vector y coordinates
            # nbFaces floating point normal vector z coordinates
            for k in 1:3
                for j in 1:size(vtk_cells, 1)
                    prntln(fl(Ns[k, j]))
                end
            end

            # ----------------- DATA SECOND PATCH -------------------------------
            # nbNodes floating point x coordinates
            # nbNodes floating point y coordinates
            # nbNodes floating point z coordinates
            for k in 1:3
                for p in points_root
                    prntln(fl(p[k]))
                end
            end

            # nbFaces floating point normal vector x coordinates
            # nbFaces floating point normal vector y coordinates
            # nbFaces floating point normal vector z coordinates
            for k in 1:3
                for j in 1:size(cells_root, 1)
                    prntln(fl(Nroot[k, j]))
                end
            end

            # ----------------- DATA THIRD PATCH -------------------------------
            # nbNodes floating point x coordinates
            # nbNodes floating point y coordinates
            # nbNodes floating point z coordinates
            for k in 1:3
                for p in points_tip
                    prntln(fl(p[k]))
                end
            end

            # nbFaces floating point normal vector x coordinates
            # nbFaces floating point normal vector y coordinates
            # nbFaces floating point normal vector z coordinates
            for k in 1:3
                for j in 1:size(cells_tip, 1)
                    prntln(fl(Ntip[k, j]))
                end
            end

        else
            error("Got invalid WOPWOP version $wopv")
        end

        close(f)



        # ------------------------------------------------------------------
        # ----------------- LIFTING-LINE COMPACT PATCH FOR LOADING ---------
        # ------------------------------------------------------------------
        # Create file
        f = open(joinpath(path,
                    this_name*"_compact"*"."*(num!=nothing ? "$(num)." : "")*wopext), "w")

        if  wopv==0.0
            # Patch declaration line
            prntln("Patch"*" "^27)
            # Number of zones
            prntln(nt(1))

            # ----------------- FIRST PATCH: LIFTING LINE ---------------------
            # imax jmax
            prnt(nt(nHS))
            prntln(nt(1))

            # imax × jmax floating point x coordinates
            # imax × jmax floating point y coordinates
            # imax × jmax floating point z coordinates
            for k in 1:3
                for j in 1:nHS
                    prntln(fl(Cs[k, j]))
                end
            end

            # imax × jmax floating point normal vector x coordinates
            # imax × jmax floating point normal vector y coordinates
            # imax × jmax floating point normal vector z coordinates
            for k in 1:3
                for j in 1:nHS
                    prntln(fl(NCs[k, j]))
                end
            end

        elseif wopv==1.0
            # Magic number
            prntln(nt(42))
            # Version number
            prnt(nt(1))
            prntln(nt(0))
            # Units for Tecplot
            prntln(st("N/m^2", 32))
            # Comments
            prntln(st("Geometry input file for PSU-WOPWOP (Format v1.0)\n"*
                      "------------------------------------------------\n"*
                      "Created by FLOWVLM (written by Eduardo Alvarez)\n"*
                      "https://github.com/byuflowlab/FLOWVLM\n"*
                      "Creation date: $(Dates.now())\n"*
                      "Units: SI\n"*
                      "Format: Structured grid, node-centered", 1024))

            # Format string
            prntln(nt(1))               # Geometry file flag
            prntln(nt(1))               # Number of zones
            prntln(nt(1))               # 1==structured, 2==unstructured
            prntln(nt(1))               # Geometry 1==constant, 2==periodic, 3==aperiodic
            prntln(nt(1))               # Normal vectors 1==node, 2==face
            prntln(nt(1))               # Floating point 1==single, 2==double
            prntln(nt(0))               # iblank values 1==included, 0==not
            prntln(nt(0))               # WOPWOP secret conspiracy

            # ----------------- FIRST PATCH: LIFTING LINE ----------------------
            # Name
            prntln(st("liftingline", 32))
            # iMax
            prntln(nt( 1 ))
            # jMax
            prntln(nt( size(lift_points, 1) ))

            # ----------------- DATA FIRST PATCH -------------------------------
            # imax × jmax floating point x coordinates
            # imax × jmax floating point y coordinates
            # imax × jmax floating point z coordinates
            for k in 1:3
                for p in lift_points
                    prntln(fl(p[k]))
                end
            end

            # imax × jmax floating point normal vector x coordinates
            # imax × jmax floating point normal vector y coordinates
            # imax × jmax floating point normal vector z coordinates
            for k in 1:3
                for j in 1:nHS
                    # NOTE: Here I assume contiguous horseshoes

                    if j == 1
                        prntln(fl(NCs[k, j]/2))
                    end

                    if j != nHS
                        prntln(fl(NCs[k, j]/2 + NCs[k, j+1]/2))
                    else
                        prntln(fl(NCs[k, j]/2))
                    end

                end
            end

        else
            error("Got invalid WOPWOP version $wopv")
        end

        close(f)
    end

  end

  return strn
end




##### CALCULATION OF SOLUTION FIELDS ###########################################
"Receives the freestream velocity function V(x,t) and the current RPM of the
rotor, and it calculates the inflow velocity field that each control point
sees in the global coordinate system"
function calc_inflow(self::Rotor, Vinf, RPM; t::FWrap=0.0, Vinds=nothing)
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
      this_Vrot = countertransform(this_Vrot, blade.invOaxis, fill(0.0, 3))

      this_Vtot = this_Vinf + this_Vrot

      # Adds any extra induced velocity
      if Vinds!=nothing
        this_Vtot += Vinds[i][j]
      end

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
                                t::FWrap=0.0, include_comps=true,
                                return_performance=false, Vref=nothing,
                                Uinds=nothing,
                                sound_spd=nothing,
                                _lookuptable::Bool=false, _Vinds=nothing,
                                tiploss_correction::Bool=false,
                                AR_to_360extrap = true)
  data = Array{FArrWrap}[]
  if include_comps
    data_Np     = FArrWrap[]
    data_Tp     = FArrWrap[]
    data_u      = FArrWrap[]
    data_v      = FArrWrap[]
    data_cl     = FArrWrap[]
    data_cd     = FArrWrap[]
    data_cn     = FArrWrap[]
    data_ct     = FArrWrap[]
    if !_lookuptable
        data_a      = FArrWrap[]
        data_ap     = FArrWrap[]
        data_phi    = FArrWrap[]
        data_alpha  = FArrWrap[]
        data_W      = FArrWrap[]
        data_F      = FArrWrap[]
        data_G      = FArrWrap[]
    end
  end
  ccbrotors, ccbsectionss, ccbopss = [], [], []

  gammas = _lookuptable ? [] : nothing

  turbine_flag = self.turbine_flag

#   turbine_flag = true  # This is a flag for ccblade to swap signs

  # Calculates inflows
  # calc_inflow(self, Vinf, RPM; t=t, Vinds=(_lookuptable ? _Vinds : nothing) )
  calc_inflow(self, Vinf, RPM; t=t, Vinds=_Vinds)

  if return_performance; coeffs = []; end;

  # Distributed load on each blade
  for blade_i in 1:self.B

    # Formats the inflow as CCBlade's Inflow
    inflow = self.sol["CCBInflow"]["field_data"][blade_i]
    inflow_x = [V[1] for V in inflow]
    inflow_y = [V[2] for V in inflow]
    inflow_x = (-1)^(!turbine_flag)*inflow_x # The negative is needed to counteract the
    occbinflow = OCCBInflow(inflow_x, inflow_y, rho) # propeller swapping sign in CCBlade

    # Generates old-CCBlade Rotor object
    occbrotor = FLOWVLM2OCCBlade(self, RPM, blade_i, turbine_flag;
                                                            sound_spd=sound_spd, AR_to_360extrap=AR_to_360extrap)
    # Convert old-CCBlade rotor to current CCBlade rotor type
    ccbrotor, ccbsections, ccbops = OCCB2CCB(occbrotor, turbine_flag,
                                                        occbinflow; pitch=0.0)
    ccboutputs = nothing
    Np, Tp, uvec, vvec, gamma = (nothing for i in 1:5)

    if _lookuptable
      Np, Tp, uvec, vvec, gamma, cn, ct, cl, cd = _calc_distributedloads_lookuptable(occbrotor,
                                                        occbinflow, turbine_flag;
                                                        tiploss_correction=tiploss_correction)
      push!(gammas, gamma)
      zeroarr = zeros(size(Np))
      ccbouputs = ccb.Outputs(Np, Tp, zeroarr, zeroarr, uvec, vvec,
                                zeroarr, zeroarr, zeroarr, cl,
                                cd, cn, ct, zeroarr, zeroarr)
    else
      # Calls CCBlade
      # NOTE TO SELF: Forces normal and tangential to the plane of rotation
      # Np, Tp, uvec, vvec = ccb.distributedloads(ccbrotor, ccbinflow, turbine_flag)
      ccboutputs = ccb.solve.(Ref(ccbrotor), ccbsections, ccbops)
      Np, Tp, a, ap, uvec, vvec, phi, alpha, W, cl, cd, cn, ct, loss, effloss = ccboutputs.Np, ccboutputs.Tp, ccboutputs.a, ccboutputs.ap, ccboutputs.u, ccboutputs.v, ccboutputs.phi, ccboutputs.alpha, ccboutputs.W, ccboutputs.cl, ccboutputs.cd, ccboutputs.cn, ccboutputs.ct, ccboutputs.F, ccboutputs.G
    end

    # Convert forces from CCBlade's c.s. to global c.s.
    ccb_Fs = [ [Np[i], Tp[i], 0.0] for i in 1:get_mBlade(self) ]
    Fs = [ _ccblade2global( get_blade(self, blade_i), ccb_F, self.CW;
                          translate=false) for ccb_F in ccb_Fs ]
    # Stores the field
    push!(data, Fs)
    if include_comps
      push!(data_Np     , Np)
      push!(data_Tp     , Tp)
      push!(data_u      , uvec)
      push!(data_v      , vvec)
      push!(data_cl     , cl)
      push!(data_cd     , cd)
      push!(data_cn     , cn)
      push!(data_ct     , ct)
        if !_lookuptable
            push!(data_a      , a)
            push!(data_ap     , ap)
            push!(data_phi    , phi)
            push!(data_alpha  , alpha)
            push!(data_W      , W)
            push!(data_F      , loss)
            push!(data_G      , effloss)
        end
      push!(ccbrotors, ccbrotor)
      push!(ccbsectionss, ccbsections)
      push!(ccbopss, ccbops)
    end

    if return_performance
      if Vref==nothing
        error("Performance coefficients requested, but no reference freestream"*
              " has been provided.")
      end
      # T, Q = ccb.thrusttorque(ccbrotor, [ccbinflow], turbine_flag)
      # eta, CT, CQ = ccb.nondim(T, Q, Vref, 2*pi*RPM/60, rho,
      #                           ccbrotor.Rtip, ccbrotor.precone, turbine_flag)
      #NOTE: Define ccboutputs
      T, Q = ccb.thrusttorque(ccbrotor, ccbsections, ccboutputs)
      eta, CT, CQ = ccb.nondim(T, Q, Vref, 2*pi*RPM/60, rho, ccbrotor)
      push!(coeffs, [eta,CT,CQ])
    end

    if Uinds!=nothing; push!(Uinds, [uvec, vvec]); end;
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
    field = Dict(
        "field_name" => "cl",
        "field_type" => "scalar",
        "field_data" => data_cl
        )
    self.sol[field["field_name"]] = field
    field = Dict(
            "field_name" => "cd",
            "field_type" => "scalar",
            "field_data" => data_cd
            )
    self.sol[field["field_name"]] = field
    field = Dict(
            "field_name" => "cn",
            "field_type" => "scalar",
            "field_data" => data_cn
            )
    self.sol[field["field_name"]] = field
    field = Dict(
            "field_name" => "ct",
            "field_type" => "scalar",
            "field_data" => data_ct
            )
    self.sol[field["field_name"]] = field
    if !_lookuptable
        field = Dict(
                "field_name" => "u",
                "field_type" => "scalar",
                "field_data" => data_u
                )
        self.sol[field["field_name"]] = field
        field = Dict(
                    "field_name" => "v",
                    "field_type" => "scalar",
                    "field_data" => data_v
                    )
        self.sol[field["field_name"]] = field
        field = Dict(
                    "field_name" => "a",
                    "field_type" => "scalar",
                    "field_data" => data_a
                    )
        self.sol[field["field_name"]] = field
        field = Dict(
                    "field_name" => "ap",
                    "field_type" => "scalar",
                    "field_data" => data_ap
                    )
        self.sol[field["field_name"]] = field
        field = Dict(
                    "field_name" => "phi",
                    "field_type" => "scalar",
                    "field_data" => data_phi
                    )
        self.sol[field["field_name"]] = field
        field = Dict(
                    "field_name" => "alpha",
                    "field_type" => "scalar",
                    "field_data" => data_alpha
                    )
        self.sol[field["field_name"]] = field
        field = Dict(
                    "field_name" => "W",
                    "field_type" => "scalar",
                    "field_data" => data_W
                    )
        self.sol[field["field_name"]] = field
        field = Dict(
                    "field_name" => "F",
                    "field_type" => "scalar",
                    "field_data" => data_F
                    )
        self.sol[field["field_name"]] = field
        field = Dict(
                    "field_name" => "G",
                    "field_type" => "scalar",
                    "field_data" => data_G
                    )
        self.sol[field["field_name"]] = field
    end
    field = Dict(
                "field_name" => "ccbrotor",
                "field_type" => "not-vtk",
                "field_data" => ccbrotors
                )
    self.sol[field["field_name"]] = field
    field = Dict(
                "field_name" => "ccbsections",
                "field_type" => "not-vtk",
                "field_data" => ccbsectionss
                )
    self.sol[field["field_name"]] = field
    field = Dict(
                "field_name" => "ccbops",
                "field_type" => "not-vtk",
                "field_data" => ccbopss
                )
    self.sol[field["field_name"]] = field
  end

  if return_performance
    return coeffs, gammas
  else
    return nothing, gammas
  end
end

"Calculates sectional aerodynamic forces in a rotor where the field
DistributedLoad has already been solved for. It also calculates the bound
circulation Gamma"
function calc_aerodynamicforces(self::Rotor, rho::FWrap; overwritegammas=nothing)
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
    if overwritegammas!=nothing
      gamma = (-1)^(self.CW) * overwritegammas[blade_i]
    else
      bladeaxis = get_blade(self, blade_i).Oaxis[2, :]    # Blade axis
      Lsgn = [sign(dot(ru, bladeaxis)) for ru in runit]   # Sign of lift
      gamma = (-1)^(self.CW) * Lsgn .* Lmag./(rho*Vmag)
    end

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

"""
  `calc_thrust_torque(rotor)`

Integrates the load distribution along every blade to return the thrust and
torque of the rotor.
"""
function calc_thrust_torque(self::Rotor;)
  # Error cases
  if !("Np" in keys(self.sol))
    error("Thrust calculation requested, but Np field not found."*
          " Call `solvefromV()` or `solvefromCCBlade` before calling this"*
          " function")
  elseif !("Tp" in keys(self.sol))
    error("Torque calculation requested, but Tp field not found."*
          " Call `solvefromV()` or `solvefromCCBlade` before calling this"*
          " function")
  end

  thrust = 0.0
  torque = 0.0

  # Radial length of every horseshoe
  lengths = Float64[2*(self._r[1]-self.hubR)]
  for i in 2:get_mBlade(self)
    push!(lengths, 2*(self._r[i]-self._r[i-1])-lengths[end] )
  end

  # Verifying the logic
  if abs((self.hubR + sum(lengths))/self.rotorR - 1) > 1e-7
    error("LOGIC ERROR: Sum of lengths don't add up the radius."*
          " Resulted in $((self.hubR + sum(lengths))), expected $self.rotorR.")
  end

  # Iterates over every blade
  for blade_i in 1:self.B
    Np = self.sol["Np"]["field_data"][blade_i]
    Tp = self.sol["Tp"]["field_data"][blade_i]

    # Iterates over every horseshoe
    for j in 1:get_mBlade(self)
      # Integrates over this horseshoe
      thrust += Np[j]*lengths[j]
      torque += Tp[j]*lengths[j]*self._r[j]
    end
  end

  return thrust, torque
end


"""
  `calc_thrust_torque_coeffs(rotor, rho, t)`

Integrates the load distribution along every blade to return the thrust and
torque coefficients of the rotor. Thrust coefficient CT calculated as
`CT=\frac{T}{\rho n^2 D^4}`, torque coefficient CQ calculated as
`CQ=\frac{Q}{\rho n^2 D^5}`
"""
function calc_thrust_torque_coeffs(self::Rotor, rho::Real)
    thrust, torque = calc_thrust_torque(self)
    n = self.RPM / 60
    D = 2*self.rotorR

    return thrust/(rho*n^2*D^4), torque/(rho*n^2*D^5)
end

function calc_thrust_torque_coeffs(self::Rotor, rho::Real, Vinf::Real, turbine_flag::Bool)
  thrust, torque = calc_thrust_torque(self)
  n = self.RPM / 60
  D = 2*self.rotorR
  if turbine_flag == false
    return thrust/(rho*n^2*D^4), torque/(rho*n^2*D^5)
  else
    q = 0.5*rho*Vinf^2
    A = pi*self.rotorR^2
    return thrust/(q*A), torque/(q*self.rotorR*A)
  end
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
  _reset(self._wingsystem; verbose=verbose, keep_Vinf=keep_Vinf,
                                                            keep_sol=keep_sol)
  _resetRotor(self; verbose=verbose, keep_RPM=keep_RPM)
end

function _resetRotor(self::Rotor; verbose=false, keep_RPM=false)
  if verbose; println("Resetting Rotor"); end;
  self.sol = Dict()
  if !keep_RPM; self.RPM = nothing; end;
end

"Generates the blade and discretizes it into lattices"
function _generate_blade(self::Rotor, n::IWrap; r::FWrap=1.0,
                          central=false, refinement=[], spl_k=5, spl_s=0.01)

  # Splines
  spline_k = min(size(self.r)[1]-1, spl_k)
  spline_bc = "error"
  # spline_s = 0.001
  spline_s = spl_s
  _spl_chord = Dierckx.Spline1D(self.r, self.chord; k=spline_k,s=spline_s)
  _spl_theta = Dierckx.Spline1D(self.r, self.theta; k=spline_k,s=spline_s)
  _spl_LE_x = Dierckx.Spline1D(self.r, self.LE_x; k=spline_k,s=spline_s)
  _spl_LE_z = Dierckx.Spline1D(self.r, self.LE_z; k=spline_k,s=spline_s)
  spl_chord(x) = Dierckx.evaluate(_spl_chord, x)
  spl_theta(x) = (-1)^(self.CW)*Dierckx.evaluate(_spl_theta, x)
  spl_theta2(x) = (-1)^(self.turbine_flag)*spl_theta(x)
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
                spl_chord(self.r[1]), spl_theta2(self.r[1]))

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
              spl_chord(this_r), spl_theta2(this_r), 1)

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
                                      elem_LE_x, elem_LE_z; figsize_factor=2/3)

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
    tht = 180/pi*atan((le-te)[3],-(le-te)[1])

    push!(vlm_r, ley)
    push!(vlm_chord, c)
    push!(vlm_theta, tht)
    push!(vlm_LE_x, lex)
    push!(vlm_LE_z, lez)
  end

  # Plots
  fig = figure("discret_verif", figsize=[7*2,5*2]*figsize_factor)
  axs = fig.subplots(2, 2)

  for (i,(lbl, cr, cchord, ctheta, cLE_x, cLE_z)) in enumerate([
              ["Element", elem_r, elem_chord, elem_theta, elem_LE_x, elem_LE_z],
              ["Lattice", vlm_r, vlm_chord, vlm_theta, vlm_LE_x, vlm_LE_z]])

    ax = axs[2*(i-1)+1]
    ax.title.set_text("Discretization Verification - $lbl")
    ax.plot(r/Rtip, chord/Rtip, "ok", label="Chord data", alpha=0.75)
    ax.plot(cr/Rtip, cchord/Rtip, "--or", label="Chord Spline", alpha=0.75)
    ax.plot(r/Rtip, LE_x/Rtip, "^k", label="LE-x data", alpha=0.75)
    ax.plot(cr/Rtip, -cLE_x/Rtip, "--^g", label="LE-x Spline", alpha=0.75)
    ax.plot(r/Rtip, LE_z/Rtip, "*k", label="LE-z data", alpha=0.75)
    ax.plot(cr/Rtip, cLE_z/Rtip, "--*b", label="LE-z Spline", alpha=0.75)
    ax.set_xlabel(L"$r/R$")
    ax.set_ylabel(L"$c/R$, $x/R$, $z/R$")
    ax.legend(loc="best", frameon=false)
    ax.grid(true, color="0.8", linestyle="--")

    ax = axs[2*(i-1)+2]
    ax.plot(r/Rtip, theta, "ok", label="Twist data", alpha=0.75)
    ax.plot(cr/Rtip, ctheta, "--^r", label="Twist Spline", alpha=0.75)
    ax.set_xlabel(L"$r/R$")
    ax.set_ylabel(L"Twist $\theta$ ($^\circ$)")
    ax.legend(loc="best", frameon=false)
    ax.grid(true, color="0.8", linestyle="--")
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
    points = gt.multidiscretize(aux_f, 0, 1,
          [(sec[1], ns[i], sec[3], false) for (i,sec) in enumerate(refinement)])
  else
    points = gt.discretize(aux_f, 0, 1, n, r; central=central)
  end

  # Normalized position of control points
  norm_CPs = [ (points[i]+points[i-1])/2 for i in 2:size(points)[1] ]

  # ------ Iterates over CPs calculating the polar of the airfoil at that CP ---

  prev_contour, next_contour = nothing, nothing # Airfoil contours
  prev_i, next_i = nothing, 1 # Indices of airfoils previous and next to cur CP

  # WARNING: The blending function only blends airfoil properties and contours,
  #       but it ignores all other parameters (Re, Ma, etc)

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
    blended_pypolar = prev_pypolar.blend(next_pypolar, weight)

    # Blends airfoil geometry
    blended_x = weight*next_contour[1] + (1-weight)*prev_contour[1]
    blended_y = weight*next_contour[2] + (1-weight)*prev_contour[2]

    # Blended Polar object
    blended_polar = ap._pyPolar2Polar(blended_pypolar; x=blended_x, y=blended_y)

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
    self._polarroot = ap._pyPolar2Polar(self.airfoils[1][2].pyPolar;
                                            x=root_x, y=root_y)
    self._polartip = ap._pyPolar2Polar(self.airfoils[end][2].pyPolar;
                                            x=tip_x, y=tip_y)
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
  fun_upper = gt.parameterize(upper[1], upper[2], fill(0.0, size(upper[1])); inj_var=1)
  fun_lower = gt.parameterize(lower[1], lower[2], fill(0.0, size(lower[1])); inj_var=1)

  # New discretization for both surfaces
  upper_points = gt.discretize(fun_upper, 0, 1, n_upper, r[1]; central=central)
  lower_points = gt.discretize(fun_lower, 0, 1, n_lower, r[1]; central=central)

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

"""Calculates the load distribution by using the airfoil lookup table on the
given inflow (this assumes that the inflow already includes all induced velocity
and it is the effective inflow).
"""
function _calc_distributedloads_lookuptable(ccbrotor::OCCBRotor,
                                            ccbinflow::OCCBInflow,
                                            turbine_flag::Bool;
                                            tiploss_correction::Bool=false)

  # check if propeller
  swapsign = turbine_flag ? 1 : -1

  # initialize arrays
  n = length(ccbrotor.r)
  Np = fill(0.0, n)
  Tp = fill(0.0, n)
  cn = fill(0.0, n)
  ct = fill(0.0, n)
  cl = fill(0.0, n)
  cd = fill(0.0, n)
  uvec = fill(0.0, n)
  vvec = fill(0.0, n)
  gamma = fill(0.0, n)

  for i in 1:n
    twist = swapsign*ccbrotor.theta[i]
    Vx = swapsign*ccbinflow.Vx[i]
    Vy = ccbinflow.Vy[i]

    aux1 = 0.5*ccbinflow.rho*(Vx*Vx+Vy*Vy)*ccbrotor.chord[i]

    # Effective angle of attack (rad)
    thetaV = atan(Vx, Vy)
    thetaeff = thetaV - twist
    # println("angles = $([twist, thetaeff]*180/pi)\tVx,Vy=$([Vx, Vy])")

    # airfoil cl/cd
    cl[i], cd[i] = occb_airfoil(ccbrotor.af[i], thetaeff)

    # Tip and hub correction factor
    if tiploss_correction
        B, Rtip, Rhub, r = ccbrotor.B, ccbrotor.Rtip, ccbrotor.Rhub, ccbrotor.r[i]

        # factortip = B/2.0*(Rtip - r)/(r*abs(thetaV))
        # factortip = B/2.0*(Rtip - r)/(r*abs(sin(thetaV)))
        # factortip = B/2.0*((Rtip-Rhub)/(r-Rhub) - 1.0)/abs(sin(max(10.0*pi/180, thetaV)))
        # Ftip = 2.0/pi*acos(exp(-factortip))

        # factortip = B/2.0*((Rtip-0.8*Rtip)/(r-0.8*Rtip) - 1.0)/abs(sin(thetaV))
        # Ftip = 2.0/pi*acos(exp(-abs(factortip)))

        # factortip = B/2.0*((Rtip-Rhub)/(r-Rhub) - 1.0)/abs(sin(thetaV))
        # Ftip = 2.0/pi*acos(exp(-factortip))
        # Ftip = 2.0/pi*acos(exp(-factortip^2))
        # Ftip = 2.0/pi*acos(exp(-factortip^10))

        # i == 1 ? println("i\tr/Rtip\tthetaV*180/pi\t(Rtip-Rhub)/(r-Rhub) - 1.0\tabs(sin(thetaV))\tfactortip\tFtip") : nothing
        # println("$i\t$(round(r/Rtip, digits=3))\t$(round(thetaV*180/pi, digits=3))\t$(round((Rtip-Rhub)/(r-Rhub) - 1.0, digits=3))"*
        #         "\t$(round(abs(sin(thetaV)), digits=3))\t$(round(factortip, digits=3))\t$(round(Ftip, digits=3))")

        # factorhub = B/2.0*(r - Rhub)/(Rhub*abs(thetaV))
        # factorhub = B/2.0*(r/Rhub - 1)/abs(sin(max(10.0*pi/180, thetaV)))
        # Fhub = 2.0/pi*acos(exp(-factorhub))

        factortip = B/2.0*(  (Rtip/r)^2 - 1  )/abs(sin(thetaV))^0.25
        Ftip = 2.0/pi*acos(exp(-factortip))

        factorhub = B/2.0*(   (r/Rhub)^0.6 - 1   )^5/abs(sin(thetaV))^0.5
        Fhub = 2.0/pi*acos(exp(-factorhub))

        F = Ftip * Fhub

        cl[i] *= F
        cd[i] *= F
    end

    # normal and tangential coefficients
    sthtV = sin(thetaV)
    cthtV = cos(thetaV)
    cn[i] = cl[i]*cthtV + cd[i]*sthtV
    ct[i] = swapsign*(cl[i]*sthtV - cd[i]*cthtV)

    # Normal and tangential forces per unit length
    Np[i] = cn[i]*aux1
    Tp[i] = ct[i]*aux1

    # Circulation
    gamma[i] = cl[i]*sqrt(Vx*Vx+Vy*Vy)*ccbrotor.chord[i]/2
  end

  # reverse sign of outputs for propellers
  # Tp *= swapsign # I already reversed this
  vvec *= swapsign

  return Np, Tp, uvec, vvec, gamma, cn, ct, cl, cd
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


function Base.deepcopy_internal(x::Rotor, stackdict::IdDict)
    if haskey(stackdict, x)
        return stackdict[x]
    end

    y = Rotor(  Base.deepcopy_internal(x.CW, stackdict),
                Base.deepcopy_internal(x.r, stackdict),
                Base.deepcopy_internal(x.chord, stackdict),
                Base.deepcopy_internal(x.theta, stackdict),
                Base.deepcopy_internal(x.LE_x, stackdict),
                Base.deepcopy_internal(x.LE_z, stackdict),
                Base.deepcopy_internal(x.B, stackdict),
                Base.deepcopy_internal(x.airfoils, stackdict),
                Base.deepcopy_internal(x.turbine_flag, stackdict),
                Base.deepcopy_internal(x.RPM, stackdict),
                Base.deepcopy_internal(x.hubR, stackdict),
                Base.deepcopy_internal(x.rotorR, stackdict),
                Base.deepcopy_internal(x.m, stackdict),
                Base.deepcopy_internal(x.sol, stackdict),
                Base.deepcopy_internal(x._wingsystem, stackdict),
                Base.deepcopy_internal(x._r, stackdict),
                Base.deepcopy_internal(x._chord, stackdict),
                Base.deepcopy_internal(x._theta, stackdict),
                Base.deepcopy_internal(x._LE_x, stackdict),
                Base.deepcopy_internal(x._LE_z, stackdict),
                Base.deepcopy_internal(x._polars, stackdict),
                Base.deepcopy_internal(x._polarroot, stackdict),
                Base.deepcopy_internal(x._polartip, stackdict))

    stackdict[x] = y
    return y
end
##### END OF ROTOR CLASS #######################################################












#
