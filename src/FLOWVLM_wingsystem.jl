# Class and methods for modeling interactions between multiple lifting surfaces

################################################################################
# WINGSYSTEM CLASS
################################################################################
"""
  `WingSystem()`
Initiates a system of wings. All method applicable to a Wing object are
applicable to a WingSystem. When solved, it will calculate the interaction
between wings.
"""
type WingSystem
  # Properties
  wings::typeof([])             # Wings in the system
  wing_names::typeof(String[])      # Names of the wings
  O::Array{Float64,1}               # Origin of local reference frame
  Oaxis::Array{Float64,2}           # Unit vectors of the local reference frame
  invOaxis::Array{Float64,2}        # Inverse unit vectors
  Vinf::Any                         # Vinf function used in current solution

  # Data storage
  sol::typeof(Dict())               # Solution fields available

  WingSystem( wings=[], wing_names=String[],
                O=[0.0,0.0,0.0],
                Oaxis=[1.0 0 0; 0 1 0; 0 0 1],
                invOaxis=[1.0 0 0; 0 1 0; 0 0 1],
                Vinf=nothing,
              sol=Dict()
      ) = new(wings, wing_names,
                O, Oaxis, invOaxis, Vinf,
              sol)
end

"""
  `addwing(self::WingSystem, wing_name::String, wing::Wing)`
Adds a wing to the system with the position and orientation of local reference
being interpreted in relation to the local reference frame of the system.
"""
function addwing(self::WingSystem, wing_name::String, wing;
                  overwrite=false)
  # Error case
  if wing_name in self.wing_names
    if overwrite
      ind = findfirst(self.wing_names, wing_name)
      deleteat!(self.wing_names, ind)
      deleteat!(self.wings, ind)
    else
      error("Wing '$wing_name' already exists in the system."*
            "Use overwrite option.")
    end
  end

  # Adds it
  push!(self.wings,  wing)
  push!(self.wing_names,  wing_name)

  # Interprets its reference frame in relation to the system's frame
  new_O, new_Oaxis = _interpret(wing.O, wing.Oaxis, self.O, self.invOaxis)
  setcoordsystem(wing, new_O, new_Oaxis)

  _reset(self)
end

"""
  `setcoordsystem(self, O, Oaxis; check=true, wings=String[])`
Redefines the local coordinate system of the system, where `O` is the new origin
and `Oaxis` is the matrix [i; j; k] of unit vectors. It transforms the
coordinates of all the wings in the system accordingly.

To change the local coordinate system of a specific wing relative to the
system's coordinate system, give the name of the wing in an array under argument
`wings`.
"""
function setcoordsystem(self::WingSystem, O::Array{Float64,1},
                            Oaxis::Array{Float64,2};
                            check=true, wings::Array{String,1}=String[])

  if check; check_coord_sys(Oaxis); end;

  # ----- Case of specific wings
  if size(wings)[1]!=0
    for wing_name in wings
      wing = get_wing(self, wing_name)
      setcoordsystem(wing, O, Oaxis; check=false)
      addwing(self, wing_name, wing; overwrite=true)
    end
    _reset(self; keep_Vinf=true)
    return
  end

  # ----- Case of the entire system
  invOaxis = inv(Oaxis)

  # Transforms the local coordinate system of each wing
  for wing in self.wings
    # Current frame in global coordinates
    curr_O, curr_Oaxis = wing.O, wing.Oaxis
    # Current frame in the system's coordinates
    local_curr_O, local_curr_Oaxis = _counter_interpret(curr_O, curr_Oaxis,
                                                      self.O, self.Oaxis)
    # New frame in global coordinates
    new_O, new_Oaxis = _interpret(local_curr_O, local_curr_Oaxis, O, invOaxis)

    setcoordsystem(wing, new_O, new_Oaxis)
  end

  # Sets the new system's reference frame
  self.O = O
  self.Oaxis = Oaxis
  self.invOaxis = invOaxis

  _reset(self; keep_Vinf=true)
end

function setcoordsystem(self::WingSystem, O::Array{Float64,1},
                            Oaxis::Array{Array{Float64,1},1};
                            check=true)
  dims = 3
  M = zeros(dims, dims)
  for i in 1:dims
    M[i, :] = Oaxis[i]
  end
  setcoordsystem(self, O, M; check=check)
end

"Sets Vinf(X,t) as the freestream of all wings in the system"
function setVinf(self::WingSystem, Vinf)
  _reset(self)
  self.Vinf = Vinf
  for wing in self.wings
    setVinf(wing, Vinf)
  end
end

"Returns the m-th control point of the system"
function getControlPoint(self::WingSystem, m::Int64)
  wing, _m = _fetch_wing(self, m)
  return getControlPoint(wing, _m)
end

"Returns the m-th horseshoe of the system in the global coordinate system"
function getHorseshoe(self::WingSystem, m::Int64; t::Float64=0.0, extraVinf...)
  wing, _m = _fetch_wing(self, m)
  return getHorseshoe(wing, _m; t=t, extraVinf...)
end

"Returns total number of lattices in the wing"
function get_m(self::WingSystem)
  m = 0
  for wing in self.wings
    m += get_m(wing)
  end
  return m
end

"Returns the wing in the system under the requested name"
function get_wing(self::WingSystem, wing_name::String)
  wing_i = findfirst(self.wing_names, wing_name)
  return get_wing(self, wing_i)
end

function get_wing(self::WingSystem, wing_names::Array{String,1})
  wings = []
  for wing_name in wing_names
    push!(wings, get_wing(self, wing_name))
  end
  return wings
end

"Returns the i-th wing in the system"
function get_wing(self::WingSystem, wing_i::Int64)
  return self.wings[wing_i]
end

"Returns the wings in the system as a dictionary Dict(wing_name => Wing)"
function get_wings(self::WingSystem)
  out = Dict()
  for (i, wing) in enumerate(self.wings)
    out[self.wing_names[i]] = wing
  end
  return out
end

##### CALCULATIONS #############################################################
"For a coordinate system 'inception2' that is incapsulated inside another
coordinate system 'inception1', it returns its interpretation in the global
coordinate system"
function _interpret(O2::Array{Float64,1}, Oaxis2::Array{Float64,2},
                    O1::Array{Float64,1}, invOaxis1::Array{Float64,2})
  new_O = countertransform(O2, invOaxis1, O1)
  new_Oaxis = zeros(3,3)
  for i in 1:3
    unit = Oaxis2[i, :]
    new_unit = countertransform(unit, invOaxis1, [0.0, 0, 0])
    new_Oaxis[i, :] = new_unit
  end
  return new_O, new_Oaxis
end

"For a coordinate system 'inception2' that is incapsulated inside another
coordinate system 'inception1', it receives its interpretation in the global
coordinate system and counterinterprets it back to the 'inception1' system."
function _counter_interpret(O2::Array{Float64,1}, Oaxis2::Array{Float64,2},
                    O1::Array{Float64,1}, Oaxis1::Array{Float64,2})
  new_O = transform(O2, Oaxis1, O1)
  new_Oaxis = zeros(3,3)
  for i in 1:3
    unit = Oaxis2[i, :]
    new_unit = transform(unit, Oaxis1, [0.0, 0, 0])
    new_Oaxis[i, :] = new_unit
  end
  return new_O, new_Oaxis
end

##### INTERNAL FUNCTIONS #######################################################
function _reset(self::WingSystem; verbose=false, keep_Vinf=false)
  if verbose; println("Resetting WingSystem"); end;
  self.sol = Dict()
  if keep_Vinf==false; self.Vinf = nothing; end;

  for wing in self.wings
    _reset(wing; verbose=verbose, keep_Vinf=keep_Vinf)
  end
end

function _addsolution(self::WingSystem, field_name::String, sol_field; t::Float64=0.0)
  self.sol[field_name] = sol_field
  prev_m = 0
  for wing in self.wings
    this_m = get_m(wing)
    _addsolution(wing, field_name, sol_field[ (prev_m+1) : prev_m+this_m ]; t=t)
    prev_m += this_m
  end
end

"Given a panel index in the system, it returns the wing the panel belongs to
and that wing's index of the panel"
function _fetch_wing(self::WingSystem, m::Int64)
  if m<=0 || m>get_m(self)
    error("Requested unexistent panel (m<=0 || m>get_m(self))")
  end

  _m = m
  # Finds the wing the panel belongs to
  wing = nothing
  for wing_i in self.wings
    if _m<=get_m(wing_i)
      wing = wing_i
      break
    else
      _m -= get_m(wing_i)
    end
  end

  if wing==nothing
    error("CRITICAL ERROR: Logic error! (wing==nothing)")
  end

  return wing, _m
end
##### END OF CLASS #############################################################
