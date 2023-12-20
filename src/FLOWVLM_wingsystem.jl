# Class and methods for modeling interactions between multiple lifting surfaces

################################################################################
# WINGSYSTEM CLASS
################################################################################
"""
    WingSystem(wings::Array{Union{Wing, WingSystem}}, wing_names::Array{String})

Initiates a system of wings. All methods applicable to a Wing object are
applicable to a WingSystem. When solved, it will calculate the interaction
between wings.
"""
mutable struct WingSystem{TF_design, TF_trajectory<:FWrap} <: AbstractWing{TF_design,TF_trajectory}
  # Properties
  wings::Array{Any,1}             # Wings in the system
  wing_names::Array{String,1}     # Names of the wings
  O::Vector{TF_trajectory}                     # Origin of local reference frame
  Oaxis::Matrix{TF_trajectory}                   # Unit vectors of the local reference frame
  invOaxis::Matrix{TF_trajectory}                # Inverse unit vectors
  Vinf::Union{Nothing,Function}                       # Vinf function used in current solution

  # Data storage
  sol::Dict{String, Any}          # Solution fields available
end

function WingSystem(; TF_design=Float64, TF_trajectory=Float64, wings=[], wing_names=String[],
    O=[0.0,0.0,0.0],
    Oaxis=[1.0 0 0; 0 1 0; 0 0 1],
    invOaxis=[1.0 0 0; 0 1 0; 0 0 1],
    Vinf=nothing,
  sol=Dict{String,Any}()
)
    O = TF_trajectory.(O)
    Oaxis = TF_trajectory.(Oaxis)
    invOaxis = TF_trajectory.(invOaxis)
    return WingSystem{TF_design,TF_trajectory}(wings, wing_names,
            O, Oaxis, invOaxis, Vinf, sol)
end

function WingSystem{TF_design,TF_trajectory}(; wings=[], wing_names=String[],
    O=[0.0,0.0,0.0],
    Oaxis=[1.0 0 0; 0 1 0; 0 0 1],
    invOaxis=[1.0 0 0; 0 1 0; 0 0 1],
    Vinf=nothing,
    sol=Dict{String,Any}()
) where {TF_design,TF_trajectory}
    O = TF_trajectory.(O)
    Oaxis = TF_trajectory.(Oaxis)
    invOaxis = TF_trajectory.(invOaxis)

    return WingSystem{TF_design,TF_trajectory}(wings, wing_names,
            O, Oaxis, invOaxis, Vinf, sol)
end

"""
    addwing(self::WingSystem, wing_name::String, wing::Union{Wing, Rotor})

Adds a wing (or rotor) to the system. The local reference frame of the wing
then is then in relation to the local reference frame of the System.
"""
function addwing(self::WingSystem, wing_name::String, wing;
                  overwrite=false, reset=true)
  # Error case
  if wing_name in self.wing_names
    if overwrite
      ind = findfirst(x->x==wing_name, self.wing_names)
      deleteat!(self.wing_names, ind)
      deleteat!(self.wings, ind)
    else
      error("Wing '$wing_name' already exists in the system."*
            " Use overwrite option.")
    end
  end

  # Adds it
  push!(self.wings,  wing)
  push!(self.wing_names,  wing_name)

  # Interprets its reference frame in relation to the system's frame
  new_O, new_Oaxis = _interpret(_get_O(wing), _get_Oaxis(wing), self.O, self.invOaxis)
  setcoordsystem(wing, new_O, new_Oaxis)

  if reset
      _reset(self)
  end
end

"""
    setcoordsystem(system::WingSystem, O::Vector, Oaxis::Matrix; wings=String[])

Redefines the local coordinate system of the system, where `O` is the new origin
and `Oaxis` is the matrix of unit vectors. It transforms the coordinates of all
wings in the system accordingly.

To change the local coordinate system of a specific wing relative to the
system's coordinate system, give the name of the wing in an array under argument
`wings`.
"""
function setcoordsystem(self::WingSystem{<:Any,TF_trajectory}, O::Vector{<:FWrap},
                            Oaxis::Matrix{<:FWrap};
                            check=true, wings::Array{String,1}=String[]) where TF_trajectory

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
    curr_O, curr_Oaxis = _get_O(wing), _get_Oaxis(wing)
    # Current frame in the system's coordinates
    local_curr_O, local_curr_Oaxis = _counter_interpret(curr_O, curr_Oaxis,
                                                      self.O, self.Oaxis)
    # New frame in global coordinates
    new_O, new_Oaxis = _interpret(local_curr_O, local_curr_Oaxis, O, invOaxis)

    setcoordsystem(wing, new_O, new_Oaxis)
  end

  # Sets the new system's reference frame
  self.O .= O
  self.Oaxis .= Oaxis
  self.invOaxis .= invOaxis

  _reset(self; keep_Vinf=true)
end

function setcoordsystem(self::WingSystem, O::Vector{<:FWrap},
                            Oaxis::Array{T,1} where {T<:AbstractArray};
                            check=true)
  dims = 3
  M = zeros(eltype(Oaxis), dims, dims)
  for i in 1:dims
    M[i, :] = Oaxis[i]
  end
  setcoordsystem(self, O, M; check=check)
end

"Sets Vinf(X,t) as the freestream of all wings in the system"
function setVinf(self::WingSystem, Vinf; keep_sol=false)
  _reset(self; keep_sol=keep_sol)
  self.Vinf = Vinf
  for wing in self.wings
    setVinf(wing, Vinf; keep_sol=keep_sol)
  end
end

"Returns the undisturbed freestream at each control point, or at the horseshoe
point indicated as `target`."
function getVinfs(self::WingSystem{TF_design,TF_trajectory}; t::FWrap=0.0, target="CP",
                              extraVinf=nothing, extraVinfArgs...) where {TF_design,TF_trajectory}

  TF_promoted = promote_type(TF_design,TF_trajectory,typeof(t))
  Vinfs = Vector{TF_promoted}[]
  for wing in self.wings
    for V in getVinfs(wing; t=t, target=target,
                                  extraVinf=extraVinf, extraVinfArgs...)
      push!(Vinfs, V)
    end
  end
  return Vinfs
end

"Returns the m-th control point of the system"
function getControlPoint(self::WingSystem, m::IWrap)
  wing, _m = _fetch_wing(self, m)
  return getControlPoint(wing, _m)
end

"Returns the m-th horseshoe of the system in the global coordinate system"
function getHorseshoe(self::WingSystem, m::IWrap; t::FWrap=0.0, extraVinf...)
  wing, _m = _fetch_wing(self, m)
  return getHorseshoe(wing, _m; t=t, extraVinf...)
end

"""
    get_m(system::WingSystem)

Returns the total number of horseshoes in the system
"""
function get_m(self::WingSystem)
  m = 0
  for wing in self.wings
    m += get_m(wing)
  end
  return m
end

"""
    get_wing(self::WingSystem, wing_name::String)

Returns the wing of name `wing_name`.
"""
function get_wing(self::WingSystem, wing_name::String)
  wing_i = findfirst(x->x==wing_name, self.wing_names)
  return get_wing(self, wing_i)
end

function get_wing(self::WingSystem, wing_names::Array{String,1})
  wings = []
  for wing_name in wing_names
    push!(wings, get_wing(self, wing_name))
  end
  return wings
end

"""
    get_wing(self::WingSystem, wing_i::Int)

Returns the i-th wing.
"""
function get_wing(self::WingSystem, wing_i::IWrap)
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
function _interpret(O2::Vector{<:FWrap}, Oaxis2::Matrix{<:FWrap},
                    O1::Vector{<:FWrap}, invOaxis1::Matrix{<:FWrap})
  new_O = countertransform(O2, invOaxis1, O1)
  new_Oaxis = zeros(eltype(Oaxis2), 3,3)
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
function _counter_interpret(O2::Vector{<:FWrap}, Oaxis2::Matrix{<:FWrap},
                    O1::Vector{<:FWrap}, Oaxis1::Matrix{<:FWrap})
  new_O = transform(O2, Oaxis1, O1)
  new_Oaxis = zeros(eltype(Oaxis2), 3,3)
  for i in 1:3
    unit = Oaxis2[i, :]
    new_unit = transform(unit, Oaxis1, [0.0, 0, 0])
    new_Oaxis[i, :] = new_unit
  end
  return new_O, new_Oaxis
end

##### INTERNAL FUNCTIONS #######################################################
function _reset(self::WingSystem; verbose=false, keep_Vinf=false, keep_sol=false)
  if verbose; println("Resetting WingSystem"); end;

  if keep_sol==false
      self.sol = Dict()
  else
      self.sol = Dict([entry for entry in self.sol if entry[1]!="Gamma" ])
  end
  if keep_Vinf==false; self.Vinf = nothing; end;

  for wing in self.wings
    _reset(wing; verbose=verbose, keep_Vinf=keep_Vinf, keep_sol=keep_sol)
  end
end

function _addsolution(self::WingSystem, field_name::String, sol_field; t::FWrap=0.0)
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
function _fetch_wing(self::WingSystem, m::IWrap)
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

function _get_O(wing)
  return wing.O
end
function _get_Oaxis(wing)
  return wing.Oaxis
end

function Base.deepcopy_internal(x::WingSystem{TF_design,TF_trajectory}, stackdict::IdDict) where {TF_design,TF_trajectory}
    if haskey(stackdict, x)
        return stackdict[x]
    end

    y = WingSystem{TF_design, TF_trajectory}( Base.deepcopy_internal(x.wings, stackdict),
                    Base.deepcopy_internal(x.wing_names, stackdict),
                    Base.deepcopy_internal(x.O, stackdict),
                    Base.deepcopy_internal(x.Oaxis, stackdict),
                    Base.deepcopy_internal(x.invOaxis, stackdict),
                    x.Vinf,
                    Base.deepcopy_internal(x.sol, stackdict))

    stackdict[x] = y
    return y
end
##### END OF CLASS #############################################################
