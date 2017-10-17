
using Dierckx

################################################################################
# ROTOR CLASS
################################################################################

type Rotor
  # NOTE TO SELF: r is the y-direction on a wing, hence, remember to build the
  #               blade from root in the direction of positive y.

  # Initialization variables (USER INPUT)
  CW::Bool                      # True for clockwise rotation
  r::Array{Float64,1}           # Radius position for the following variables
  chord::Array{Float64,1}       # Chord length
  theta::Array{Float64,1}       # Angle of attack from the rotor's axis
  LE_x::Array{Float64,1}        # x-position of leading edge
  LE_z::Array{Float64,1}        # z-position of leading edge (Height from plane
                                #                                  of rotation)
  B::Int64                      # Number of blades

  # Properties
  hubR::Float64                 # Hub radius
  rotorR::Float64               # Rotor radius

  # Data storage
  _wingsystem::WingSystem       # Rotor assembly

  Rotor(
          CW, r, chord, theta, LE_x, LE_z, B,
          hubR=r[1], rotorR=r[end],
          _wingsystem=WingSystem()
        ) = new(
          CW, r, chord, theta, LE_x, LE_z, B,
          hubR, rotorR,
          _wingsystem
        )
end

"Initializes the propeller "
function initialize(self::Rotor, n::Int64)
  # Checks for arguments consistency
  _check(self)

  # Generates blade
  blade = _generate_blade(self, n)

  # Generates full rotor
  init_angle = 0
  d_angle = 2*pi/self.B
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
    this_Oaxis = self.CW ? [0 0 -1; 0 1 0; 1 0 0] : [0 0 1; 0 1 0; -1 0 0]
    this_Oaxis = this_Oaxis*[
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
end

function setcoordsystem(self::Rotor, args...)
  setcoordsystem(self._wingsystem, args...)
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

##### INTERNAL FUNCTIONS #############################################################
"Checks for consistency in internal variables"
function _check(self::Rotor)
  # Matching dimensions
  nr = size(self.r)[1]
  for (label, arr) in [("chord",self.chord), ("theta",self.theta),
                        ("LE_x",self.LE_x), ("LE_z",self.LE_z)]
    if size(arr)[1] != nr
      error("Invalid dimensions in array $label. Expected $nr, found $(size(arr)[1]).")
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
  _spl_chord = Dierckx.Spline1D(self.r, self.chord;
                                  k=spline_k,s=spline_s)
  _spl_theta = Dierckx.Spline1D(self.r, -self.theta;
                                  k=spline_k,s=spline_s)
  _spl_LE_x = Dierckx.Spline1D(self.r, self.LE_x;
                                  k=spline_k,s=spline_s)
  _spl_LE_z = Dierckx.Spline1D(self.r, self.LE_z;
                                  k=spline_k,s=spline_s)
  spl_chord(x) = Dierckx.evaluate(_spl_chord, x)
  spl_theta(x) = (-1)^(self.CW==false)*Dierckx.evaluate(_spl_theta, x)
  spl_LE_x(x) = Dierckx.evaluate(_spl_LE_x, x)
  spl_LE_z(x) = Dierckx.evaluate(_spl_LE_z, x)


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
  end

  # Verifies correct discretization
  if abs((cumlen-l)/l)>0.0000000001
    error("Critical logic error! cumlen!=l ($cumlen!=$l)")
  end

  return blade
end













#
