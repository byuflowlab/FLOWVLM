# ------------- OLD CCBlade TYPES ----------------------------------------------

# pretabulated cl/cd data
struct OCCBAirfoilData
    cl::Dierckx.Spline1D
    cd::Dierckx.Spline1D
end

"""
    Rotor(r, chord, theta, af, Rhub, Rtip, B, precone)

Define rotor geometry.

**Arguments**
- `r::Array{Float64, 1}`: radial locations (m)
- `chord::Array{Float64, 1}`: chord lengths (m)
- `theta::Array{Float64, 1}`: total twist including pitch (rad)
- `af::Array{AirfoilData, 1}`: airfoils
- `Rhub::Float64`: hub radius (along blade length)
- `Rtip::Float64`: tip radius (along blade length)
- `B::Int64`: number of blades
- `precone::Float64`: precone angle (rad)
"""
struct OCCBRotor
    r#::Array{Float64, 1}
    chord#::Array{Float64, 1}
    theta#::Array{Float64, 1}
    af#::Array{AirfoilData, 1}
    Rhub#::Float64
    Rtip#::Float64
    B#::Int64
    precone#::Float64
end

# operating point for the turbine/propeller
struct OCCBInflow
    Vx#::Array{Float64, 1}
    Vy#::Array{Float64, 1}
    rho#::Float64
end

"""
    af_from_data(alpha, cl, cd)

Create an AirfoilData object directly from alpha, cl, and cd arrays.

af_from_aerodynfile calls this function indirectly.  Uses a cubic B-spline
(if the order of the data permits it).  A small amount of smoothing of
lift and drag coefficients is also applied to aid performance
for gradient-based optimization.
"""
function occb_af_from_data(alpha, cl, cd; spl_k=3)

    k = min(length(alpha)-1, spl_k)  # can't use cubic spline if number of entries in alpha is small

    # 1D interpolations for now.  ignoring Re dependence (which is very minor)
    afcl = Dierckx.Spline1D(alpha*pi/180.0, cl; k=k, s=0.1)
    afcd = Dierckx.Spline1D(alpha*pi/180.0, cd; k=k, s=0.001)
    af = OCCBAirfoilData(afcl, afcd)

    return af
end

"""
private

evalute airfoil spline at alpha
"""
function occb_airfoil(af::OCCBAirfoilData, alpha::Real)

    cl = af.cl(alpha)
    cd = af.cd(alpha)

    return cl, cd
end

# ------------- FLOWVLM <-> CCBlade FUNCTIONS ----------------------------------
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
function FLOWVLM2OCCBlade(self,#::Rotor,
                          RPM, blade_i::IWrap, turbine_flag::Bool;
                          sound_spd=nothing, AR_to_360extrap=false, CDmax = 1.3)


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
  af = OCCBAirfoilData[]
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
      if Ma>=1
        error("Mach correction requested on Ma = $Ma >= 1.0")
      end
      alpha, cl = ap.get_cl(polar)
      this_polar = ap.Polar(ap.get_Re(polar), alpha, cl/sqrt(1-Ma^2),
                              ap.get_cd(polar)[2], ap.get_cm(polar)[2];
                                            ap._get_nonpypolar_args(polar)...)
    else
      this_polar = polar
    end

    # 3D corrections
    this_polar = ap.correction3D(this_polar, r_over_R, c_over_r, tsr)

    # 360 extrapolation
    if AR_to_360extrap
        # use a linear interpolation instead of Splines; then don't expose this flag in FLOWUnsteady
        c_spline1D = Spline1D(self._r / self.rotorR, self._chord; k=1)
        c_75 = c_spline1D(0.75)
        AR = c_75 / self.rotorR
        this_polar = ap.extrapolate(this_polar, CDmax, AR=AR)
    else
        this_polar = ap.extrapolate(this_polar, CDmax)
    end

    # Makes sure the polar is injective for easing the spline
    this_polar = ap.injective(this_polar)

    # Converts to CCBlade's AirfoilData object
    alpha, cl = ap.get_cl(this_polar)
    _, cd = ap.get_cd(this_polar)
    ccb_polar = occb_af_from_data(alpha, cl, cd; spl_k=5)

    push!(af, ccb_polar)
  end

  rotor = OCCBRotor(self._r, self._chord,
                    (-1)^(turbine_flag)*(-1)^self.CW*self._theta*pi/180, af,
                    Rhub, Rtip, self.B, precone)

  return rotor
end

"""
Convert a rotor object from the old-CCBlade format to the new CCBlade type.
"""
function OCCB2CCB(orotor::OCCBRotor, turbine::Bool, oinflow::OCCBInflow;
                                                                pitch::Real=0.0)
    rotor = ccb.Rotor(orotor.Rhub, orotor.Rtip, orotor.B,
                                    turbine, pitch, orotor.precone)

    airfoil_funs = [(alpha, Re, Mach) -> occb_airfoil(af, alpha) for af in orotor.af]

    sections = ccb.Section.(orotor.r, orotor.chord, orotor.theta, airfoil_funs)

    ops = ccb.OperatingPoint.(oinflow.Vx, oinflow.Vy, oinflow.rho)

    return rotor, sections, ops
end


"""Receives a vector in the global coordinate system and transforms it into
CCBlade's coordinate system relative to `blade`. NOTE: This function only
rotates `V` into the new axis without translating it unless otherwise indicated.
(for definition of axes see notebook entry 20171202)
"""
function _global2ccblade(blade::Wing, V::FArrWrap, CW::Bool;
                                                        translate::Bool=false)
  # V in FLOWVLM Rotor's blade c.s.
  V_vlm = transform(V, blade.Oaxis, translate ? blade.O : fill(0.0, 3))

  # CCBlade c.s. transformation matrix
  ccb_Oaxis = _ccbladeOaxis(blade, CW)

  # V in CCBlade's c.s.
  V_ccb = transform(V_vlm, ccb_Oaxis, fill(0.0, 3))

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
  V_vlm = countertransform(V, inv(ccb_Oaxis), fill(0.0, 3))

  # V in global c.s.
  V_glob = countertransform(V_vlm, blade.invOaxis, translate ? blade.O : fill(0.0, 3))

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
