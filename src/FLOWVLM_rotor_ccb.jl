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

# ------------- More robust rotation correction and extrapolation functions ----
# """
#     fitliftslope(aoa,cl,tol::Real=0.05,allnegfit::Bool=false)
# Returns the lift slope and zero lift angle of attack as well as the data used
# to compute these parameters.  The tolerance `tol` (in terms of cl) for the fit
# and whether to include all data points below zero angle of attack `allnegfit`
# may be specified. Returns liftslope,zeroliftangle,aoafit,clfit
# """
# function fitliftslope(aoa,cl,tol::Real=0.05,allnegfit::Bool=false,center=0.0)

#   if length(aoa) != length(cl)
#     error("aoa and cl must have the same length")
#   end

#   # split into positive and negative angles of attack
#   idxcenter = argmin(abs.(aoa .- center))
#   idxneg = findall(aoa .< aoa[idxcenter])
#   sort!(idxneg,rev = true)
#   idxpos = findall(aoa .> aoa[idxcenter])
#   sort!(idxpos)

#   # initialize arrays of data to fit
#   aoafit = [aoa[idxcenter]]
#   clfit = [cl[idxcenter]]
#   # create indexes for increasing and decreasing aoa
#   ipos = 1
#   ineg = 1
#   # initialize flag to stop adding data > 0
#   upperflag = false
#   if ipos > length(idxpos)
#     upperflag = true
#   end
#   # initialize flag to stop adding data < 0
#   lowerflag = false
#   if ineg > length(idxneg)
#     lowerflag = true
#   end
#   # compute lift slope and zero lift angle of attack
#   success = false
#   liftslope = 2*pi
#   zeroliftangle = 0.0
#   for i = 1:(length(idxpos)+length(idxneg))
#     # add positive angle of attack point to fit if closest to zero aoa
#     if upperflag == false
#       if lowerflag == true || abs(aoa[idxpos[ipos]]) <= abs(aoa[idxneg[ineg]])
#         # add point
#         push!(aoafit,aoa[idxpos[ipos]])
#         push!(clfit,cl[idxpos[ipos]])
#         # order points
#         order = sortperm(aoafit)
#         aoafit = aoafit[order]
#         clfit = clfit[order]
#         # create least squares fit
#         lsqsol = hcat(aoafit, -ones(length(aoafit)))\clfit
#         # if liftslope not defined from previous iteration, assign it here
#         if ipos == 1
#             liftslope = lsqsol[1]
#         end
#         # get slope and error of added point
#         yerror = lsqsol[1]*(aoa[idxpos[ipos:end]].-lsqsol[2]/lsqsol[1]).-cl[idxpos[ipos:end]]
#         # slope = (clfit[end] - clfit[end-1])/(aoafit[end] - aoafit[end-1])
#         # check if on linear portion of lift curve
#         if any(yerror .< tol) #slope > 0.5*liftslope ||
#           liftslope = lsqsol[1]
#           zeroliftangle = lsqsol[2]/lsqsol[1]
#           success = true
#         else # if not on linear portion remove from fit
#           aoafit = aoafit[1:(end-1)]
#           clfit = clfit[1:(end-1)]
#           upperflag = true
#         end
#         #  increment counter and check for more points greater than aoa=0
#         ipos = ipos+1
#         if ipos > length(idxpos)
#           upperflag = true
#         end
#       end
#     end
#     # add negative angle of attack point to fit if closest to zero
#     if lowerflag == false
#       if upperflag == true || abs(aoa[idxneg[ineg]]) <= abs(aoa[idxpos[ipos]])
#         # add point
#         pushfirst!(aoafit,aoa[idxneg[ineg]])
#         pushfirst!(clfit,cl[idxneg[ineg]])
#         # create least squares fit
#         lsqsol = hcat(aoafit, -ones(length(aoafit)))\clfit
#         # get error of added point
#         yerror = lsqsol[1]*(aoa[idxneg[ineg:end]].-lsqsol[2]/lsqsol[1]).-cl[idxneg[ineg:end]]
#         # check if on linear portion of lift curve
#         if any(abs.(yerror) .< tol) || allnegfit
#           liftslope = lsqsol[1]
#           zeroliftangle = lsqsol[2]/lsqsol[1]
#           success = true
#         else # if not on linear portion remove from fit
#           aoafit = aoafit[2:end]
#           clfit = clfit[2:end]
#           lowerflag = true
#         end
#         #  increment counter and check for more points less than aoa=0
#         ineg = ineg+1
#         if ineg > length(idxneg)
#           lowerflag = true
#         end
#       end
#     end
#     # break if no more points to add
#     if upperflag && lowerflag
#       break
#     end
#   end
#   if success == false
#     liftslope,zeroliftangle,aoafit,clfit = fitliftslope(aoa,cl,tol,allnegfit,aoa[idxcenter]+1.0*pi/180)
#   end

#   return liftslope,zeroliftangle,aoafit,clfit
# end


# """
#     correction3D(cl, cd, cr, rR, tsr, alpha, phi=alpha, alpha0, alpha_max_corr=30*pi/180)

# Apply rotational corrections (3D stall delay) to airfoil data using DUSelig (lift) and Eggers (drag) models.

# **Arguments**
# - `cl::Array{Float64}`: lift coefficient before correction
# - `cd::Array{Float64}`: drag coefficient before correction
# - `cr::Float64`: local chord / local radius
# - `rR::Float64`: local radius / tip radius
# - `tsr::Float64`: local tip speed ratio (Omega r / Vinf)
# - `alpha::Array{Float64}`: local angle of attack (radians)
# - `phi::Float64`: local inflow angles (defaults to angle of attack is precomputing since it is only known for on-the-fly computations)
# - `alpha0::Float64`: zero-lift angle of attack (radians)
# - `alpha_max_corr::Float64`: angle of attack for maximum correction (tapers off to zero by 90 degrees) (radians)

# **Returns**
# - `cl::Float64`: lift coefficient after correction
# - `cd::Float64`: drag coefficient after correction
# """
# function correction3D(cl, cd, cr, rR, tsr, alpha, phi=alpha, alpha_max_corr=30*pi/180)

#     # if cd[1] == cd[2] == cd[3]
#     #     return cl, cd
#     # end

#     m, alpha0 ,_ ,_ = fitliftslope(alpha,cl)

#     for i=1:length(alpha)
#         # Du-Selig correction for lift
#         Lambda = tsr / sqrt(1 + tsr^2)
#         expon = 1.0 / (Lambda * rR)
#         fcl = 1.0/m*(1.6*cr/0.1267*(1.0-cr^expon)/(1.0+cr^expon)-1)

#         # linear lift line
#         cl_linear = m*(alpha[i] - alpha0)

#         # adjustment for max correction
#         amax = atan(1/0.12) - 5*pi/180  # account for singularity in Eggers (not pi/2)
#         if abs(alpha[i]) >= amax
#             adj = 0.0
#         elseif abs(alpha[i]) > alpha_max_corr
#             adj = ((amax-abs(alpha[i]))/(amax-alpha_max_corr))^2
#         else
#             adj = 1.0
#         end

#         # increment in cl
#         deltacl = adj*fcl*(cl_linear - cl[i])
#         cl[i] += deltacl

#         # Eggers correction for drag
#         deltacd = deltacl * (sin(phi[i]) - 0.12*cos(phi[i]))/(cos(phi[i]) + 0.12*sin(phi[i]))  # note that we can actually use phi instead of alpha as is done in airfoilprep.py b/c this is done at each iteration
#         cd[i] += deltacd
#     end

#         return cl, cd

# end


# """
#     extrapolate(alpha, cl, cd, cr75, nalpha=50)

# Viterna extrapolation.  Follows Viterna paper and somewhat follows NREL version of AirfoilPrep, but with some modifications for better robustness and smoothness.

# **Arguments**
# - `alpha::Vector{Float64}`: angles of attack (radians)
# - `cl::Vector{Float64}`: correspnding lift coefficients
# - `cd::Vector{Float64}`: correspnding drag coefficients
# - `cr75::Float64`: chord/Rtip at 75% Rtip
# - `nalpha::Int64`: number of discrete points (angles of attack) to include in extrapolation

# **Returns**
# - `alpha::Vector{Float64}`: angle of attack from -pi to pi (radians)
# - `cl::Vector{Float64}`: correspnding extrapolated lift coefficients
# - `cd::Vector{Float64}`: correspnding extrapolated drag coefficients
# """
# function extrapolate(alpha, cl, cd, cr75, nalpha=60)

#     # if isapprox(alpha[1],-180.0) || isapprox(alpha[1],-pi)
#     #     return alpha, cl, cd
#     # end

#     # estimate cdmax
#     AR = 1.0 / cr75
#     cdmaxAR = 1.11 + 0.018*AR
#     cdmax = max(maximum(cd), cdmaxAR)

#     # find clmax
#     i_ps = argmax(cl)  # positive stall
#     cl_ps = cl[i_ps]
#     cd_ps = cd[i_ps]
#     a_ps = alpha[i_ps]

#     # and clmin
#     i_bs = alpha .< a_ps  # before stall
#     i_ns = argmin(cl[i_bs])  # negative stall
#     cl_ns = cl[i_bs][i_ns]
#     cd_ns = cd[i_bs][i_ns]
#     a_ns = alpha[i_bs][i_ns]

#     # coefficients in method
#     B1pos = cdmax
#     A1pos = B1pos/2.0 * ones(nalpha)
#     sa = sin(a_ps)
#     ca = cos(a_ps)
#     A2pos = (cl_ps - cdmax*sa*ca)*sa/ca^2
#     B2pos = (cd_ps - cdmax*sa^2)/ca * ones(nalpha)

#     B1neg = cdmax
#     A1neg = B1neg/2.0
#     sa = sin(a_ns)
#     ca = cos(a_ns)
#     A2neg = (cl_ns - cdmax*sa*ca)*sa/ca^2 * ones(nalpha)
#     B2neg = (cd_ns - cdmax*sa^2)/ca * ones(nalpha)

#     # angles of attack to extrapolate to
#     apos = range(alpha[end], pi, length=nalpha+1)
#     apos = apos[2:end]  # don't duplicate point
#     aneg = range(-pi, alpha[1], length=nalpha+1)
#     aneg = aneg[1:end-1]  # don't duplicate point

#     # high aoa adjustments
#     adjpos = ones(nalpha)
#     idx = findall(apos .>= pi/2)
#     adjpos[idx] .= -0.7
#     A1pos[idx] .*= -1
#     B2pos[idx] .*= -1

#     # idx = findall(aneg .<= -alpha[end])

#     adjneg = ones(nalpha)
#     idx = findall(aneg .<= -pi/2)
#     adjneg[idx] .= 0.7
#     A2neg[idx] .*= -1
#     B2neg[idx] .*= -1

#     # extrapolate
#     clpos = @. adjpos * (A1pos*sin(2*apos) + A2pos*cos(apos)^2/sin(apos))
#     cdpos = @. B1pos*sin(apos)^2 + B2pos*cos(apos)
#     clneg = @. adjneg * (A1neg*sin(2*aneg) + A2neg*cos(aneg)^2/sin(aneg))
#     cdneg = @. B1neg*sin(aneg)^2 + B2neg*cos(aneg)

#     # # override region between -alpha_high and alpha_low (if it exists)
#     # idx = findall(-alpha[end] .<= aneg .<= alpha[1])
#     # @. clneg[idx] = -cl[end]*0.7 + (aneg[idx]+alpha[end])/(alpha[1]+alpha[end])*(cl[1]+cl[end]*0.7)
#     # @. cdneg[idx] = cd[1] + (aneg[idx]-alpha[1])/(-alpha[end]-alpha[1])*(cd[end]-cd[1])


#     # override with linear variation at ends
#     idx = findall(apos .>= pi-a_ps)
#     @. clpos[idx] = (apos[idx] - pi)/a_ps*cl_ps*0.7
#     idx = findall(aneg .<= -pi-a_ns)
#     @. clneg[idx] = (aneg[idx] + pi)/a_ns*cl_ns*0.7

#     # concatenate
#     alphafull = [aneg; alpha; apos]
#     clfull = [clneg; cl; clpos]
#     cdfull = [cdneg; cd; cdpos]

#     # don't allow negative drag
#     cdfull = max.(cdfull, 0.0001)
#     return alphafull, clfull, cdfull
# end
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
                          sound_spd=nothing, AR_to_360extrap=true, CDmax = 1.3)


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
    if sound_spd!==nothing
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

    #OLD
    # #3D corrections
    # this_polar = ap.correction3D(this_polar, r_over_R, c_over_r, tsr)

    # # 360 extrapolation
    # if AR_to_360extrap
    #     # use a linear interpolation instead of Splines; then don't expose this flag in FLOWUnsteady
    #     c_spline1D = Spline1D(self._r / self.rotorR, self._chord; k=1)
    #     c_75 = c_spline1D(0.75)
    #     AR = c_75 / self.rotorR
    #     this_polar = ap.extrapolate(this_polar, CDmax, AR=AR)
    # else
    #     this_polar = ap.extrapolate(this_polar, CDmax)
    # end


    #NEW
    #360 extrapolation
    #get coeffs, note: alpha comes out in degrees
    alpha, cl = ap.get_cl(this_polar)

    cd = ap.get_cd(this_polar)[2]

    c_spline1D = Spline1D(self._r / self.rotorR, self._chord; k=1)
    c_75 = c_spline1D(0.75)
    #run extrapolation
    #note: convert alpha to radians on input, output in radians
    if !isapprox(alpha[1],-180.0)
        #convert to radians for corrections
        alpha *= pi/180

        #extrapolate
        alpha, cl, cd = ccb.viterna(alpha, cl, cd, c_75)
        #run 3D corrections
        #note: alpha still in radians
        if tsr !== nothing
            #rotation corrections
            cl, cd = ccb.rotation_correction(ccb.DuSeligEggers(),cl, cd, c_over_r, r_over_R, tsr, alpha)
        end
        #convert alpha back to degrees
        alpha *= 180/pi
    end
    #reconstruct polar for injective function
    #note: alpha back in degrees
    this_polar = ap.Polar(ap.get_Re(this_polar), alpha, cl, cd, zeros(length(alpha)); ap._get_nonpypolar_args(this_polar)...)
    # Makes sure the polar is injective for easing the spline
    this_polar = ap.injective(this_polar)

    # Converts to CCBlade's AirfoilData object
    #note: alpha back in degrees
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
