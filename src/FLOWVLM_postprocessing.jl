# All sort of functions for postprocessing that are not relevant for the solver
# or geometry definition (calculation of aerodynamic forces, fluid domain,
# visualization, etc).



################################################################################
# WING POSTPROCESSING
################################################################################
"""
  `calculate_field(wing, field_name [, qinf, S, l, r_cg])`
Calculates the field `field_name` on a Wing or WingSystem. See `FIELDS` for
a list of implemented fields.
"""
function calculate_field(wing, field_name::String;
                        rhoinf=nothing, qinf="automatic",
                        S="automatic", l="automatic", r_cg="automatic",
                        t::FWrap=0.0, lifting_interac::Bool=true)

  # --- ERROR CASES
  ## Unknown field
  if false==(field_name in keys(FIELDS))
    error("Invalid solution field '$field_name'."*
          " Valid fields: $(keys(FIELDS))")
  end
  ## Missing dependent field
  dependents, field_type = FIELDS[field_name]
  for dep in dependents
    if false==(dep in keys(wing.sol))
      error("Field '$field_name' requested, but '$dep' hasn't been solved yet")
    end
  end

  ######## AERODYNAMIC FORCES ###########################
  if field_name in ["Ftot", "L", "S", "D"]
    if rhoinf==nothing
      error("$(field_name) requested but rhoinf is missing")
    end

    F, SS, D, L = _calculate_forces(wing, rhoinf; t=t,
                                            lifting_interac=lifting_interac)
    _addsolution(wing, "Ftot", F)
    _addsolution(wing, "S", SS)
    _addsolution(wing, "D", D)
    _addsolution(wing, "L", L)

  ######## FORCE COEFFICIENTS ###########################
  elseif field_name in ["CFtot", "CL", "CS", "CD"]
    if rhoinf==nothing && qinf!="automatic"
      println("$(field_name) requested with a given qinf, but rhoinf is missing"
                *". Given qinf will be ignored.")
      Vinf = _aveVinf(wing; t=t)
      _rhoinf = 1.0
      _qinf = (1.0/2)*_rhoinf*dot(Vinf, Vinf)
    elseif qinf=="automatic"
      Vinf = _aveVinf(wing; t=t)
      _rhoinf = 1.0
      _qinf = (1.0/2)*_rhoinf*dot(Vinf, Vinf)
    else
      _rhoinf = rhoinf
      _qinf = qinf
    end
    _S = S=="automatic" ? _S = planform_area(wing) : S

    CF, CS, CD, CL = _calculate_force_coeffs(wing, _rhoinf, _qinf, _S; t=t,
                                              lifting_interac=lifting_interac)
    _addsolution(wing, "CFtot", CF)
    _addsolution(wing, "CS", CS)
    _addsolution(wing, "CD", CD)
    _addsolution(wing, "CL", CL)

  ######## AERODYNAMIC FORCE COEFFICIENT PER UNIT SPAN NORMALIZED ##########
  elseif field_name in ["Cftot/CFtot", "Cd/CD", "Cs/CS", "Cl/CL"]
    if rhoinf==nothing && qinf!="automatic"
      println("$(field_name) requested with a given qinf, but rhoinf is missing"
                *". Given qinf will be ignored.")
      Vinf = _aveVinf(wing; t=t)
      _rhoinf = 1.0
      _qinf = (1.0/2)*_rhoinf*dot(Vinf, Vinf)
    elseif qinf=="automatic"
      Vinf = _aveVinf(wing; t=t)
      _rhoinf = 1.0
      _qinf = (1.0/2)*_rhoinf*dot(Vinf, Vinf)
    else
      _rhoinf = rhoinf
      _qinf = qinf
    end
    _S = S=="automatic" ? _S = planform_area(wing) : S

    Cf, Cs, Cd, Cl = _calculate_force_coeffs(wing, _rhoinf, _qinf, _S;
                                              per_unit_span=true, t=t,
                                              lifting_interac=lifting_interac)

    # Converts them into scalars
    s_Cf, s_Cs, s_Cd, s_Cl = [], [], [],[]
    scalars = [s_Cf, s_Cs, s_Cd, s_Cl]
    aveVinf = _aveVinf(wing; t=t)
    for (i, f) in enumerate([Cf, Cs, Cd, Cl])
      for elem in f
        if i==3     # Keeps the sign for the drag direction
          sgn = sign(dot(elem, aveVinf))
        else
          sgn = 1
        end
        push!(scalars[i], sgn*norm(elem))
      end
    end

    # Calculates overalls
    info = fields_summary(wing)

    # Determines the span of the wing
    _min, _max = _span_eff(wing)
    span = _max-_min

    _addsolution(wing, "Cftot/CFtot", s_Cf/(info["CFtot"]/span))
    _addsolution(wing, "Cs/CS", s_Cs/(info["CS"]/span))
    _addsolution(wing, "Cd/CD", s_Cd/(info["CD"]/span))
    _addsolution(wing, "Cl/CL", s_Cl/(info["CL"]/span))

  ######## PANEL AREA ###########################
  elseif field_name=="A"
    _addsolution(wing, field_name, _calculate_areas(wing))


  ######## MOMENTS ###########################
  elseif field_name in ["Mtot", "M_L", "M_M", "M_N"]

    # Center of gravity
    _r_cg = r_cg=="automatic" ? get_CG(wing) : r_cg

    Mtot, M_L, M_M, M_N = _calculate_moments(wing, _r_cg)
    _addsolution(wing, "Mtot", Mtot)
    _addsolution(wing, "M_L", M_L)
    _addsolution(wing, "M_M", M_M)
    _addsolution(wing, "M_N", M_N)


  ######## MOMENTS COEFFICIENTS ###########################
  elseif field_name in ["CMtot", "CM_L", "CM_M", "CM_N"]
    if qinf in ["automatic", nothing, 0.0]
      error("$(field_name) requested but qinf is missing")
    end

    _l = l=="automatic" ? get_barc(wing) : l
    _S = S=="automatic" ? planform_area(wing) : S

    _addsolution(wing, "CMtot", wing.sol["Mtot"]/(qinf*_S*_l))
    _addsolution(wing, "CM_L", wing.sol["M_L"]/(qinf*_S*_l))
    _addsolution(wing, "CM_M", wing.sol["M_M"]/(qinf*_S*_l))
    _addsolution(wing, "CM_N", wing.sol["M_N"]/(qinf*_S*_l))

  ######## ERROR CASE ###################################
  else
    error("Calculation of $(field_name) has not been implemented yet!")
  end
end




"Aerodynamic force calculated through Kutta-Joukowski theorem.
Give it `per_unit_span=true` to calculate the force per unit length of span.
If `lifting_interac`==false, the force doesn't include the induced velocity,
hence, no interaction between lifting surfaces is calculated and induced drag
will be zero, but it is an efficient way of calculating the lift of an
individual lifting surface"
function _calculate_forces(wing, rhoinf::FWrap;
                          t::FWrap=0.0, per_unit_span::Bool=false,
                          lifting_interac::Bool=true)

  F = []
  m = get_m(wing)

  # -------------- CALCULATES TOTAL FORCE F
  # Iterates over horseshoes
  for i in 1:m
    Ap, A, B, Bp, CP, infDA, infDB, Gamma = getHorseshoe(wing, i)
    force = zeros(3)

    # Iterates over bound vortices of the horseshoe
    for BV in [[A,B], [Ap,A], [B,Bp]]
      # Midway distance
      X = (BV[1] + BV[2])/2
      # Freestream velocity (undisturbed+induced)
      V = wing.Vinf(X,t)
      # V = wing.sol["Vinf"][i]
      # NOTE: This line is the bottle neck, here the problem becomes order
      #       O((4*N)^2), with N the number of lattices. It is necessary
      #       to calculate the induced velocity if interaction between
      #       lifting surfaces is wanted. If not, induced drag can be
      #       efficiently calculated in the far field (Trefftz plane)
      if lifting_interac
        V += Vind(wing, X; t=t, ign_col=true)
      end
      # Vinf x (B-A)
      crss = cross(V, BV[2]-BV[1])

      if per_unit_span
        l = norm(cross(V/norm(V),BV[2]-BV[1]))
      else
        l = 1
      end

      # Force
      force += rhoinf * Gamma * crss / l
    end

    push!(F, force)
  end

  # Case of calculating drag at the far field
  if !lifting_interac
    aveVinf = _aveVinf(wing; t=t)
    Ftrefftz = calculate_force_trefftz(wing, aveVinf, rhoinf;
                                          per_unit_span=per_unit_span)
    for i in 1:m
      F[i] += Ftrefftz[i]
    end
  end

  # -------------- DECOMPOSES F INTO COMPONENTS
  S,D,L = _decompose(wing, F; t=t)

  return F, S, D, L
end


function _calculate_force_coeffs(wing, rhoinf::FWrap, qinf::FWrap,
                                  S::FWrap; t::FWrap=0.0,
                                  per_unit_span::Bool=false,
                                  lifting_interac::Bool=true)
  F,SS,D,L = _calculate_forces(wing, rhoinf; t=t, per_unit_span=per_unit_span,
                                    lifting_interac=lifting_interac)
  aux = qinf*S

  return F/aux, SS/aux, D/aux, L/aux
end

"Calculates total moment and decomposition into roll, pitch, and yaw"
function _calculate_moments(wing, r_cg)
  m = get_m(wing)

  # Direction of rolling, pitching, and yawing moments
  unit_L = -wing.Oaxis[1, :]
  unit_M = wing.Oaxis[2, :]
  unit_N = -wing.Oaxis[3, :]
  units = [unit_L, unit_M, unit_N]

  Mtot = []
  M_L = []
  M_M = []
  M_N = []
  for i in 1:m
      # Force on this panel
      Ftot_i = wing.sol["Ftot"][i]
      # Center of this panel
      r_i = get_r(wing, i)
      # Moment on this panel
      Mtot_i = cross( (r_i-r_cg)  , Ftot_i )
      push!(Mtot, Mtot_i)
      # Decomposes the moment
      for (j, M) in enumerate([M_L, M_M, M_N])
        this = dot(Mtot_i, units[j])*units[j]
        push!(M, this)
      end
  end

  return Mtot, M_L, M_M, M_N
end

"Calculates the force induced by the wake sheet assuming it planar in the
direction of Vinf. This method works well only for an isolated lifting surface
and should not be used on interacting lifting surfaces."
function calculate_force_trefftz(wing, Vinf::FArrWrap, rho::FWrap;
                                        per_unit_span::Bool=false)
  HSs = getHorseshoes(wing)
  m = get_m(wing)

  DVinf = Vinf/norm(Vinf)

  Vinds = zeros(FWrap, m)   # Vortex sheet induced velocity at each bound vortex
  for i in 1:m  # Iterates over semi-infinite vortices calculating Vind

    if i==1 || HSs[i-1][3]!=HSs[i][2]
      gamma = HSs[i][8]
    else
      gamma = HSs[i][8] - HSs[i-1][8]
    end
    headVortex = HSs[i][2]      # Head of semi-infinite vortex at A

    bgamma = nothing
    if i==m || HSs[i][3]!=HSs[i+1][2]
      bgamma = -HSs[i][8]
      bheadVortex = HSs[i][3]
    end

    for j in 1:m    # Iterates over HS

      # Determines distance from semi-infinite vortex
      meanAB = (HSs[j][3]+HSs[j][2])/2 # Position half-way A and B of this HS
      H = meanAB - headVortex          # Vector between head vortex and AB mean
      count_proj = H - DVinf*(dot(H, DVinf))  # Counter projection to inf vortex
      h = norm(count_proj)            # Normal distance to semi-infinite vortex

      # Determines sign of induced velocity
      AB = HSs[j][3]-HSs[j][2]        # Direction from A to B
      sgn = sign(dot(cross(DVinf, AB), cross(DVinf, H)))

      if h>1e-8
        Vinds[j] += sgn*gamma/(4*pi*h)
      end

      if bgamma!=nothing # Case this was the end of a wing
        H = meanAB - bheadVortex
        count_proj = H - DVinf*(dot(H, DVinf))
        h = norm(count_proj)
        sgn = sign(dot(cross(DVinf, AB), cross(DVinf, H)))
        if h>1e-8
          Vinds[j] += sgn*bgamma/(4*pi*h)
        end
      end

    end
  end


  F = FArrWrap[zeros(3) for i in 1:m]
  for i in 1:m  # Iterates over bound vortex calculating induced force
    _, A, B, _, _, _, _, Gamma  = HSs[i]
    AB = B-A
    count_proj = AB - DVinf*(dot(AB, DVinf))  # Counter projection about inf vortex
    l = norm(count_proj)        # Length normal to inf vortex direction
    if per_unit_span
      unitspan = l
    else
      unitspan = 1
    end
    F[i] += rho*Gamma*Vinds[i]*l*DVinf / unitspan
  end

  return F
end

"Returns the average Vinf from all control points"
function _aveVinf(wing; t::FWrap=0.0)
  Vinf = zeros(3)
  for i in 1:get_m(wing)
    Vinf += wing.Vinf(getControlPoint(wing, i), t)
  end
  Vinf = Vinf/get_m(wing)
  return Vinf
end

"Decomposes a force field into sideslip, drag, and lift components"
function _decompose(wing, F; t::FWrap=0.0)
  Vinf = _aveVinf(wing; t=t)     # Vinf used for decomposing forces
  # Vinf = mean(wing.sol["Vinf"])

  # Unit vectors
  s_hat = wing.Oaxis[2, :]                                  # Sideslip
  d_hat = Vinf - dot(Vinf, s_hat)*s_hat; d_hat=d_hat/norm(d_hat); # Drag
  l_hat = cross(d_hat, s_hat); l_hat=l_hat/norm(l_hat);     # Lift
  decs_hat = [s_hat, d_hat, l_hat]

  # Checks successful decomposition
  if false==check_coord_sys(decs_hat, raise_error=false)
    error("Failed to decomposed forces."*
          " Possible cause is that drag is aligned with sideslip"*
          "\n\ts_hat=$s_hat\td_hat=$d_hat\tl_hat=$l_hat")
  end

  # Decomposes
  S, D, L = [], [], []
  for (i,dec) in enumerate([S,D,L])
    dec_hat = decs_hat[i]
    for this_force in F
      push!(dec, dot(this_force, dec_hat)*dec_hat)
    end
  end
  return S,D,L

end

"Calculates the area of each panel"
function _calculate_areas(wing)
  A = [get_A(wing, i) for i in 1:get_m(wing)]
  return A
end

"Returns a dictionary summarizing key parameters on a Wing or WingSystem"
function fields_summary(wing; drag_trick=false)
  dict = Dict()

  # Iterates over each field
  for field_name in keys(wing.sol)
    if drag_trick && field_name in ["D", "CD"]
      val = 0.0
    elseif FIELDS[field_name][2]=="vector"
      val = [0.0,0.0,0.0]
    else
      val = 0.0
    end

    # Adds each panel's value together
    for elem in wing.sol[field_name]
      val += (drag_trick && field_name in ["D", "CD"]) ? norm(elem) : elem
    end

    # Saves it
    if FIELDS[field_name][2]=="vector"
      dict[field_name] = norm(val)
    else
      dict[field_name] = val
    end

  end
  # Special information
  dict["control_points"] = get_m(wing)
  return dict
end

"Returns the span effective of a Wing or WingSystem"
function _span_eff(wing; min=Inf, max=-Inf)
  if typeof(wing)==Wing
    this_max = maximum(wing._ywingdcr)
    this_min = minimum(wing._ywingdcr)
    _max = max>this_max ? max : this_max
    _min = min<this_min ? min : this_min
    return _min, _max
  elseif typeof(wing) in [typeof(Wing[]), Array{Any,1}]
    _min, _max = min, max
    for w in wing
      _min, _max = _span_eff(w; min=_min, max=_max)
    end
    return _min, _max
  elseif typeof(wing)==WingSystem
    return _span_eff(wing.wings)
  else
    error("Logic error!")
  end
end

##### END OF WING POSTPROCESSING ###############################################
