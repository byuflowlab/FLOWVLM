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
                        S="automatic", l="automatic", r_cg="automatic")

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

    F, SS, D, L = _calculate_forces(wing, rhoinf)
    _addsolution(wing, "Ftot", F)
    _addsolution(wing, "S", SS)
    _addsolution(wing, "D", D)
    _addsolution(wing, "L", L)

  ######## FORCE COEFFICIENTS ###########################
  elseif field_name in ["CFtot", "CL", "CS", "CD"]
    if rhoinf==nothing && qinf!="automatic"
      println("$(field_name) requested with a given qinf, but rhoinf is missing"
                *". Given qinf will be ignored.")
      Vinf = _aveVinf(wing)
      _rhoinf = 1.0
      _qinf = (1.0/2)*_rhoinf*dot(Vinf, Vinf)
    elseif qinf=="automatic"
      Vinf = _aveVinf(wing)
      _rhoinf = 1.0
      _qinf = (1.0/2)*_rhoinf*dot(Vinf, Vinf)
    else
      _rhoinf = rhoinf
      _qinf = qinf
    end
    _S = S=="automatic" ? _S = planform_area(wing) : S

    CF, CS, CD, CL = _calculate_force_coeffs(wing, _rhoinf, _qinf, S)
    _addsolution(wing, "CFtot", CF)
    _addsolution(wing, "CS", CS)
    _addsolution(wing, "CD", CD)
    _addsolution(wing, "CL", CL)

  ######## AERODYNAMIC FORCE COEFFICIENT PER UNIT SPAN NORMALIZED ##########
  elseif field_name in ["Cftot/CFtot", "Cd/CD", "Cs/CS", "Cl/CL"]
    if rhoinf==nothing && qinf!="automatic"
      println("$(field_name) requested with a given qinf, but rhoinf is missing"
                *". Given qinf will be ignored.")
      Vinf = _aveVinf(wing)
      _rhoinf = 1.0
      _qinf = (1.0/2)*_rhoinf*dot(Vinf, Vinf)
    elseif qinf=="automatic"
      Vinf = _aveVinf(wing)
      _rhoinf = 1.0
      _qinf = (1.0/2)*_rhoinf*dot(Vinf, Vinf)
    else
      _rhoinf = rhoinf
      _qinf = qinf
    end
    _S = S=="automatic" ? _S = planform_area(wing) : S

    Cf, Cs, Cd, Cl = _calculate_force_coeffs(wing, _rhoinf, _qinf, S;
                                              per_unit_span=true)

    # Converts them into scalars
    s_Cf, s_Cs, s_Cd, s_Cl = [], [], [],[]
    scalars = [s_Cf, s_Cs, s_Cd, s_Cl]
    for (i, f) in enumerate([Cf, Cs, Cd, Cl])
      for elem in f
        push!(scalars[i], norm(elem))
      end
    end

    # Calculates overalls
    info = fields_summary(wing)

    # Determines the span of the wing
    span = maximum(wing._ywingdcr)-minimum(wing._ywingdcr)

    _addsolution(wing, "Cftot/CFtot", s_Cf/(info["CFtot"]/span))
    _addsolution(wing, "Cs/CS", s_Cs/(info["CS"]/span))
    _addsolution(wing, "Cd/CD", s_Cd/(info["CD"]/span))
    _addsolution(wing, "Cl/CL", s_Cl/(info["CL"]/span))

  ######## ERROR CASE ###################################
  else
    error("Calculation of $(field_name) has not been implemented yet!")
  end
end




"Aerodynamic force calculated through Kutta-Joukowski theorem.
Give it `per_unit_span=true` to calculate the force per unit length of span."
function _calculate_forces(wing, rhoinf::Float64;
                          t::Float64=0.0, per_unit_span::Bool=false)

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
      V = wing.Vinf(X,t) + Vind(wing, X; t=t, ign_col=true)
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

  # -------------- DECOMPOSES F INTO COMPONENTS
  S,D,L = _decompose(wing, F)

  return F, S, D, L
end


function _calculate_force_coeffs(wing, rhoinf::Float64, qinf::Float64,
                                  S::Float64; t::Float64=0.0,
                                  per_unit_span::Bool=false)
  F,SS,D,L = _calculate_forces(wing, rhoinf; t=t, per_unit_span=per_unit_span)
  aux = qinf*S

  return F/aux, SS/aux, D/aux, L/aux
end

"Returns the average Vinf from all control points"
function _aveVinf(wing; t::Float64=0.0)
  Vinf = zeros(3)
  for i in 1:get_m(wing)
    Vinf += wing.Vinf(getControlPoint(wing, i), t)
  end
  Vinf = Vinf/get_m(wing)
  return Vinf
end

"Decomposes a force field into sideslip, drag, and lift components"
function _decompose(wing, F)
  Vinf = _aveVinf(wing)     # Vinf used for decomposing forces

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

"Returns a dictionary summarizing key parameters on a Wing or WingSystem"
function fields_summary(wing)
  dict = Dict()

  # Iterates over each field
  for field_name in keys(wing.sol)
    if FIELDS[field_name]=="vector"
      val = [0.0,0.0,0.0]
    else
      val = 0.0
    end

    # Adds each panel's value together
    for elem in wing.sol[field_name]
      val += elem
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

##### END OF WING POSTPROCESSING ###############################################
