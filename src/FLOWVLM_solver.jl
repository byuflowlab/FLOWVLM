
module VLMSolver

# Criteria of colinearity
const col_crit = 1/10^12  # NOTE: Anything less than 1/10^15 reaches float precision.
global n_col = 0          # Number of colinears found
# global colinears = []     # Store colinear points here
global mute_warning = false
function _mute_warning(booln::Bool)
  global mute_warning = booln
end


# A horseshoe is defined as a 5-segments vortex by the array
# HS = [Ap, A, B, Bp, CP, infDA, infDB, Gamma], with
#
#   * `Ap::Array{Float64,1}`    : A-side trailing edge.
#   * `A::Array{Float64,1}`     : A-side of the bound vortex.
#   * `B::Array{Float64,1}`     : B-side of the bound vortex.
#   * `Bp::Array{Float64,1}`    : B-side trailing edge.
#   * `CP::Array{Float64,1}`    : Control point of the lattice associated to the HS.
#   * `infDA::Array{Float64,1}` : Direction of A-side semi-infinite vortex.
#   * `infDB::Array{Float64,1}` : Direction of B-side semi-infinite vortex.
#   * `Gamma::Float64 or nothing`: Circulation of the horseshoe.
#
# infDA and infDB must be unitary vectors pointing from the trailing edge in
# the direction of infinite.


################################################################################
# SOLVER
################################################################################
"""
  `VLMsolve(HSs, Vinf[, t=0])`
Solves for the circulation of a collection of horseshoe vortices with the
boundary condition of no-through flow.

  # Arguments
  * `HSs::Array{Float64,1}` : Collection of horseshoes.
  * `Vinf::function(X,t)`   : Undisturbed freestream.

  # Optional arguments
  * `t::Float64`            : Time argument for evaluation of Vinf.
"""
function solve(HSs::Array{Array{Any,1},1}, Vinf; t::Float64=0.0)

  n = size(HSs)[1]    # Number of horseshoes
  G = zeros(n, n)     # Geometry matrix
  Vn = zeros(n)       # Normal velocity matrix

  # ------------ BUILD MATRICES G AND Vn ------------
  # Iterates over control points
  for i in 1:n
    Api, Ai, Bi, Bpi, CPi, infDAi, infDBi, Gamma = HSs[i]

    # Calculates the normal of the panel
    crss = cross(CPi-Ai, Bi-Ai)
    nhat = crss/sqrt(dot(crss,crss))

    # Iterates over horseshoes
    for j in 1:n
      HS = HSs[j]
      GeomFac = V(HS, CPi)
      Gij = (1/4/pi)*GeomFac
      Gijn = dot(Gij, nhat)
      G[i,j] = Gijn
    end

    # Normal component of Vinf
    this_Vinf = Vinf(CPi, t)
    Vinfn = dot(this_Vinf, nhat)
    Vn[i] = -Vinfn
  end

  # ------------ SOLVES FOR GAMMA ------------
  invG = inv(G)
  Gamma = invG * Vn
  return Gamma
end


"""
Returns the induced velocity at `C` by horseshoe `HS`.
It returns the geometric factor if `Gamma`==nothing.
"""
function V(HS::Array{Any,1}, C; ign_col::Bool=false)
  Ap, A, B, Bp, CP, infDA, infDB, Gamma = HS
  VApinf = _V_Ainf_in(Ap, infDA, C, Gamma; ign_col=ign_col)
  VApA = _V_AB(Ap, A, C, Gamma; ign_col=ign_col)
  VAB = _V_AB(A, B, C, Gamma; ign_col=ign_col)
  VBBp = _V_AB(B, Bp, C, Gamma; ign_col=ign_col)
  VBpinf = _V_Ainf_out(Bp, infDB, C, Gamma; ign_col=ign_col)
  V = VApinf + VApA + VAB + VBBp + VBpinf
  return V
end
##### END OF SOLVER ############################################################

################################################################################
# CALCULATION
################################################################################
"""
Returns the induced velocity of the bound vortex AB on point `C`.
Give gamma=nothing to return the geometric factor (Fac1*Fac2).
"""
function _V_AB(A::Array{Float64,1}, B, C, gamma; ign_col::Bool=false)

  r0 = B-A
  r1 = C-A
  r2 = C-B
  crss = cross(r1,r2)
  magsqr = dot(crss, crss)

  # Checks colinearity
  if _check_collinear(magsqr, col_crit; ign_col=ign_col)
    if ign_col==false && n_col==1 && mute_warning==false
      println("\n\t magsqr:$magsqr \n\t A:$A \n\t B:$B \n\t C:$C")
    end
    return [0.0,0.0,0.0]
  end

  F1 = crss/magsqr
  aux = r1/sqrt(dot(r1,r1)) - r2/sqrt(dot(r2,r2))
  F2 = dot(r0, aux)

  if gamma==nothing
    return F1*F2
  else
    return (gamma/4/pi)*F1*F2
  end
end


"""
Returns the induced velocity on point `C` by the semi-infinite vortex `Ainf.`
The vortex outgoes from `A` in the direction of `infD`.
Give `gamma=nothing` to return the geometric factor `(Fac1*Fac2)`.
"""
function _V_Ainf_out(A::Array{Float64,1},
                      infD::Array{Float64,1}, C, gamma;
                      ign_col::Bool=false)
  AC = C-A

  unitinfD = infD / sqrt(dot(infD, infD))
  AAp = dot(unitinfD, AC)*unitinfD
  Ap = AAp + A

  boundAAp = _V_AB(A, Ap, C, gamma; ign_col=ign_col)

  ApC = C - Ap
  crss = cross(infD, ApC)
  mag = sqrt(dot(crss, crss))

  # Checks colinearity
  if _check_collinear(mag, col_crit; ign_col=ign_col)
    if ign_col==false && n_col==1 && mute_warning==false
      println("\n\t magsqr:$magsqr \n\t A:$A \n\t infD:$infD \n\t C:$C")
    end
    return [0.0,0.0,0.0]
  end

  h = mag/sqrt(dot(infD, infD))
  n = crss/mag
  F = n/h

  if gamma==nothing
    return F + boundAAp
  else
    return (gamma/4/pi)*F + boundAAp
  end
end


"""
Returns the induced velocity on point `C` by the semi-infinite vortex `Ainf.`
The vortex incomes from the direction of `infD` to `A`.
Give `gamma=nothing` to return the geometric factor `(Fac1*Fac2)`.
"""
function _V_Ainf_in(A::Array{Float64,1}, infD::Array{Float64,1}, C,
               gamma; ign_col::Bool=false)
  aux = _V_Ainf_out(A, infD, C, gamma; ign_col=ign_col)
  return (-1)*aux
end

"Checks collinearity"
function _check_collinear(magsqr, col_crit; ign_col::Bool=false)
  if magsqr<col_crit
    if ign_col==false
      if n_col==0 && mute_warning==false
        warn("Requested induced velocity on a point colinear with vortex."*
                " Returning 0")
      end
      global n_col += 1
      # push!(colinears, C)
      if n_col%10==0 && mute_warning==false
        warn("Number of colinears found: $n_col")
      end
    end
    return true
  end
  return false
end
##### END OF CALCULATIONS ######################################################

end # END OF MODULE
