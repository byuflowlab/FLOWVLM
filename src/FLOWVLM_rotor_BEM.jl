#=
  TODO
  - [ ] Code your own 3D correction since airfoilprepy doesn't agree with its
        validation case. Include the 3D corrections recommended in McCrink, M.
        H., & Gregory, J. W. (2017), *Blade Element Momentum Modeling of
        Low-Reynolds Electric Propulsion Systems*.
  - [ ] Add Glauert's Mach correction.
=#

import Roots



"""
  Returns the inflow angle `phi` (in RADIANS) of the blade element.

  **Arguments**
  * `cl::Real`        : Lift coefficient function of AOA in DEGREES.
  * `cd::Real`        : Drag coefficient function of AOA in DEGREES.
  * `roR::Real`       : Radial position r/R.
  * `sigma::Real`     : Annulus' solidity B*c/(2*pi*r).
  * `lambda_r::Real`  : Annulus' tip-speed ratio omega*r/Vinf.
  * `B::Signed`       : Number of blades.
"""
function calc_phi(cl, cd, roR::Real, sigma::Real, lambda_r::Real, B::Signed)

  # Residual function
  function R(phi)
    a, ap = calc_indfactors(phi, cl(180/pi*phi), cd(180/pi*phi), roR, sigma, B)

    return sin(phi)/(1+a) - cos(phi)/(lambda_r*(1-ap))
  end

  # NOTE: Defining the range as [0, pi/2] hard codes the solution to exclude
  #       reversed flow.
  phi = Roots.fzero(R, [0, pi/2])
  a, ap = calc_indfactors(phi, cl(180/pi*phi), cd(180/pi*phi), roR, sigma, B)

  return phi, a, ap
end

"""
  Returns both axial and tangential induction factors `(a, ap)` at a blade
element of radial position `roR` with chord/position ratio `cor` (chord length c
over the local radial location r), overall tip-speed ratio `lambda` (velocity of
the rotor's tip over freestream omega*R/Vinf), equivalent local tip-speed ratio
`lambda_r` (omega*r/Vinf), an assumed local angle `phi` (rad), and lift and
drag coefficients `cl`, `cd` of the element at `phi`.

  **Arguments**
  * `phi::Real`       : Assumed induced angle in RADIANS.
  * `cl::Real`        : Lift coefficient at `phi`.
  * `cd::Real`        : Drag coefficient at `phi`.
  * `roR::Real`       : Radial position r/R.
  * `sigma::Real`     : Annulus' solidity B*c/(2*pi*r).
  * `B::Signed`       : Number of blades.
"""
function calc_indfactors(phi::Real, cl::Real, cd::Real, roR::Real, sigma::Real,
                                                                      B::Signed)
  cphi = cos(phi)
  sphi = sin(phi)

  # Normal and tangential force coefficients
  cn = cl*cphi - cd*sphi
  ct = cl*sphi + cd*cphi

  # Tip loss correction factor
  F = Prandtl_F(B, phi, roR)

  # Induction factors
  a = 1/( 4*F*sphi^2/(sigma*cn) - 1 )
  ap = 1/( 4*F*sphi*cphi/(sigma*ct) + 1 )

  return a, ap
end


"Prandtl's tip loss correction factor, with `B` number of blades,
`phi` local inflow angle (rad), and `roR` radial position."
function Prandtl_F(B, phi, roR)
    f = B/(2*abs(sin(phi))) * (1/roR - 1)
    F = 2/pi * acos(exp(-f))
    return F
end
