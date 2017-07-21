



################################################################################
# CALCULATION WING PROPERTIES
################################################################################
"""
Returns the area of the m-th panel, with `self` either Wing or WingSystem.
For derivation, see notebook entry 20170619.
"""
function get_A(self, m::Int64)
  Ap, A, B, Bp, CP, infDA, infDB, Gamma = getHorseshoe(self, m, nothing)
  CPAB = (A+B)/2
  lcp = norm(CPAB-CP)/(pm-pn)
  h = sqrt( (A[2]-B[2])^2 + (A[3]-B[3])^2  )
  Apanel = lcp*h
  return Apanel
end

"""
Returns the center of gravity of the m-th panel, with `self` either Wing or
WingSystem.
"""
function get_r(self, m::Int64)
  Ap, A, B, Bp, CP, infDA, infDB, Gamma = getHorseshoe(self, m, nothing)
  CPAB = (A+B)/2
  r = (CPAB+CP)/2
  return r
end

"""
Returns the center of gravity of the wing, with `self` either Wing or
WingSystem.
"""
function get_CG(self)
  sum_A = 0.0
  sum_rA = zeros(3)
  # Iterates over each panel
  for i in 1:get_m(self)
    r = get_r(self, i)
    A = get_A(self, i)
    sum_A += A
    sum_rA += r*A
  end
  r_cg = sum_rA/sum_A
  return r_cg
end

"""
  Returns the mean chord length of a Wing or WingSystem
"""
function get_barc(self)
  sum_lA = 0.0
  sum_A = 0.0

  # Iterates over each panel
  for i in 1:get_m(self)
    Ap, A, B, Bp, CP, infDA, infDB, Gamma = getHorseshoe(self, i)
    CPAB = (A+B)/2

    # Mean chord of this panel
    lcp = norm(CPAB-CP)/(pm-pn)
    # Area of this panel
    Apanel = get_A(self, i)

    sum_A += Apanel
    sum_lA += lcp*Apanel
  end
  barc = sum_lA/sum_A
  return barc
end



################################################################################
# UTILITIES
################################################################################
"""
    `simpleWing(b, ar, tr, alpha, lambda, gamma; alpha_tip=alpha, n=20, r=2.0)`
  Generates a single-section wing.

  # Arguments
  *   b       :   (float) span.
  *   ar      :   (float) aspect ratio defined as the span over the tip's chord.
  *   tr      :   (float) taper ratio (tip chord / root chord).
  *   alpha   :   (float) angle of attack in degrees.
  *   lambda  :   (float) sweep in degrees.
  *   gamma   :   (float) dihedral in degrees.
  *   (OPTIONALS)
  *   alpha_tip : (float) angle of attack at the tip. Use this to generate twist
  *   n       :   (int)   number of horseshoes per side of the wing.
  *   r       :   (float) horseshoes' expansion ratio
"""
function simpleWing(b::Float64, ar::Float64, tr::Float64,
                    alpha::Float64, lambda::Float64, gamma::Float64;
                    alpha_tip=nothing,
                    n::Int64=20, r::Float64=2.0)
  cr = 1/tr
  c_tip = b/ar
  c_root = cr*c_tip
  alpha_t = alpha_tip==nothing ? alpha : alpha_tip

  y_tip = b/2
  x_tip = y_tip*tan(lambda*pi/180)
  z_tip = y_tip*tan(gamma*pi/180)

  wing = Wing(x_tip, -y_tip, z_tip, c_tip, alpha_t)
  addchord(wing, 0.0, 0.0, 0.0, c_root, alpha, n, r)
  addchord(wing, x_tip, y_tip, z_tip, c_tip, alpha_t, n, 1/r)

  return wing
end

"""
 Saves the wing domain in VTK legacy format.

 Generates the following files:
 * `[filename]_dom.vtk`
"""
function savedomain(self::Wing, filename::String;
                  save_horseshoes::Bool=true,
                  path::String="", comment::String="",
                  num=nothing)
  aux = num!=nothing ? ".$num" : ""
  ext = "_vlm"*aux*".vtk"

  if path !=""
    _path = string(path, (path[end]!="/" ? "/" : ""))
  else
    _path = ""
  end

  # HEADER
  n = get_m(self)+1 # Number of chords
  f = open(string(_path, filename, ext), "w")
  header = "# vtk DataFile Version 4.0" # File version and identifier
  header = string(header, "\n", "[MyVLM domain - lattices] ", comment) # Title
  header = string(header, "\n", "ASCII") # File format
  header = string(header, "\n", "DATASET UNSTRUCTURED_GRID")
  write(f, header)

  # POINTS
  nle = n
  nte = n
  ncp = n-1
  if save_horseshoes
    nhs = n-1
  else
    nhs = 0
  end
  write(f, string("\n", "POINTS ", nle+nte+ncp+nhs*6, " float"))
  ## Leading edge
  for i in 1:nle
    LE = getLE(self, i)
    line1 = string(LE[1], " ", LE[2], " ", LE[3])
    write(f, string("\n", line1))
  end
  ## Trailing edge
  for i in 1:nte
    TE = getTE(self, i)
    line2 = string(TE[1], " ", TE[2], " ", TE[3])
    write(f, string("\n", line2))
  end
  ## Control points (m)
  for i in 1:ncp
    CP = getControlPoint(self, i)
    line = string(CP[1], " ", CP[2], " ", CP[3])
    write(f, string("\n", line))
  end
  ## Horseshoes
  x_vor_end = maximum(self._xtwingdcr)*1.25
  for i in 1:nhs
    Ap, A, B, Bp, CP, infDA, infDB, Gamma  = getHorseshoe(self, i)

    # Calculates how long to extend the semi-infinite vortices
    # factor = (x_vor_end - Ap[1])/infDA[1]
    factor = (x_vor_end - self._xtwingdcr[i])/infDA[1]
    Apinf = Ap + factor*infDA
    # factor = (x_vor_end - Bp[1])/infDB[1]
    factor = (x_vor_end - self._xtwingdcr[i+1])/infDB[1]
    Bpinf = Bp + factor*infDB

    for point in [Apinf, Ap, A, B, Bp, Bpinf]
      line = string(point[1], " ", point[2], " ", point[3])
      write(f, string("\n", line))
    end
  end

  # CELLS
  write(f, string("\n\n", "CELLS ", n-1+ncp+nhs,
                  " ", (n-1)*5 + ncp*2 + nhs*7))
  ## Lattices
  for i in 0:n-2
    line = string(4, " ", i, " ", i+nle, " ", i+nle+1, " ", i+1)
    write(f, string("\n", line))
  end
  ## Control points
  for i in nle+nte:nle+nte+ncp-1
    line = string(1, " ", i)
    write(f, string("\n", line))
  end
  ## Horseshoes
  aux1 = nle+nte+ncp+nhs*6
  count = 1
  line = string(6, " ")
  for i in nle+nte+ncp:aux1-1
    line = string(line, " ", i)
    if count%6==0
      write(f, string("\n", line))
      line = string(6, " ")
    end
    count += 1
  end


  # CELL TYPES
  write(f, string("\n\n", "CELL_TYPES ", n-1+ncp+nhs))
  ## Lattices
  for i in 0:n-2
    write(f, string("\n", 9))
  end
  ## Control points
  for i in 1:ncp
    write(f, string("\n", 1))
  end
  ## Horseshoes
  for i in 1:nhs
    write(f, string("\n", 4))
  end


  # FIELDS
  initiated = false
  for field_name in keys(self.sol)
    if false==(field_name in keys(FIELDS))
      error(string("CRITICAL ERROR: field ", field_name, " not found in FIELDS"))
    end

    if initiated==false
      write(f, string("\n\n", "CELL_DATA ", n-1+ncp+nhs))
      initiated = true
    end

    if FIELDS[field_name][2]=="vector"
      write(f, string("\n\n", "VECTORS ", field_name," float"))
      for i in 1:n-1+ncp+nhs
            vect = self.sol[field_name][(i-1)%ncp+1]
            line = string(vect[1], " ", vect[2], " ", vect[3])
            write(f, string("\n", line))
      end
    elseif FIELDS[field_name][2]=="scalar"
      write(f, string("\n\n", "SCALARS ", field_name," float"))
      write(f, string("\n", "LOOKUP_TABLE default"))
      for i in 1:n-1+ncp+nhs
            if field_name=="Gamma"
              sclr = self.sol[field_name][(i-1)%ncp+1][6]
            else
              sclr = self.sol[field_name][(i-1)%ncp+1]
            end
            try
              line = string(isnan(sclr) ? -1 : sclr)
            catch
              error(string(sclr))
            end

            write(f, string("\n", line))
      end
    else
      error(string("CRITICAL ERROR: field type ", FIELDS[field_name][2],
            " hasn't been implemented!"))
    end
  end

  close(f)
end




################################################################################
# CALCULATIONS
################################################################################
"""
  `V(self::Wing or WingSystem, P::Float64[])`
Calculates the induced velocity at point P (Gamma field must have been solved)
"""
function V(self, P; ign_col::Bool=false)
  if false==("Gamma" in keys(self.sol))
    error("ERROR: Gamma field not found!")
  end

  m = get_m(self)

  V_tot = zeros(3)
  # Iterates over every horseshoe in the wing
  for i in 1:m
    this_HS = getHorseshoe(self, i, nothing)
    V_this = VLMSolver.V(this_HS, P; ign_col=ign_col)
    V_tot += V_this
  end
  return V_tot
end



################################################################################
# ALGEBRA
################################################################################
"""
  Rotates and translates the vector V.

  Receives the i', j', k' unit vectors of an euclidean system with origin T, and
  returns V'.
"""
function transform(V::typeof(Float64[]),
                    ip::typeof(Float64[]), jp::typeof(Float64[]),
                    kp::typeof(Float64[]), T::typeof(Float64[]))
  M = zeros(3,3)
  units = [ip, jp, kp]
  for i in 1:3
    for j in 1:3
      M[i,j] = units[i][j]
    end
  end
  return M*(V-T)
end

function transform(Vs::Array{Array{Float64,1},1},
                    ip::typeof(Float64[]), jp::typeof(Float64[]),
                    kp::typeof(Float64[]), T::typeof(Float64[]))
  out = Array{Float64,1}[]
  for V in Vs
    push!(out, transform(V, ip, jp, kp, T))
  end
  return out
end

"""
  Rotates and translates back a vector V' that had been rotated and translated
  into the system (i', j', k') with origin T, and
  returns the original V.
"""
function countertransform(Vp::typeof(Float64[]),
                            ip::typeof(Float64[]), jp::typeof(Float64[]),
                            kp::typeof(Float64[]), T::typeof(Float64[]))
  M = zeros(3,3)
  units = [ip, jp, kp]
  for i in 1:3
    for j in 1:3
      M[i,j] = units[i][j]
    end
  end
  return inv(M)*Vp + T
end

function countertransform(Vps::Array{Array{Float64,1},1},
                    ip::typeof(Float64[]), jp::typeof(Float64[]),
                    kp::typeof(Float64[]), T::typeof(Float64[]))
  out = Array{Float64,1}[]
  for Vp in Vps
    push!(out, countertransform(Vp, ip, jp, kp, T))
  end
  return out
end
