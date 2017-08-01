"Useful tools for postprocessing potential flows"
module PP


################################################################################
# FLUID DOMAIN CLASS
################################################################################
"""
    `FluidDomain(P_min, P_max, NDIVS)`

Generates a meshed fluid domain to calculate fields.

  # Arguments
  *   `P_min`       : (Float64[3]) minimum point of the fluid domain.
  *   `P_max`       : (Float64[3]) maximum point of the fluid domain.
  *   `NDIVS`       : (Float64[3]) Number of divisions on the fluid domain on
                      each coordinate.
"""
type FluidDomain
  # User inputs
  P_min::Array{Float64,1}  # minimum point of the fluid domain
  P_max::Array{Float64,1}  # maximum point of the fluid domain
  NDIVS::Array{Int64,1}    # Number of divisions of the fluid domain in
                            # each coordinate.

  # Properties
  D::typeof(zeros(0,0,0,0))     # Fluid domain in the global reference frame
  sol::typeof(Dict())           # Solutions fields
  field_types::Dict{Any,Any}  # Type of each field
  O::Array{Float64,1}           # Origin of local reference frame
  Oaxis::Array{Float64,2}       # Unit vectors of the local reference frame
  invOaxis::Array{Float64,2}    # Inverse unit vectors


  # Data storage
  _D::typeof(zeros(0,0,0,0))    # Fluid domain in the local reference frame

  FluidDomain(P_min, P_max, NDIVS,
              D=_generatedomain(P_min, P_max, NDIVS),
                sol=Dict(),
                field_types=Dict(),
                O=[0.0,0.0,0.0],
                Oaxis=[1.0 0 0; 0 1 0; 0 0 1],
                invOaxis=[1.0 0 0; 0 1 0; 0 0 1],
              _D = D,
      ) = new(P_min, P_max, NDIVS,
              D, sol, field_types, O, Oaxis, invOaxis,
              _D)
end

"""
  `calculate(fluiddomain, fields)`
Calculates for the requested fields, with fields in the following format:

fields = [ field1, field2, ...]
with fieldi = Dict( "field_name" => field_name
                    "field_type" => "scalar" or "vector"
                    "field_function" => field_function
                   )
and field_function a function that receives a position vector and returns
either a scalar or a vector as indicated.
"""
function calculate(self::FluidDomain, fields::Array{Dict{String,Any},1};
              verbose::Bool=true)
  D = self.D
  NDIVS = self.NDIVS

  for field in fields
    field_type = field["field_type"]
    field_name = field["field_name"]
    field_function = field["field_function"]

    if field_name in keys(self.sol) && verbose
      warn("Overwritting existing field: $field_name")
    end

    # Initiates the field
    N = NDIVS + 1
    if field_type=="vector"
      F = zeros(N[1], N[2], N[3], 3 )
    elseif field_type=="scalar"
      F = zeros(N[1], N[2], N[3])
    else
      error("Unkown field type $field_type")
    end

    # Calculates the field
    for i in 1:N[1]
      for j in 1:N[2]
        for k in 1:N[3]

          X = D[i,j,k, 2:4]       # Position vector
          Val = field_function(X) # Value of the field at this position

          if field_type=="vector"
            F[i,j,k,:] = Val[:]
          else
            F[i,j,k] = Val
          end

        end
      end
    end

    self.sol[field_name] = F
    self.field_types[field_name] = field_type
  end
end


"""
  `setcoordsystem(self, O, Oaxis; check=true)`

Redefines the local coordinate system of the domain, where `O` is the new origin
and `Oaxis` is the matrix [i; j; k] of unit vectors
"""
function setcoordsystem(self::FluidDomain, O::Array{Float64,1},
                            Oaxis::Array{Float64,2};
                            check=true)

  if check; check_coord_sys(Oaxis); end;

  self.O = O
  self.Oaxis = Oaxis
  self.invOaxis = inv(Oaxis)
  _reset(self)
  _calculate_D(self::FluidDomain)
end


function setcoordsystem(self::FluidDomain, O::Array{Float64,1},
                            Oaxis::Array{Array{Float64,1},1};
                            check=true)
  dims = 3
  M = zeros(dims, dims)
  for i in 1:dims
    M[i, :] = Oaxis[i]
  end
  setcoordsystem(self, O, M; check=check)
end

"""
 Saves the fluid domain in VTK legacy format.

 Generates the following files:
 * `[filename]_fdom.vtk`
"""
function save(self::FluidDomain, filename::String;
                    num::Int64=-1, time=nothing,
                    path::String="", comment::String="")
  if num!=-1
    nt = ".$(num)"
  else
    nt = ""
  end
  ext = "_fdom$(nt).vtk"

  if path !=""
    _path = string(path, (path[end]!="/" ? "/" : ""))
  else
    _path = ""
  end


  # HEADER
  f = open(string(_path, filename, ext), "w")
  header = "# vtk DataFile Version 4.0" # File version and identifier
  header = string(header, "\n", "[MyVLM fluid domain] ", comment) # Title
  header = string(header, "\n", "ASCII") # File format
  header = string(header, "\n", "DATASET UNSTRUCTURED_GRID")
  write(f, header)

  # TIME
  if time!=nothing
    line0 = "\nFIELD FieldData 1"
    line1 = "\nSIM_TIME 1 1 double"
    line2 = "\n$(time)"
    write(f, line0*line1*line2)
  end

  # POINTS
  NDIVS = self.NDIVS
  N = NDIVS+1
  n = N[1]*N[2]*N[3] # Number of nodes
  write(f, string("\n", "POINTS ", n, " float"))
  for i in 1:N[1]
    for j in 1:N[2]
      for k in 1:N[3]
        P = self.D[i,j,k,2:4]
        line = string(P[1], " ", P[2], " ", P[3])
        write(f, string("\n", line))
      end
    end
  end

  # CELLS
  # write(f, string("\n\n", "CELLS ", n, " ", n*2))
  # for i in 0:n-1
  #   line = string(1, " ", i)
  #   write(f, string("\n", line))
  # end
  write(f, string("\n\n", "CELLS ", NDIVS[1]*NDIVS[2]*NDIVS[3],
                  " ", NDIVS[1]*NDIVS[2]*NDIVS[3]*9))
  m = 1
  for i in 1:N[1]
    if i%N[1]==0
      m += N[3]*N[2]
    else
      for j in 1:N[2]
        if j%N[2]==0
          m += N[3]
        else
          for k in 1:N[3]
            if k%N[3]==0
              nothing
            else
              n0 = m-1
              n1 = n0 + N[3]*N[2]
              n2 = n1 + N[3]
              n3 = n0 + N[3]
              n4 = n0 + 1
              n5 = n1 + 1
              n6 = n2 + 1
              n7 = n3 + 1
              line = string(8, " ", n0, " ", n1, " ", n2, " ", n3, " ", n4, " ", n5,
                                 " ", n6, " ", n7)
              write(f, string("\n", line))
            end
            m += 1
          end
        end
      end
    end
  end

  # CELL TYPES
  # write(f, string("\n\n", "CELL_TYPES ", n))
  write(f, string("\n\n", "CELL_TYPES ", NDIVS[1]*NDIVS[2]*NDIVS[3]))
  ## Control points
  for i in 1:NDIVS[1]*NDIVS[2]*NDIVS[3]
    # write(f, string("\n", 1))
    write(f, string("\n", 12))
  end

  # FIELDS
  initiated = false
  for field_name in keys(self.sol)
    field_type = self.field_types[field_name]

    if initiated==false
      write(f, string("\n\n", "POINT_DATA ", n))
      initiated = true
    end

    if field_type=="vector"
      write(f, string("\n\n", "VECTORS ", field_name," float"))
      for i in 1:N[1]
        for j in 1:N[2]
          for k in 1:N[3]
            vect = self.sol[field_name][i,j,k,:]
            line = string(vect[1], " ", vect[2], " ", vect[3])
            write(f, string("\n", line))
          end
        end
      end
    elseif field_type=="scalar"
      write(f, string("\n\n", "SCALARS ", field_name," float"))
      write(f, string("\n", "LOOKUP_TABLE default"))
      for i in 1:N[1]
        for j in 1:N[2]
          for k in 1:N[3]
            sclr = self.sol[field_name][i,j,k]
            line = string(sclr)
            write(f, string("\n", line))
          end
        end
      end
    else
      error(string("CRITICAL ERROR: field type $field_type"*
                    " hasn't been implemented yet!"))
    end
  end

  close(f)
end


##### INTERNAL FUNCTIONS #######################################################
"""
  Generates the nodes of the domain.
"""
function _generatedomain(P_min::Array{Float64,1}, P_max::Array{Float64,1},
                          NDIVS::Array{Int64,1})
  for i in 1:3
    if P_min[i]>P_max[i]
      error("Invalid P_min and P_max!")
    end
  end

  N = NDIVS + 1 # Number of nodes on each direction
  D = zeros(N[1], N[2], N[3], 4)
  xlen = P_max[1]-P_min[1]
  ylen = P_max[2]-P_min[2]
  zlen = P_max[3]-P_min[3]
  dx = xlen/N[1]
  dy = ylen/N[2]
  dz = zlen/N[3]
  x = P_min[1]
  n = 1
  for i in 1:N[1]
    y = P_min[2]
    for j in 1:N[2]
      z = P_min[3]
      for k in 1:N[3]
        D[i,j,k, :] = [n,x,y,z] # Here it creates a node
        n += 1
        z += dz
      end
      y += dy
    end
    x += dx
  end
  return D
end

function _reset(self::FluidDomain)
  self.sol=Dict()
  self.field_types=Dict()
end

"Calculates the domain in the global reference frame"
function _calculate_D(self::FluidDomain)
  _D = self._D  # Domain in local ref frame

  N = self.NDIVS + 1 # Number of nodes on each direction
  D = zeros(N[1], N[2], N[3], 4)  # Domain in global ref frame

  for i in 1:N[1]
    for j in 1:N[2]
      for k in 1:N[3]
        _P = _D[i,j,k, :]         # Point in local ref frame
        n = _P[1]                  # Index of the point
                                  # Point in global ref frame
        P = countertransform(_P[2:end], self.invOaxis, self.O)

        D[i,j,k, :] = vcat(n, P) # Saves it
      end
    end
  end

  self.D = D
end
##### END OF FLUIDDOMAIN CLASS #################################################




################################################################################
# ALGEBRA
################################################################################
"""
Rotates and translates the vector V.

Receives the i', j', k' unit vectors of an euclidean system with origin T, and
returns V'. (In this version, the unit vectors have been organized as a matrix
M)
"""
function transform(V::typeof(Float64[]),
                    M::Array{Float64,2}, T::typeof(Float64[]))
  return M*(V-T)
end

function transform(Vs::Array{Array{Float64,1},1},
                    M::Array{Float64,2}, T::typeof(Float64[]))
  out = Array{Float64,1}[]
  for V in Vs
    push!(out, transform(V, M, T))
  end
  return out
end

"""
Rotates and translates back a vector V' that had been rotated and translated
into the system (i', j', k') with origin T, and returns the original V.
To ease repetitive computation, instead of giving the unit vectors, give the
inverse of their matrix.
"""
function countertransform(Vp::typeof(Float64[]),
                          invM::Array{Float64,2}, T::typeof(Float64[]))
  return invM*Vp + T
end

function countertransform(Vps::Array{Array{Float64,1},1},
                          invM::Array{Float64,2}, T::typeof(Float64[]))
  out = Array{Float64,1}[]
  for Vp in Vps
    push!(out, countertransform(Vp, invM, T))
  end
  return out
end

"Checks that the unit vectors given as the matrix M=[i;j;k] define a coordinate
system"
function check_coord_sys(M::Array{Float64,2}; raise_error::Bool=true)
  # Checks normalization
  for i in 1:size(M)[1]
    if abs(norm(M[i,:])-1) > 0.00000001
      println(M)
      if raise_error
        error("Not unitary axis: $(M[i,:])")
      else
        return false
      end
    end
  end

  # Checks ortogonality
  for i in size(M)[1]
    xi = M[i, :]
    xip1 = M[(i%size(M)[1])+1, :]
    proj = abs(dot(xi, xip1))
    if proj>0.00000001
      if raise_error
        error("Non-ortogonal system")
      else
        return false
      end
    end
  end
  return true
end

function check_coord_sys(M::Array{Array{Float64,1},1}; raise_error::Bool=true)
  dims = 3
  newM = zeros(dims,dims)
  for i in 1:dims
    newM[i, :] = M[i]
  end
  return check_coord_sys(newM; raise_error=raise_error)
end

end # END OF MODULE
