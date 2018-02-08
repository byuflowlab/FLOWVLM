import ForwardDiff
import Snopt
import JLD

vlm_path = "../"
include(vlm_path*"src/FLOWVLM.jl")
vlm = FLOWVLM

save_path = "../temps/opt_bwing03/"  # Output folder
run_name = "opt_bwing"
paraview = true                   # Calls paraview when done
prompt = true                     # Whether to prompt the user
verbose = true                    # Prints verbose on function calls
verbosetype = 2                   # 1:Short, 2:Shorter 3:Long
stepsverbose = 1                  # Steps between verbose


# ------------------- SIMULATION PARAMETERS ------------------------------------
# Wing Parameters (Bertin's wing)
b = 98*0.0254               # (m) span
AR = 5.0                    # Span over tip chord
Sref=b^2/AR                 # Reference area for coefficients
n = 100                      # Lattices in semi-span
nc = 50                     # Chord positions along semi-span

# Freestream
AOA = 4.2                   # (deg) wing angle of attack
magVinf = 163*0.3048        # (m/s) freestream
Vinf(X,t) = magVinf*[1,0,0]
rhoinf = 9.093/10^1         # (kg/m^3) air density

# Chords
pos = 1.0*[i for i in 0:1/(nc-1):1]   # Position of each chord along semi-span
clen = 1.0*[1 for i in 1:nc]      # Length of each chord as a fraction of tip chord
twist = 0.0*[1 for i in 1:nc]     # (deg) twist at each chord
sweep = 40.0*[1 for i in 1:nc-1]  # (deg) sweep of each section
dihed = 0.0*[1 for i in 1:nc-1]   # (deg) dihedral of each section

# ------------------- OPTIMIZATION PARAMETERS ----------------------------------
fobj_nin = nc-1             # Number of input variables on objective (same than x0)
fcomp_nin = fobj_nin        # Number of input variables on computation function
fcomp_nout = 2              # Dimension ouput of computation function to differentiate

x0 = twist[2:end]                     # Initial guess
lb = -10.0*[1 for i in 1:fobj_nin]      # Lower bounds
ub = 15.0*[1 for i in 1:fobj_nin]       # Upper bounds

options = Dict{String, Any}(          # SNOPT options
        "Derivative option"            => 1,
        "Verify level"                 => 1,
        "Major optimality tolerance"   => 1e-6
    )





# ------------------- OPTIMIZATION SETUP ---------------------------------------
# Constraints
# (See SNOPTfun)

# Storage
global fcalls = 0           # Current function calls
Xs = []                     # Optimization path
fs = []                     # Objective along path
gs = []                     # Gradient along path

# Creates save path
if save_path!=nothing
  vlm.vtk.create_path(save_path, prompt)
  _, this_file = splitdir(@__FILE__)
  run(`cp $this_file $(joinpath(save_path,this_file))`)   # Saves this run file
end

"Computation functions: compute here objective and some constrains"
function funs(x; output_vlm=true)
# Induced drag as a function of twist distribution on Bertin's wing
    this_twist = vcat(twist[1], x)   # Forces the root to be fixed

    # Generates the wing
    wing = vlm.complexWing(b, AR, n, pos, clen, this_twist, sweep, dihed)

    # Positions it at the angle of attack
    vlm.setVinf(wing, Vinf)
    M = vlm.vtk.rotation_matrix(0.0, -AOA, 0.0)
    vlm.setcoordsystem(wing, zeros(3), M)

    # Solves the lattice
    vlm.solve(wing, Vinf)

    # Calculates induced drag
    vlm.calculate_field(wing, "CFtot"; S=Sref)

    # Saves the wing
    if output_vlm && save_path!=nothing
      vlm.calculate_field(wing, "Cftot/CFtot"; S=Sref)
      vlm.save(wing, run_name; path=save_path, num=fcalls)
    end

    info = vlm.fields_summary(wing)
    CD = info["CD"]
    CL = info["CL"]
    return [CD, CL] # NOTE: Function output must be an array for ForwardDiff
end

"Gradients of computation functions"
function grads(x)
    if size(x,1)!=fcomp_nin
      error("Invalid number of variables. Expected $fcomp_nin, got $(size(x,1))")
    end

    # Function to differentiate
    fun(x) = funs(x; output_vlm=false) # fun[1] is CD, fun[2] is CL

    # nin = size(x,1)                       # Number of variables
    nout = fcomp_nout                       # Dimensions of function output
    nin = fcomp_nin

    dfdx = zeros(nout, nin)
    cfg = ForwardDiff.JacobianConfig(nothing, x,  ForwardDiff.Chunk{nin}())
    ForwardDiff.jacobian!(dfdx, fun, x, cfg)

    return dfdx
end

# SNOPT-wrapped function
function SNOPTfun(x)
    if size(x,1)!=fobj_nin
      error("Invalid number of variables. Expected $fobj_nin, got $(size(x,1))")
    end

    # Evaluates computation function
    CD, CL = funs(x)

    # Evaluates gradients
    gradseval = grads(x)
    dCDdx = gradseval[1,:]
    dCLdx = gradseval[2,:]

    f = CD                # Optimization objective
    g = dCDdx             # Optimization gradient
    fail = false          # Fail flag

    # Constrains
    con1 = -CL + 0.25                         # Minimum CL
    conxmins = [-(xi+AOA)-10 for xi in x]     # Min local AOA
    conxmaxs = [(xi+AOA)-15 for xi in x]      # Max local AOA
    ncons = 1 + 2*fobj_nin                    # Number of constraints
    c = vcat(con1, conxmins, conxmaxs)

    # Constrains gradients
    dcdx = zeros(ncons, fobj_nin)
    dcdx[1,:] = -dCLdx   # con1
    for i in 1:fobj_nin # conxmins
      for j in 1:fobj_nin
        dcdx[1+i, j] = -Int64(j==i)
      end
    end
    for i in 1:fobj_nin # conxmaxs
      for j in 1:fobj_nin
        dcdx[1+fobj_nin+i, j] = Int64(j==i)
      end
    end

    # Saves path
    global fcalls += 1
    push!(Xs, deepcopy(x));
    push!(fs, deepcopy(f));
    push!(gs, deepcopy(g));

    # Verbose
    if verbose
      if verbosetype==1
        if fcalls==1
          println("-----------------------------------------------------------")
          println("#fcall\tx\t\t\tf(x)\tc(x)")
          println("-----------------------------------------------------------")
        end
        if fcalls%stepsverbose==0
          println("$fcalls\t$(x)\t$(f)\t$(c)")
        end

      elseif verbosetype==2
        if fcalls==1
          println("-----------------------------------------------------------")
          println("#fcall\tx\t\t\tf(x)\tc(x)")
          println("-----------------------------------------------------------")
        end
        if fcalls%stepsverbose==0
          println("$fcalls\t$(round.(x,3))\t$(round.(f,3))\t$(round.(c,3))")
        end

      elseif verbosetype==3
        if fcalls%stepsverbose==0
          println("-----------------------------------------------------------")
          println("\t\tFUNCTION CALL #$fcalls")
          println("-----------------------------------------------------------")
          println("\tx=\n\t\t$x")
          println("\tf=\n\t\t$f")
          println("\tg=\n\t\t$g")
          println("\tc=\n\t\t$c")
          for i in 1:size(dcdx,1)
            println("\tdc$(i)dx=\n\t\t$(dcdx[i,:])")
          end
        end
      else
        error("Invalid verbose type $verbosetype.")
      end
    end

    return f, c, g, dcdx, fail
end


# ------------------- RUN OPTIMIZATION -----------------------------------------
xopt, fopt, info = Snopt.snopt(SNOPTfun, deepcopy(x0), lb, ub, options)




# ------------------- POSTPROCESSING -------------------------------------------
# Report
println("***********************************************************")
println("\t\tSUMMARY")
println("***********************************************************")
println(info)
println("\tFunction calls: $fcalls")
println("\tx0: $x0")
println("\txopt: $xopt")
println("\tfopt: $fopt")


# Saves optimization path
if save_path!=nothing
  JLD.save(joinpath(save_path,run_name*".jld"), "Xs", Xs, "fs", fs, "gs", gs)
end

# Calls paraview
if save_path!=nothing && paraview
  run(`paraview --data=$(save_path*run_name)_vlm...vtk`)
end
