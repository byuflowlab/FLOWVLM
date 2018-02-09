#   Induced drag minimization respect to twist distribution using Bertin's wing
# as the baseline, subject to a specific lift coefficient and local angle of
# attacks not exceeding separation angles.


import ForwardDiff
import Snopt
import JLD

include("functions.jl")

modulepath, this_file = splitdir(@__FILE__)

vlm_path = joinpath(modulepath,"../../")
include(vlm_path*"src/FLOWVLM.jl")
vlm = FLOWVLM

save_path = joinpath(modulepath, "../../temps/opt_bwingPROF00/")  # Output folder
run_name = "opt_bwing"
paraview = true                   # Calls paraview when done
prompt = false                     # Whether to prompt the user
verbose = true                    # Prints verbose on function calls
verbosetype = 2                   # 1:Short, 2:Shorter, 3:Long
stepsverbose = 1                  # Steps between verbose


# ------------------- SIMULATION PARAMETERS ------------------------------------
# Wing Parameters (Bertin's wing)
b = 98*0.0254               # (m) span
AR = 5.0                    # Span over tip chord
Sref=b^2/AR                 # Reference area for coefficients
n = 100                     # Lattices in semi-span
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
sweep = 45.0*[1 for i in 1:nc-1]  # (deg) sweep of each section
dihed = 0.0*[1 for i in 1:nc-1]   # (deg) dihedral of each section

# Calculations
trefftz = true              # Calculates induced drag at the Trefftz plane

# ------------------- OPTIMIZATION PARAMETERS ----------------------------------
fobj_nin = nc-1             # Number of input variables on objective (same than x0)
fcomp_nin = fobj_nin        # Number of input variables on computation function
fcomp_nout = 2              # Dimension ouput of computation function to differentiate

x0 = twist[2:end]                     # Initial guess
lb = -(10.0+AOA)*ones(fobj_nin)       # Lower bounds
ub = (15.0-AOA)*ones(fobj_nin)        # Upper bounds

options = Dict{String, Any}(          # SNOPT options
        "Scale option"                 => 1,
        "Derivative option"            => 1,
        "Verify level"                 => 0,
        "Major optimality tolerance"   => 1e-7,
        "Minor optimality tolerance"   => 1e-8,
        "Major feasibility tolerance"  => 1e-7,
        "Minor feasibility tolerance"  => 1e-8,
        "Scale tolerance"              => .95,
        "Print frequency"              => stepsverbose
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
  # Saves this run file
  run(`cp $(joinpath(modulepath,this_file)) $(joinpath(save_path,this_file))`)
end

"Computation functions: compute here objective and some constrains"
function funs(x; output_vlm=true, output_wing=nothing)
println("\t\tRabbit1.1: Generate wing")
prev_t = time()
# Induced drag as a function of twist distribution on Bertin's wing
    this_twist = vcat(twist[1], x)   # Forces the root to be fixed

    # Generates the wing
    wing = vlm.complexWing(b, AR, n, pos, clen, this_twist, sweep, dihed)
println("\t\t\ttime: $(time()-prev_t) s\n")
println("\t\tRabbit1.2: Set coords")
prev_t = time()

    # Positions it at the angle of attack
    vlm.setVinf(wing, Vinf)
    M = vlm.vtk.rotation_matrix(0.0, -AOA, 0.0)
    vlm.setcoordsystem(wing, zeros(3), M)

println("\t\t\ttime: $(time()-prev_t) s\n")
println("\t\tRabbit1.3: Solve VLM")
prev_t = time()
    # Solves the lattice
    vlm.solve(wing, Vinf)

println("\t\t\ttime: $(time()-prev_t) s\n")
println("\t\tRabbit1.4: Aerodynamic properties1")
prev_t = time()
    # Calculates induced drag
    vlm.calculate_field(wing, "CFtot"; S=Sref, lifting_interac=!trefftz)
println("\t\t\ttime: $(time()-prev_t) s\n")
println("\t\tRabbit1.4: Aerodynamic properties2")
prev_t = time()
    vlm.calculate_field(wing, "Cftot/CFtot"; S=Sref, lifting_interac=!trefftz)

println("\t\t\ttime: $(time()-prev_t) s\n")
println("\t\tRabbit1.5: Output VTK")
prev_t = time()
    # Saves the wing
    if output_vlm && save_path!=nothing
      vlm.save(wing, run_name; path=save_path, num=fcalls)
    end

println("\t\t\ttime: $(time()-prev_t) s\n")
println("\t\tRabbit1.6: Tail")
prev_t = time()
    # Outputs the wing
    if output_wing!=nothing && typeof(output_wing)==Array{Any,1}
      push!(output_wing, wing)
    end

    info = vlm.fields_summary(wing)
    CD = info["CD"]
    CL = info["CL"]
println("\t\t\ttime: $(time()-prev_t) s\n")
    return [CD, CL] # NOTE: Function output must be an array for ForwardDiff
end

# Function to differentiate
fun(x) = funs(x; output_vlm=false) # fun[1] is CD, fun[2] is CL

"Gradients of computation functions"
function grads(x)
    if size(x,1)!=fcomp_nin
      error("Invalid number of variables. Expected $fcomp_nin, got $(size(x,1))")
    end

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
println("Rabbit: SNOPT functionc all")
prev_t = time()
    if size(x,1)!=fobj_nin
      error("Invalid number of variables. Expected $fobj_nin, got $(size(x,1))")
    end

    # Evaluates computation function
println("\tRabbit1: Call compt function")
prev_t = time()
    CD, CL = funs(x)
println("\t\ttime: $(time()-prev_t) s\n")

    # Evaluates gradients
println("\tRabbit2: Call grad function")
prev_t = time()
    gradseval = grads(x)
println("\t\ttime: $(time()-prev_t) s\n")
    dCDdx = gradseval[1,:]
    dCLdx = gradseval[2,:]

    f = CD*1e3                # Optimization objective
    g = dCDdx*1e3             # Optimization gradient
    fail = false          # Fail flag

    # Constrains
    con1 = -CL + 0.232                        # Minimum CL
    # conxmins = [-(xi+AOA)-10 for xi in x]     # Min local AOA
    # conxmaxs = [(xi+AOA)-15 for xi in x]      # Max local AOA
    # ncons = 1 + 2*fobj_nin                    # Number of constraints
    # c = vcat(con1, conxmins, conxmaxs)
    ncons = 1                    # Number of constraints
    c = [con1]

    # Constrains gradients
    dcdx = zeros(ncons, fobj_nin)
    dcdx[1,:] = -dCLdx   # con1
    # for i in 1:fobj_nin # conxmins
    #   for j in 1:fobj_nin
    #     dcdx[1+i, j] = -Int64(j==i)
    #   end
    # end
    # for i in 1:fobj_nin # conxmaxs
    #   for j in 1:fobj_nin
    #     dcdx[1+fobj_nin+i, j] = Int64(j==i)
    #   end
    # end

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

println("\ttime: $(time()-prev_t) s\n")
    return f, c, g, dcdx, fail
end


# ------------------- RUN OPTIMIZATION -----------------------------------------
println("***********************************************************")
println("RUNNING OPTIMIZATION $save_path$run_name")
println("***********************************************************")
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
  for fl in ["snopt-summary.out", "snopt-print.out"]
    cp(fl, joinpath(save_path, fl))
  end
end

# Compares with Bertin's wing
compareBertins(x0, xopt, funs)

# Calls paraview
if save_path!=nothing && paraview
  run(`paraview --data=$(save_path*run_name)_vlm...vtk`)
end
