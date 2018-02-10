# Induced drag minimization respect to chord length distribution and twist
# distribution using Bertin's wing as the baseline, subject to a specific lift
# coefficient.

import ForwardDiff
import Snopt
import JLD

modulepath, this_file = splitdir(@__FILE__)

vlm_path = joinpath(modulepath,"../../")
include(vlm_path*"src/FLOWVLM.jl")
vlm = FLOWVLM

include("functions.jl")

println("Defining parameters...")

save_path = joinpath(modulepath, "../../temps/opt_bwtwistchord00/")  # Output folder
run_name = "opt_bwtwistchord"
paraview = true                   # Calls paraview when done
save_horseshoes = false           # If saving VTKs, whether to save horseshoes
vis_des_space = true              # Plots design space when done (time consuming)
anim_des_space = false            # Animates the plot instead (VERY time consuming)
prompt = true                     # Whether to prompt the user
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
clen = 1.0*ones(nc)         # Length of each chord as a fraction of tip chord
twist = 0.0*ones(nc)        # (deg) twist at each chord
sweep = 45.0*ones(nc-1)     # (deg) sweep of each section
dihed = 0.0*ones(nc-1)      # (deg) dihedral of each section

# Calculations
trefftz = true              # Calculates induced drag at the Trefftz plane

# ------------------- OPTIMIZATION PARAMETERS ----------------------------------
fobj_nin = 2*nc             # Number of input variables on objective (same than x0)
fcomp_nin = fobj_nin        # Number of input variables on computation function
fcomp_nout = 2              # Dimension ouput of computation function to differentiate

x0 = vcat(twist, clen)      # Initial guess
lb = vcat(-(10+AOA)*ones(nc), 0.01*ones(nc))   # Lower bounds
ub = vcat((15-AOA)*ones(nc), 2.0*ones(nc))     # Upper bounds

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

println("Preparing optimization...")
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
function funs(x; output_vlm=true, output_wing=nothing, distrcoeffs=true)
# Induced drag as a function of chord length distribution on Bertin's wing
    this_twist = x[1:nc]
    this_clen = x[nc+1:end]

    # Generates the wing
    wing = vlm.complexWing(b, AR, n, pos, this_clen, this_twist, sweep, dihed)

    # Positions the wing at the right angle of attack
    vlm.setVinf(wing, Vinf)
    M = vlm.vtk.rotation_matrix(0.0, -AOA, 0.0)
    vlm.setcoordsystem(wing, zeros(3), M)

    # Solves the lattice
    vlm.solve(wing, Vinf)

    # Calculates induced drag
    vlm.calculate_field(wing, "CFtot"; S=Sref, lifting_interac=!trefftz)
    if distrcoeffs
      vlm.calculate_field(wing, "Cftot/CFtot"; S=Sref, lifting_interac=!trefftz)
    end

    # Saves the wing
    if output_vlm && save_path!=nothing
      vlm.save(wing, run_name; path=save_path, num=fcalls,
                                            save_horseshoes=save_horseshoes)
    end

    # Outputs the wing
    if output_wing!=nothing && typeof(output_wing)==Array{Any,1}
      push!(output_wing, wing)
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
    fun(x) = funs(x; output_vlm=false, distrcoeffs=false) # fun[1] is CD, fun[2] is CL

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

    f = CD*1e3                # Optimization objective
    g = dCDdx*1e3             # Optimization gradient
    fail = false              # Fail flag

    # Constrains
    con1 = -CL + 0.232
    ncons = 1                 # Number of nonlinear constraints
    c = [con1]

    # Constrains gradients
    dcdx = zeros(ncons, fobj_nin)
    dcdx[1,:] = -dCLdx   # con1

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

plot_opt(fs, gs)
if save_path!=nothing
  tight_layout()
  savefig(joinpath(save_path, run_name)*"_conv.png")
end


println("Plotting chord length distribution...")
fig = figure("chorddistr")
plot(1:nc, clen, "-k", label="Baseline")
plot(1:nc, xopt, "--r", label="Optimized")
xlim([1, nc])
ylim([minimum(lb), maximum(ub)])
xlabel("Chord number")
ylabel(L"c / c_{tip baseline}")
legend(loc="best")
grid(true, color="0.8", linestyle="--")
title("Semi-span chord length distribution")
if save_path!=nothing
  tight_layout()
  savefig(joinpath(save_path, run_name)*"_distr.png")
end


println("Saving results...")
# Saves optimization path
if save_path!=nothing
  # Saves variables
  JLD.save(joinpath(save_path,run_name*".jld"), "Xs", Xs, "fs", fs, "gs", gs,
                                                    "xopt", xopt, "fopt", fopt)
  # NOTE: For loading the variables, do:
  #       `Xs, fs, gs = JLD.load(save_path*run_name*".jld", "Xs", "fs", "gs");`
  # Saves SNOPT outputs
  for fl in ["snopt-summary.out", "snopt-print.out"]
    cp(fl, joinpath(save_path, fl))
  end
end

# Compares with Bertin's wing
println("Plotting comparison...")
compareBertins(x0, xopt, funs)
if save_path!=nothing
  tight_layout()
  savefig(joinpath(save_path, run_name)*"_bertincomp.png")
end

# Calls paraview
println("Calling paraview...")
if save_path!=nothing && paraview
  run(`paraview --data=$(save_path*run_name)_vlm...vtk`)
end

# Plots the design space
if vis_des_space
  if anim_des_space && save_path!=nothing
    println("Creating design space animation...")
    design_space_animation(save_path, run_name, funs, Xs, xopt, lb, ub;
                    verbose=true, ndiscr=20, alive=true, lbl_add_val=false,
                    output_vlm=false, output_wing=nothing, distrcoeffs=false)
  else
    println("Visualizing design space...")
    design_space(funs, Xs, xopt, lb, ub; ndiscr=25, alive=true, lbl_add_val=false,
                    output_vlm=false, output_wing=nothing, distrcoeffs=false)
    if save_path!=nothing
      tight_layout()
      savefig(joinpath(save_path, run_name)*"_path.png")
    end
  end
end
