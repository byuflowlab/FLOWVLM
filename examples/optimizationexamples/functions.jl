# A bunch of functions for processing the optimization examples

using PyPlot


"Receives both initial guess `x0` and optimized configuration `xopt` from the
optimization run, and compares this with Bertin's wing analytical (*Aerodynamics
for Engineers* Example 7.2, pp. 343) and experimental data (Weber and Brebner,
1958, *Low-speed tests on 45-deg swept-back wings, part I*, Tables 3 and 4)."
function compareBertins(x0, xopt, compfuns; verbose=true)

    # Generates x0 and xopt wings
    wings = Any[]
    if verbose
      println("Generating wings...")
      @time w0_CD, w0_CL     = compfuns(x0; output_vlm=false, output_wing=wings)
      @time wopt_CD, wopt_CL = compfuns(xopt; output_vlm=false, output_wing=wings)
    else
      w0_CD, w0_CL     = compfuns(x0; output_vlm=false, output_wing=wings)
      wopt_CD, wopt_CL = compfuns(xopt; output_vlm=false, output_wing=wings)
    end
    wing0 = wings[1]
    wingopt = wings[2]


    # -------------------------------------------------
    # --- AERODYNAMIC COEFFICIENTS AT ALPHA=4.2
    # Weber's experimental data (Table 4)
    web_CL = 0.238
    web_CD = 0.005

    println("\nDRAG COEFFICIENT")
    println("\tBertin's experimental CD:\t$(web_CD)")
    println("\tBaseline FLOWVLM CD:\t\t$(w0_CD)")
    println("\tOptimized FLOWVLM CD:\t\t$(wopt_CD)")
    println("\nLIFT COEFFICIENT")
    println("\tBertin's experimental CL:\t$(web_CL)")
    println("\tBaseline FLOWVLM CL:\t\t$(w0_CL)")
    println("\tOptimized FLOWVLM CL:\t\t$(wopt_CL)")


    # -------------------------------------------------
    # --- LIFT DISTRIBUTION AT ALPHA=4.2
    alpha = 4.2
    y2b_0 = 2*wing0._ym/b
    y2b_opt = 2*wingopt._ym/b
    ClCL_0 = wing0.sol["Cl/CL"]
    ClCL_opt = wingopt.sol["Cl/CL"]

    # Weber's data (Table 3)
    web_2yb = [0.0, 0.041, 0.082, 0.163, 0.245, 0.367, 0.510, 0.653, 0.898, 0.949]
    web_Cl = [0.235, 0.241, 0.248, 0.253, 0.251, 0.251, 0.251, 0.246, 0.192, 0.171]
    web_ClCL = web_Cl/web_CL


    fig = figure("BertinWingComparison",figsize=(7*2,5*1))
    subplot(121)
    plot(web_2yb, web_ClCL, "ok", label="Experimental")
    plot(y2b_0, ClCL_0, "--.b", label=L"FLOWVLM $x_0$")
    plot(y2b_opt, ClCL_opt, ".-.r", label=L"FLOWVLM $x_{opt}$")

    xlim([0,1]);
    xlabel(L"$\frac{2y}{b}$");
    ylim([minimum([minimum(ClCL_0), minimum(ClCL_opt), minimum(web_ClCL), 0]),
          maximum([maximum(ClCL_0), maximum(ClCL_opt), maximum(web_ClCL)])*1.1]);
    ylabel(L"$\frac{Cl}{CL}$");
    grid(true, color="0.8", linestyle="--")
    legend(loc="best");
    title(L"Spanwise lift distribution at $\alpha=4.2^\circ$");

    # --- DRAG DISTRIBUTION AT ALPHA=4.2
    CdCD_0 = wing0.sol["Cd/CD"]
    CdCD_opt = wingopt.sol["Cd/CD"]

    # Weber's data (Table 3)
    web_2yb = [0.0, 0.041, 0.082, 0.163, 0.245, 0.367, 0.510, 0.653, 0.898, 0.949]
    web_Cd = [0.059, 0.025, 0.016, 0.009, 0.007, 0.006, 0.006, 0.004, -0.002, -0.007]
    web_CdCD = web_Cd/web_CD

    subplot(122)
    plot(web_2yb, web_CdCD, "ok", label="Experimental")
    plot(y2b_0, CdCD_0, "--.b", label=L"FLOWVLM $x_0$")
    plot(y2b_opt, CdCD_opt, ".-.r", label=L"FLOWVLM $x_{opt}$")

    xlim([0,1]);
    xlabel(L"$\frac{2y}{b}$");
    # ylim([minimum([minimum(CdCD_0), minimum(CdCD_opt), minimum(web_CdCD), 0]),
    #       maximum([maximum(CdCD_0), maximum(CdCD_opt), maximum(web_CdCD)])*1.1]);
    ylim([minimum([minimum(web_CdCD), 0]),
          maximum([maximum(web_CdCD)])*1.1]);
    ylabel(L"$\frac{Cd}{CD}$");
    grid(true, color="0.8", linestyle="--")
    legend(loc="best");
    title(L"Spanwise drag distribution at $\alpha=4.2^\circ$");
    # # =================================================
end






using PyPlot
import Roots

"Plots a design space of only two dimensions"
function plot_space(f, xmin,xmax, ymin,ymax, n;
                    Xs=nothing, fs=nothing, cons=nothing,
                    xlbl=L"x", ylbl=L"y", zlbl=L"f(x,y)",
                    title_str="Objective function")
  x = linspace(xmin, xmax, n)
  y = linspace(ymin, ymax, n)

  xgrid = repmat(x',n,1)
  ygrid = repmat(y,1,n)

  z = zeros(n,n)
  for i in 1:n
      for j in 1:n
          z[j,i] = f([x[i],y[j]])[1]
      end
  end
#     zmin, zmax = minimum(z), maximum(z)


  # ----------------- Objective plots --------------------
  fig = figure("opt1", figsize=(7*2,5))

  # Surface plot
  subplot(121)
  ax1 = fig[:add_subplot](1,2,1, projection = "3d")
  title(title_str)
  ## Surface plot of objective
  ax1[:plot_surface](xgrid, ygrid, z,
              rstride=2,edgecolors="k", cstride=2,
              cmap=ColorMap("coolwarm"), alpha=0.5,
              linewidth=0.05)
  ## Optimization path
  if Xs!=nothing
      ax1[:plot]([this_X[1] for this_X in Xs],
                  [this_X[2] for this_X in Xs],
                  [f(x)[1] for x in Xs], "--ok")
      ax1[:scatter3D]([Xs[end][1]], [Xs[end][2]],
          [f(Xs[end])[1]], marker="*", c="y", label="Optimum")
  end
  xlabel(xlbl)
  ylabel(ylbl)
  zlabel(zlbl)
  xlim([xmin, xmax])
  ylim([ymin, ymax])
#     zlim([zmin, zmax])

  # Contour plot
  subplot(122)
  ax2 = fig[:add_subplot](1,2,2)
  cp = ax2[:contour](xgrid, ygrid, z, 15,
                  colors="black", linewidth=2.0)
  ax2[:clabel](cp, inline=1, fontsize=10)
  if Xs!=nothing
      ax2[:plot]([this_X[1] for this_X in Xs],
                  [this_X[2] for this_X in Xs], "-ok")
      ax2[:plot]([Xs[end][1]], [Xs[end][2]], "*y", label="Optimum")
  end
  xlabel(xlbl)
  ylabel(ylbl)
  xlim([xmin, xmax])
  ylim([ymin, ymax])
  tight_layout()

  # Constraints cons[i](x) < 0
  if cons!=nothing
      plot_flag = false
      ax2[:plot]([], []) # Dummy to match coloring
      for (i,con) in enumerate(cons)  # Iterates over constrains
          line_x, line_y, line_z = [], [], []
          for xi in x  # Iterates over x finding y cons[i](x,y) = 0
              wrap_con(yi) = con([xi,yi])[1]
              try
                  yi = Roots.fzero(wrap_con, ymin, ymax)
                  push!(line_x, xi)
                  push!(line_y, yi)
                  push!(line_z, f([xi,yi])[1])
              catch e
                  nothing
              end
          end
          if size(line_x)[1]!=0
              if size(line_x)[1]==1 # Case c(X)=c(X1)
                  line_x, line_y, line_z = [], [], []
                  for yi in y
                      wrap_con(xi) = con([xi,yi])[1]
                      try
                          xi = Roots.fzero(wrap_con, xmin, xmax)
                          push!(line_x, xi)
                          push!(line_y, yi)
                          push!(line_z, f([xi,yi])[1])
                      catch e
                          nothing
                      end
                  end
                  if size(line_x)[1]!=0
                      plot_flag=true
                      ax1[:plot](line_x, line_y, line_z, label="Cons. #$i")
                      ax2[:plot](line_x, line_y, label="Cons. #$i")
                  end
              else
                  plot_flag=true
                  ax1[:plot](line_x, line_y, line_z, label="Cons. #$i")
                  ax2[:plot](line_x, line_y, label="Cons. #$i")
              end
          end
      end
      if plot_flag
          ax1[:legend](loc="best")
          ax2[:legend](loc="best")
      end
  end
end

"Give the path of the optimizer and it will plot it"
function plot_opt(fs, gs; fig_name="opt_path")
  fig2 = figure(fig_name, figsize=(7*2,5*1))
  subplot(121)
  title("Convergence on f(X)")
  plot([i for i in 1:size(fs)[1]], fs, "--ok")
  ylabel("f(x)")
  xlabel("Iteration")
  grid(true, color="0.8", linestyle="--")
  subplot(122)
  title("Convergence on grad(f(X))")
  plot([i for i in 1:size(gs)[1]], [norm(this_g) for this_g in gs], "--ok")
  ylabel("|grad(f(x))|")
  xlabel("Iteration")
  grid(true, color="0.8", linestyle="--")
end

"Plots the design space of the pipeline problem at a fixed d value"
function design_space(d, lb, ub;
                      Xs=nothing, fs=nothing, cons=nothing)
  xmin, ymin = lb[1], lb[2]
  xmax, ymax = ub[1], ub[2]
  n = 100
  wrap_objf(x) = objf([x[1],x[2],d])
  plot_space(wrap_objf, xmin,xmax, ymin,ymax, n;
                  Xs=Xs, fs=fs, cons=cons,
                  xlbl=L"Flow velocity $V$",
                  ylbl=L"Pipe diameter $D$",
                  zlbl=L"Cost$(V,D,d)$",
          title_str="Objective function at d=$(round(d,6))")
end

function report(xopt,fopt)
  design_space(xopt[3], lb, ub; Xs=Xs, fs=fs, cons=cons)
  plot_opt(fs)
  println("***************************************************")
  println("*       CONSTRAINTS")
  println("***************************************************")
  print_constraints(xopt, cons)
  println("\n***************************************************")
  println("*       ANALYSIS FUNCTIONS")
  println("***************************************************")
  print_summary(xopt[1],xopt[2],xopt[3])
  println("\n***************************************************")
  println("*       OPTIMUM")
  println("***************************************************")
  println("\tFunction calls:\t\t$(size(Xs)[1])")
  println("\tFlow velocity V:\t$(round(xopt[1],1)) (ft/s)")
  println("\tPipe diameter D:\t$(round(xopt[2],2)) (ft)")
  println("\tParticle diameter d:\t$(round(xopt[3],4)) (ft)")
  println("---> Cost:\t\t\$$(round(fopt[1],0))")
end;
