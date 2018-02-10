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

    # Elliptic load
    ell_2yb = linspace(0, 1, 101)
    ell_ClCL = maximum(ClCL_opt) * sqrt.(1-ell_2yb.^2)


    fig = figure("BertinWingComparison",figsize=(7*2,5*1))
    subplot(121)
    plot(web_2yb, web_ClCL, "ok", label="Experimental")
    plot(y2b_0, ClCL_0, "--.b", label=L"FLOWVLM $x_0$")
    plot(y2b_opt, ClCL_opt, ".r", label=L"FLOWVLM $x_{opt}$")
    plot(ell_2yb, ell_ClCL, "--k", label="Elliptic distribution")

    xlim([0,1]);
    xlabel(L"$\frac{2y}{b}$");
    ylim([minimum([minimum(ClCL_0), minimum(ClCL_opt), minimum(web_ClCL), 0]),
          maximum([maximum(ClCL_0), maximum(ClCL_opt), maximum(web_ClCL)])*1.1]);
    ylabel(L"$\frac{Cl}{CL}$");
    grid(true, color="0.8", linestyle="--")
    legend(loc="best");
    title(L"Spanwise lift distribution at $\alpha=4.2^\circ$")

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


"Plots a design space of only two dimensions"
function plot_space(f, xmin,xmax, ymin,ymax, n;
                    Xs=nothing, xopt=nothing, x_i=1, y_i=2,
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
      ax1[:plot]([x[x_i] for x in Xs], [x[y_i] for x in Xs],
                                      [f(x)[1] for x in Xs], "--.k")
      # ax1[:scatter3D]([x[x_i] for x in Xs], [x[y_i] for x in Xs],
      #                 [f(x)[1] for x in Xs], marker="o", c="k", label="Optimum")
      ax1[:scatter3D]([xopt[x_i]], [xopt[y_i]], [f(xopt)[1]],
                                      s=100, marker="*", c="r", label="Optimum")
  end
  xlabel(xlbl)
  ylabel(ylbl)
  zlabel(zlbl)
  xlim([xmin, xmax])
  ylim([ymin, ymax])

  # Contour plot
  subplot(122)
  ax2 = fig[:add_subplot](1,2,2)
  cp = ax2[:contour](xgrid, ygrid, z, 15,
                  colors="black", linewidth=2.0)
  ax2[:clabel](cp, inline=1, fontsize=10)
  if Xs!=nothing
      ax2[:plot]([x[x_i] for x in Xs], [x[y_i] for x in Xs], "-ok")
      ax2[:plot]([xopt[x_i]], [xopt[y_i]], "*r", label="Optimum")
  end
  xlabel(xlbl)
  ylabel(ylbl)
  xlim([xmin, xmax])
  ylim([ymin, ymax])
  tight_layout()
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

"""
Receives the computation function `compfun`, optimization path `Xs`, the
optimum `xopt`, and linear constraints `lb`,`ub`, and generates a three
dimensional plot on the two variables of maximum variation showing the
optimization path.

NOTE: It expects to find the objective in the first output of `compfun`.
"""
function design_space(compfun, Xs,
                        xopt::AbstractArray, lb::AbstractArray,
                        ub::AbstractArray; x_i::Int64=-1, y_i::Int64=-1,
                        ndiscr::Int64=15, lbl="\nxopt", lbl_add_val=true,
                        alive=false,
                        compfun_args...)

  # ---------- Chooses what variables to put in axes x and y -------------------
  x_is = [x_i, y_i]
  # Case of automatic choice
  if -1 in x_is
    n = size(xopt, 1) # Number of variables
    needed = size([i for i in x_is if i==-1],1) # Number of variables needed

    # Calculates variation range of each variable
    variations = [
    ( maximum([X[i] for X in Xs]) - minimum([X[i] for X in Xs]) )/(ub[i]-lb[i])
                                          for i in 1:n ]

    # Finds the two variables with max variation
    max1_i = indmax(variations)
    max2_i = indmax(vcat( variations[1:max1_i-1], variations[max1_i+1:end] ))

    # Saves selection
    if needed==1
      if x_i==-1; x_is[1]=max1_i; else; x_is[2]=max1_i; end;
    else
      x_is[1] = max1_i
      x_is[2] = max2_i
    end
  end

  _x_i, _y_i = x_is

  # Ranges for surface plotting
  xmin, ymin = lb[_x_i], lb[_y_i]
  xmax, ymax = ub[_x_i], ub[_y_i]

  # Wraps compfun as a 2D function
  fncalls = 0
  function wrap_compfun(x)
    fncalls += 1
    prev_t = time()
    if size(x,1)==2 # Case of generating surface and contour plots
      compfun_x = deepcopy(xopt)
      compfun_x[_x_i] = x[1]
      compfun_x[_y_i] = x[2]
      out = compfun(compfun_x; compfun_args...)[1]
    else  # Case of plotting the path
      compfun_x = deepcopy(xopt)
      compfun_x[_x_i] = x[_x_i]
      compfun_x[_y_i] = x[_y_i]
      out = compfun(compfun_x; compfun_args...)[1]
    end
    if alive; println("\tFunction call #$fncalls: $(time()-prev_t) (s)"); end;
    return out
  end

  ttl = lbl * (lbl_add_val ? "=$(signif.(xopt,3))" : "")

  plot_space(wrap_compfun, xmin, xmax, ymin, ymax, ndiscr;
                      Xs=Xs, xlbl="x$(_x_i)", ylbl="x$(_y_i)",
                      zlbl="f(x$(_x_i),y$(_y_i))",
                      title_str="Design space at"*ttl,
                      x_i=_x_i, y_i=_y_i, xopt=xopt)
end

function design_space_animation(save_path::String, run_name::String,
                                compfun, Xs, xopt, args...; verbose=true,
                                first=1, last=-1, ext=".png",
                                optargs...)

  # Points to iterate through
  _Xs = Xs[ first : (last==-1 ? size(Xs,1) : last) ]
  if last==-1 && Xs[end]!=xopt
    push!(_Xs, xopt)
  end

  # Iterates through the optimization path
  for (i, this_x) in enumerate(_Xs)
    if verbose; println("Plotting iteration $i..."); end;
    this_Xs = _Xs[1:i]
    design_space(compfun, this_Xs, this_x, args...;
                        lbl="ite $i\n x", optargs...)

    this_num = ( i<10 ? "000" : (i<100 ? "00" : (i<1000 ? "0" : "")) )*"$i"
    savefig(joinpath(save_path, run_name)*"."*this_num*ext)
    clf()
  end

  if verbose; println("\tDone plotting"); end;
end









#
