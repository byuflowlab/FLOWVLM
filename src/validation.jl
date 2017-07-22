# Verification and validation cases for FLOWVLM

include("FLOWVLM.jl")
vlm = FLOWVLM

using PyPlot


"""
Planar, 45-deg swept-back wing (Bertin's planar wing in Example 7.2, pp. 343
of Bertin's *Aerodynamics for Engineers*). Verified with Bertin's hand
calculations and validated with experimental data from Weber and Brebner (1958),
*Low-speed tests on 45-deg swept-back wings, part I*, Tables 3 and 4.
"""
function planarWing()

  # -------------------------------------------------
  # --- LIFT AND DRAG AS A FUNCTION OF ALPHA
  magVinf = 163*0.3048 # m/s
  rhoinf = 9.093/10^1

  alphas = [2.1, 4.2, 6.3, 8.4, 10.5];
  b=98*0.0254
  lambda = 45.0
  ar = 5.0
  tr = 1.0
  gamma = 0.0
  n=4*2^4
  twist = 0.0

  wing = nothing
  CLs, CDs = [], [];
  for alpha in alphas
    function Vinf(X, t)
      return magVinf*[ cos(alpha*pi/180), 0.0, sin(alpha*pi/180)]
    end
    wing = vlm.simpleWing(b, ar, tr, twist, lambda, gamma; n=n, r=1.0)
    vlm.solve(wing, Vinf)
    vlm.calculate_field(wing, "CFtot"; S=b^2/ar)

    info = vlm.fields_summary(wing)
    push!(CLs, info["CL"]);
    push!(CDs, info["CD"]);
  end

  # Weber's experimental data (Table 4)
  data = [[2.1, 4.2, 6.3, 8.4, 10.5],
            [0.121, 0.238, 0.350, 0.456, 0.559],
            [nothing, 0.005, 0.012, 0.022, 0.035]] # [alpha, Cl, Cd]

  fig = figure("validation_planar_wing",figsize=(7*2,5*2))
  # LIFT
  subplot(221)
  plot([0, maximum(alphas)], 0.0601*[0, maximum(alphas)], "-k",
        label="Bertin's hand calculation")
  plot(data[1], data[2], "ok",
        label="Weber's experimental data")
  plot(alphas, CLs, "or", label="MyVLM")

  xlim([0,12]);
  xlabel("Angle of attack (deg)");
  ylim([0,  maximum([ 0.0601*maximum(alphas), maximum(CLs)  ])*1.1]);
  ylabel("CL");
  grid(true, color="0.8", linestyle="--")
  legend(loc="best");
  title("Planar, 45-deg swept-back wing");

  # DRAG
  subplot(222)
  plot(data[1], data[3], "ok",
        label="Weber's experimental data")
  # plot(alphas, CDs*0.035/0.00035, "or", label="MyVLM")
  plot(alphas, CDs, "or", label="MyVLM")

  xlim([0,12]);
  xlabel("Angle of attack (deg)");
  ylim( [0,  maximum(  [ maximum(data[3][2:end]), maximum(CDs) ]  )*1.1]  );
  ylabel("CD");
  grid(true, color="0.8", linestyle="--")
  legend(loc="best");
  title("Planar, 45-deg swept-back wing");



  # =================================================

  # -------------------------------------------------
  # --- LIFT DISTRIBUTION AT ALPHA=4.2
  alpha = 4.2
  function Vinf(X,t)
    return magVinf*[ cos(alpha*pi/180), 0.0, sin(alpha*pi/180)]
  end
  wing = vlm.simpleWing(b, ar, tr, twist, lambda, gamma; n=n, r=1.0)
  vlm.solve(wing, Vinf)
  vlm.calculate_field(wing, "CFtot"; S=b^2/ar)
  vlm.calculate_field(wing, "Cftot/CFtot"; S=b^2/ar)
  y2b = 2*wing._ym/b
  ClCL = wing.sol["Cl/CL"]

  # Weber's data (Table 3)
  web_2yb = [0.0, 0.041, 0.082, 0.163, 0.245, 0.367, 0.510, 0.653, 0.898, 0.949]
  web_Cl = [0.235, 0.241, 0.248, 0.253, 0.251, 0.251, 0.251, 0.246, 0.192, 0.171]
  web_CL = 0.238
  web_ClCL = web_Cl/web_CL

  subplot(223)
  plot(y2b, ClCL, "or", label="MyVLM")
  plot(web_2yb, web_ClCL, "ok", label="Weber's experimental data")

  xlim([0,1]);
  xlabel(L"$\frac{2y}{b}$");
  ylim([minimum([minimum(ClCL), minimum(web_ClCL), 0]),
        maximum([maximum(ClCL), maximum(web_ClCL)])*1.1]);
  ylabel(L"$\frac{Cl}{CL}$");
  grid(true, color="0.8", linestyle="--")
  legend(loc="best");
  title(L"Spanwise lift distribution at $\alpha=4.2^\circ$");

  # --- DRAG DISTRIBUTION AT ALPHA=4.2
  CdCD = wing.sol["Cd/CD"]

  # Weber's data (Table 3)
  web_2yb = [0.0, 0.041, 0.082, 0.163, 0.245, 0.367, 0.510, 0.653, 0.898, 0.949]
  web_Cd = [0.059, 0.025, 0.016, 0.009, 0.007, 0.006, 0.006, 0.004, -0.002, -0.007]
  web_CD = 0.005
  web_CdCD = web_Cd/web_CD

  subplot(224)
  plot(y2b, CdCD, "or", label="MyVLM")
  plot(web_2yb, web_CdCD, "ok", label="Weber's experimental data")

  xlim([0,1]);
  xlabel(L"$\frac{2y}{b}$");
  ylim([minimum([minimum(CdCD), minimum(web_CdCD), 0]),
        maximum([maximum(CdCD), maximum(web_CdCD)])*1.1]);
  ylabel(L"$\frac{Cd}{CD}$");
  grid(true, color="0.8", linestyle="--")
  legend(loc="best");
  title(L"Spanwise drag distribution at $\alpha=4.2^\circ$");
  # # =================================================

  return wing

  # function Vinf(X, t)
  #   return magVinf*[ 0, 0.0, 1]
  # end
  #
  # wing = vlm.simpleWing(b, ar, cr, wing_alpha, lambda, gamma; n=n, r=1.0)
  # kp = -[ cos(alpha*pi/180), 0.0, -sin(alpha*pi/180)]
  # ip = [ sin(alpha*pi/180), 0.0, cos(alpha*pi/180)]
  # jp = [0.0,1.0,0.0]
  # T = [100.0,50,200]
  # vlm.setcoordsystem(wing, T, [ip,jp,kp])
  # vlm.solve(wing, Vinf)
  # vlm.calculate_field(wing, "D"; rhoinf=rhoinf)
  # vlm.calculate_field(wing, "CD"; S=b^2/ar)
  # vlm.save(wing, "test"; save_horseshoes=true)
  # return wing

end
