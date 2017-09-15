# Verification and validation cases for FLOWVLM

include("../src/FLOWVLM.jl")
vlm = FLOWVLM

using PyPlot
using JLD

################################################################################
# VALIDATION CASES
################################################################################
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
  b = 98*0.0254
  lambda = 45.0
  ar = 5.0
  tr = 1.0
  gamma = 0.0
  n = 4*2^4
  twist = 0.0
  r = 20.0
  central = false

  wing = nothing
  CLs, CDs = [], [];
  for alpha in alphas
    function Vinf(X, t)
      return magVinf*[ cos(alpha*pi/180), 0.0, sin(alpha*pi/180)]
    end
    wing = vlm.simpleWing(b, ar, tr, twist, lambda, gamma;
            n=n, r=r, central=central)
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
  plot(alphas, CLs, "or", label="FLOWVLM")

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
  # plot(alphas, CDs*0.035/0.00035, "or", label="FLOWVLM")
  plot(alphas, CDs, "or", label="FLOWVLM")

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
  wing = vlm.simpleWing(b, ar, tr, twist, lambda, gamma;
                          n=n, r=r, central=central)
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
  plot(y2b, ClCL, "or", label="FLOWVLM")
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
  plot(y2b, CdCD, "or", label="FLOWVLM")
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
end


"""
Validation of predicted moments compared to results on a Warren 12 wing
(tappered, swept, flat wing)
For details, see Validation 3 in *Great OWL Publishing - Surfaces; Vortex
Lattice Module*, pp. 109.
"""
function warren12()
  magVinf = 1.0
  rhoinf = 9.093/10^1
  qinf = (1/2)*rhoinf*magVinf^2
  r_cg = [0.0, 0.0, 0.0]
  S = 2.83*0.092903 # m^2
  n = 100
  barc = 1.0*0.3048 # m

  b=2.83*0.3048 # m
  lambda = 53.54 # deg
  # ar = 2.83
  ar = 5.66
  tr = 0.5/1.5
  gamma = 0.0
  twist = 0.0
  wing = vlm.simpleWing(b, ar, tr, twist, lambda, gamma; n=n, r=1.0)

  cls = []
  cms = []
  alphas = [i for i in 1:2:12]
  for alpha in alphas
    Vinf(X,t) = magVinf*[ cos(alpha*pi/180), 0.0, sin(alpha*pi/180)]

    vlm.solve(wing, Vinf)
    vlm.calculate_field(wing, "Ftot"; rhoinf=rhoinf)
    vlm.calculate_field(wing, "CFtot"; S=S)
    vlm.calculate_field(wing, "Mtot"; r_cg = r_cg)
    vlm.calculate_field(wing, "CMtot"; S=S, l="automatic", qinf=qinf)
    info = vlm.fields_summary(wing)
    push!(cls, info["CL"]);
    push!(cms, info["CMtot"]);
  end

  # Curves reported in the reference document
  CLalpha = 2.743*pi/180 # 1/deg
  CMalpha = -3.10*pi/180 # 1/deg

  fig = figure("validation_warren12",figsize=(7*2,5))

  # LIFT
  subplot(121)
  plot([0, maximum(alphas)], CLalpha*[0, maximum(alphas)], "-k",
        label="Published data")
  plot(alphas, cls, "or", label="FLOWVLM")

  xlim([0,maximum(alphas)+1]);
  xlabel("Angle of attack (deg)");
  ylim([0,  maximum([ CLalpha*maximum(alphas), maximum(cls)  ])*1.1]);
  ylabel("CL");
  grid(true, color="0.8", linestyle="--")
  legend(loc="best");
  title("Warren 12");

  # MOMENT
  subplot(122)
  plot([0, maximum(alphas)], CMalpha*[0, maximum(alphas)], "-k",
        label="Published data")
  plot(alphas, -cms, "or", label="FLOWVLM")

  xlim([0,maximum(alphas)+1]);
  xlabel("Angle of attack (deg)");
  ylim( [minimum(  [ CMalpha*maximum(alphas), minimum(-cms) ]  )*1.1, 0]  );
  ylabel("CM");
  grid(true, color="0.8", linestyle="--")
  legend(loc="best");
  title("Warren 12");

  return wing
end


"""
Validation on different wing properties (taper, sweep, and twist) compared
with experimental data from Anderson's 1940 *Determination of the
characteristics of tapered wings*. The data was captured on a variety of cross-
sectional airfoil geometries, but the current formulation of FLOWVLM doesn't
account for airfoil geometry, hence only the symmetric case (wing 00-xx-xx)
can be expected to match.

Since all other validation cases are on untwisted wings, this study is
particularly useful to validate the effects of twist.

Give it an array of the wings to calculate.
"""
function twist(; to_calculate = [1,2], save=false, n=30)
  # Sweep angles
  lambdas = [0.0, 0.0, 15.0, 30.0, 30.0, 15.0, 15.0, 15.0, 15.0]
  # Washout angles
  washouts = -[0.0, 0.0, 0.0, 0.0, 8.5, 8.5, 0.0, 3.45, 3.45]
  # Taper ratios
  trs = [(2/3)/(4/3), 2/4, 2/4, 2/4, 2/4, 2/4, 2/4, 2/4, 2/8]
  # Aspect ratio (b^2/S)
  AR = 6.0
  # Aspect ratio as defined in simpleWing() (b/c_tip)
  ars = 1./([2/3, 2/3, 2/3, 2/3, 2/3, 2/3, 2/3, 2/3, 2/5]/AR)

  # Since this wing align the quarter chords for lambda=0, a correction must be
  # made
  lambdas_corr = lambdas + atan(  ( 1/AR )*( (1-trs)./(1+trs) )  )*180/pi

  magVinf = 1.0
  b = 30.0
  S = b^2/AR

  # Angles of attack it will compute
  alphas = [i for i in -4:15]
  deleteat!(alphas, 5) # Gets ride of the null angle of attack

  wings = []
  wings_CLs = []
  wings_CDs = []
  for i in to_calculate
    wing = vlm.simpleWing(b, ars[i], trs[i], 0.0, lambdas_corr[i], 0.0;
                          twist_tip=washouts[i], n=n, r=1.0)

    cls = []
    cds = []
    for alpha in alphas
      Vinf(X,t) = magVinf*[ cos(alpha*pi/180), 0.0, sin(alpha*pi/180)]
      vlm.solve(wing, Vinf)
      vlm.calculate_field(wing, "CFtot"; S=S)
      info = vlm.fields_summary(wing)
      push!(cls, alpha<0 ? -info["CL"] : info["CL"]);
      push!(cds, info["CD"]);
    end
    push!(wings_CLs, cls)
    push!(wings_CDs, cds)

    if save
      push!(wings, wing)
      vlm.save(wing, "twist$i", save_horseshoes=true)
    end
  end

  ############# EXPERIMENTAL DATA
  exp_CLs = [ [[-4, -2, 4, 8, 12, 16],[-0.3, -0.15, 0.3, 0.6, 0.92, 1.22]], # Wing 1 [alphas, CLs]
              [[-4, -2, 4, 8, 12, 16],[-0.2, -0.05, 0.4, 0.7, 1.0, 1.3]], # Wing 2
              [[-4, -2, 4, 8, 12, 16],[-0.18, 0.0, 0.45, 0.75, 1.075, 1.33]], # Wing 3
              [[],[]], # Wing 4
              [[],[]], # Wing 5
              [[-4, -2, 4, 8, 12, 16],[-0.4, -0.25, 0.2, 0.5, 0.8, 1.1]], # Wing 6
              [[-4, -2, 4, 8, 12, 16],[-0.275, -0.1, 0.375, 0.68, 0.95, 1.25]], # Wing 7
              [[-4, -2, 4, 8, 12, 16],[-0.4, -0.225, 0.23, 0.54, 0.85, 1.125]], # Wing 8
              [[-4, -2, 4, 8, 12, 16],[-0.37, -0.2, 0.25, 0.55, 0.85, 1.125]], # Wing 9
            ]
  exp_CDs = [ [[-4, -2, 4, 8, 12, 16],[0.0175, 0.01, 0.0175, 0.0325, 0.0575, 0.0975]], # Wing 1 [alphas, CDs]
              [[-4, -2, 4, 8, 12, 16],[0.0125, 0.0075, 0.02, 0.04, 0.07, 0.12]], # Wing 2
              [[-4, -2, 4, 8, 12, 16],[0.01, 0.009, 0.02, 0.04, 0.07, 0.115]], # Wing 3
              [[],[]], # Wing 4
              [[],[]], # Wing 5
              [[-4, -2, 4, 8, 12, 16],[0.02, 0.0175, 0.015, 0.025, 0.048, 0.08]], # Wing 6
              [[-4, -2, 4, 8, 12, 16],[0.015, 0.01, 0.019, 0.035, 0.0625, 0.12]], # Wing 7
              [[-4, -2, 4, 8, 12, 16],[0.02, 0.0125, 0.0125, 0.025, 0.05, 0.085]], # Wing 8
              [[-4, -2, 4, 8, 12, 16],[0.0185, 0.01, 0.015, 0.025, 0.05, 0.09]], # Wing 9
            ]


  ############# PLOTS
  y_data = [wings_CLs, wings_CDs]
  data_exp = [exp_CLs, exp_CDs]
  y_labels = Dict(L"C_L"=>1, L"C_D"=>2)
  n_w = size(to_calculate)[1]
  fig = figure("validation_twisted",figsize=(7*2,5*n_w))

  n_p = 1
  for (i,key) in enumerate(to_calculate)

    for (y_lbl, y_lbl_i) in y_labels
      subplot(100*n_w + 20 + n_p)
      title("Wing #$key")

      xprmntl = data_exp[y_lbl_i][key]
      exp_x, exp_y = xprmntl
      plot(exp_x, exp_y, "--ok", label="Experimental")

      y = y_data[y_lbl_i][i]
      plot(alphas, y, "or", label="FLOWVLM")

      xlim([minimum(alphas),maximum(alphas)*1.25]);
      xlabel("Angle of attack (deg)");
      y_min = minimum([minimum(y), minimum(exp_y)])
      y_max = maximum([maximum(y), maximum(exp_y)])
      y_low = y_min - (y_max-y_min)*0.25
      y_up = y_max + (y_max-y_min)*0.25
      ylim([y_low, y_up]);
      ylabel(y_lbl);
      grid(true, color="0.8", linestyle="--")
      legend(loc="best");
      n_p += 1
    end
  end

  return wings
end


"""
Comparison with experimental data on planar wings at low reynolds numbers from
Ananda 2015, *Measured aerodynamic characteristics of wings at low Reynolds
numbers*.
It is expected that the induced drag calculated through FLOWVLM won't match
quite well their measured drag due to predominant viscous effects, hence,
this study is useful for seeing the discrepancy that is introduced by not
accounting for viscous drag.
"""
function planarwing_lowreynolds(; save_w=false)
  save = save_w
  save_fd = save && true
  _save_name = "lowRE"

  b = [7.0, 10.5, 14.0, 17.5, 6.95, 10.43, 13.91, 6.75, 10.13, 13.5,]*0.0254 #m
  AR = [2.0, 3.0, 4.0, 5.0, 2.0, 3.0, 4.0, 2.0, 3.0, 4.0]     # aspect ratios
  tr = [1.0, 1.0, 1.0, 1.0, 0.75, 0.75, 0.75, 0.5, 0.5, 0.5]  # taper ratios
  lambda = atan(   (  (1-tr)./(1+tr)  )./AR   ) * 180 / pi    # sweep
  bar_c = b./AR   # mean chord

  # Aspect ratio as defined in simpleWing() (b/c_tip)
  c_tips = 2*b./AR.*(tr./(1+tr))
  ars = b./c_tips

  mu = 1.983/10^5     # Pa*s , dynamic viscosity of air
  rhoinf = 1.204         # kg/m^3 , density of air
  nu = mu/rhoinf

  n = 50
  RE = [80*10^3]

  # -------------------------- Figure 6
   alphas_6 = [i for i in -15:15]
  deleteat!(alphas_6, 16)
  cls = []
  cds = []
  wing = nothing
  for alpha in alphas_6
    # println("Solving alpha=$alpha")
    i = 2
    re = 80*10^3
    magVinf = re*nu/bar_c[i]
    Vinf(X,t) = magVinf*[ cos(alpha*pi/180), 0.0, sin(alpha*pi/180)]
    S = b[i]*bar_c[i]

    # wing = tls.simpleWing(b[i], AR[i], tr[i], 0.0, lambda[i], 0.0; n=n, r=1.0)
    wing = vlm.simpleWing(b[i], ars[i], tr[i], 0.0, lambda[i], 0.0; n=n, r=1.0)
    vlm.solve(wing, Vinf)
    vlm.calculate_field(wing, "CFtot"; S=S)
    info = vlm.fields_summary(wing)
    push!(cls, info["CL"]*(alpha<0 ? -1 : 1))
    push!(cds, info["CD"])
    # println("\tDone!")
  end

  fig = figure("validation_planar_wing_low_re",figsize=(7*2,5))
  # LIFT
  subplot(121)
  plot(alphas_6, cls, "or", label="MyVLM")

  xlim([-20, 30]);
  xlabel("Angle of attack (deg)");
  ylim([-1, 1]);
  ylabel(L"$C_L$");
  grid(true, color="0.8", linestyle="--")
  #legend(loc="best");
  title(L"Flat Plate, $AR=3, tr=1, Re=30\times 10^3$");

  # DRAG
  subplot(122)
  plot(alphas_6, cds, "or", label="MyVLM")

  xlim([-20, 30]);
  xlabel("Angle of attack (deg)");
  ylim([0, 0.5]);
  ylabel(L"$C_D$");
  grid(true, color="0.8", linestyle="--")
  #legend(loc="best");
  title(L"Flat Plate, $AR=3, tr=1, Re=30\times 10^3$");

  return wing
end

"""
Validation of interaction between lifting surfaces. This case compares predicted
values of canard interference on wing with experimental results reported in
Blair's 1973 *Canard-wing lift interference related to maneuvering aircraft at
subsonic speeds*. The wing is has a low-aspect ratio of 2.5, and 4.7 for the
canard. The canard is located fore of the wing.

RESULTS: FAILED. Both curves are completely off, and there is no way to know
where it fails; it could be that the interaction is wrongly modeled on the
solver, or that the aspects ratios scape the limits of a VLM, or the Mach number
is too high (0.7), or simply that I modeled the geometry wrong or I'm using the
wrong planform areas. Hence, this was a waste of time. I needed a simpler
validation case.
"""
function canard_wing_interaction(; body=false, n=2,
                                  save_w=false, save_fdom=false)
  file_name = "cwinter"

  # Case
  Re_i = 1
  w_i = 1
  pivot_z_i = 1
  pivot_pos_i = 2
  dc_i = 3
  iden = "$Re_i$w_i$pivot_z_i$pivot_pos_i$dc_i"

  in2m = 0.0254
  rhoinf = 1.177            # kg/m^3, density of air at 1 atm and 300K
  muinf = 1.846/10^5        # kg/ms, dynamic viscosity of air at 1 atm and 300K

  # Experimental setup
  M = [0.7, 0.9]            # Mach number
  Re = [2.58, 2.91]*10^6    # Reynolds number
  # alphas = [i for i in -4:5:24]*pi/180      # Angles of attack
  alphas = [i for i in -4:1:24]*pi/180      # Angles of attack
  bar_c = 9.18*in2m         # Wing's mean geometric chord
  w_sweep = [60, 44]*pi/180 # Wing's leading edge sweep
  c_sweep = 51.7*pi/180     # Canard's leading edge sweep
  w_pos = [18.76, 20.48]*in2m             # Longitudinal position of wing
  # c_pos = [1.61, 1.14]*bar_c*in2m       # Longitudinal position of canard
  w_S = 160*in2m^2            # Wing's reference planform area
  c_S = 25.6*in2m^2           # Canard's reference planform area
  pivot = 2.67*in2m         # Distance of canard's pivot point from LE at root
  pivot_pos = [11.04, 15.36]*in2m         # Longitudinal position of pivot
  pivot_z = [0, 1.69]*in2m      # Vertical position of pivot
  dc = [0, 2.5, 5, 7.5, 10]*pi/180        # Canard deflections

  # Geometric setup
  ## Body
  b_l = 38*in2m             # Length of body
  b_width = 1.5*in2m        # Width of body
  b_eight = 4*in2m          # height of body
  b_sweep = 70*pi/180       # Body's ficticious sweep
  b_xroot = 0.0
  b_yroot = 0.0
  b_zroot = 0.0
  b_croot = b_l
  b_ytip = body ? b_width : 0.0
  b_xtip = b_ytip*tan(b_sweep)
  b_ztip = 0.0
  b_ctip = b_croot - b_xtip
  ## Wing
  w_ytip = 10*in2m
  w_ztip = 0.0
  w_ctip = 2.67*in2m
  w_yroot = b_ytip
  w_zroot = 0.0
  w_croot = 11.73*in2m
  ## Canard
  c_lambda = c_sweep
  c_ytip = 5.50*in2m
  c_ctip = 1.07*in2m
  c_yroot = body ? b_width : 0.0
  c_croot = 5.33*in2m

  # VLM parameters
  b_n = Int(ceil(n*5))
  w_n = Int(ceil(n*20))
  c_n = Int(ceil(n*15))
  b_r = 1.0
  w_r = 4.0
  c_r = 4.0


  c_CLs, m_CLs = [], []
  for (i,alpha) in enumerate(alphas)

    # Case-dependent experimental setup
    magVinf = Re[Re_i]*muinf/(rhoinf*bar_c)
    Vinf(X,t) = magVinf*[cos(alpha), 0, sin(alpha)]

    # Case-dependent geometric setup
    ## Wing
    w_O = [w_pos[w_i], 0, 0]
    w_Oaxis = [1.0 0 0; 0 1 0; 0 0 1]
    w_lambda = w_sweep[w_i]
    w_xroot = 0 - (b_width-w_yroot)*tan(w_lambda)
    w_xtip = 0 + (w_ytip-b_width)*tan(w_lambda)
    ## Canard
    defl = dc[dc_i]
    pivot_x = pivot_pos[pivot_pos_i]
    c_O = [pivot_x, 0, pivot_z[pivot_z_i]]
    c_Oaxis = [cos(defl) 0 -sin(defl); 0 1 0; sin(defl) 0 cos(defl)]
    c_xroot = -pivot - (b_width-c_yroot)*tan(c_lambda)
    c_zroot = 0.0
    c_xtip = -pivot + (c_ytip-b_width)*tan(c_lambda)
    c_ztip = 0.0

    # Building the WingSystem
    ## Body
    if body
      _body = vlm.Wing(b_xtip, -b_ytip, b_ztip, b_ctip, 0.0)
      vlm.addchord(_body, b_xroot, b_yroot, b_zroot, b_croot, 0.0, b_n; r=b_r)
      vlm.addchord(_body, b_xtip, b_ytip, b_ztip, b_ctip, 0.0, b_n; r=1/b_r)
    end
    ## Wing
    wing = vlm.Wing(w_xtip, -w_ytip, w_ztip, w_ctip, 0.0)
    vlm.addchord(wing, w_xroot, -w_yroot, w_zroot, w_croot, 0.0, w_n; r=w_r)
    if body
      vlm.addchord(wing, w_xroot, w_yroot, w_zroot, w_croot, 0.0, 2*b_n; r=1.0)
    end
    vlm.addchord(wing, w_xtip, w_ytip, w_ztip, w_ctip, 0.0, w_n; r=1/w_r)
    vlm.setcoordsystem(wing, w_O, w_Oaxis)
    ## Canard
    canard = vlm.Wing(c_xtip, -c_ytip, c_ztip, c_ctip, 0.0)
    vlm.addchord(canard, c_xroot, -c_yroot, c_zroot, c_croot, 0.0, c_n; r=c_r)
    if body
      vlm.addchord(canard, c_xroot, c_yroot, c_zroot, c_croot, 0.0, 2*b_n; r=1.0)
    end
    vlm.addchord(canard, c_xtip, c_ytip, c_ztip, c_ctip, 0.0, c_n; r=1/c_r)
    vlm.setcoordsystem(canard, c_O, c_Oaxis)
    ## System
    system = vlm.WingSystem()
    vlm.addwing(system, "Wing", wing)
    vlm.addwing(system, "Canard", canard)
    if body; vlm.addwing(system, "Body", _body); end;

    # Solves
    vlm.solve(system, Vinf)
    vlm.calculate_field(system, "CFtot"; S=w_S+c_S)
    vlm.calculate_field(canard, "CFtot"; S=c_S)
    vlm.save(system, file_name; save_horseshoes=true, num=i)

    # Fluid domain
    if save_fdom
      P_max = [body ? b_l*2 : b_l*5/4, w_ytip*5/2, w_ytip*3/4]
      fdom = vlm.PP.FluidDomain(
              [0.0, 0.0, 0.0],                        # P_min
              P_max, # P_max
              [10, 5, 1]*2^3,                         # NDVIS
                        )
      vlm.PP.setcoordsystem(fdom,
              P_max.*[-2/15, -1/2, -1/3] + [body ? 0 : c_O[1], 0, 0 ],
              [1.0 0 0; 0 1 0; 0 0 1])
      V(X) = vlm.Vind(system, X) + Vinf(X,0)
      vlm.PP.calculate(fdom,
                [Dict("field_name"=>"U",
                      "field_type"=>"vector",
                      "field_function"=>V)]
                )
      vlm.PP.save(fdom, file_name; num=i)
    end

    c_info = vlm.fields_summary(canard)
    m_info = vlm.fields_summary(system)
    push!(c_CLs, c_info["CL"]);
    push!(m_CLs, m_info["CL"]);
  end

  # # Weber's experimental data (Table 4)
  # data = [[2.1, 4.2, 6.3, 8.4, 10.5],
  #           [0.121, 0.238, 0.350, 0.456, 0.559],
  #           [nothing, 0.005, 0.012, 0.022, 0.035]] # [alpha, Cl, Cd]

  # --------- PLOTS
  fig = figure("canard_wing",figsize=(7*2,7*1))
  suptitle("Effect of canard on wing Lambda=$(round(w_sweep[w_i]*180/pi,1))"*
        " for z=$(round(pivot_z[pivot_z_i],1)) and x=$(round(pivot_pos[pivot_pos_i],1))"*
        ", at Re=$(Re[Re_i])", fontsize="x-large")

  # _c_CLs = [ c_CLs[i]*(alphas[i]<0 ? -1:1) for i in 1:length(alphas)]
  _c_CLs = c_CLs
  _m_CLs = [ m_CLs[i]*(alphas[i]<0 ? -1:1) for i in 1:length(alphas)]

  # CANARD BALANCE
  subplot(121)
  # plot(data[1], data[2], "ok",
  #       label="Weber's experimental data")
  plot(alphas*180/pi, _c_CLs, "or", label="FLOWVLM")
  xlim([-4,24])
  xticks(-4:4:24)
  xlabel("Angle of attack (deg)")
  ylim([-0.3, 0.7])
  yticks(-0.3:0.1:0.7)
  ylabel("CL")
  grid(true, color="0.8", linestyle="--")
  legend(loc="best")
  title("Canard balance")

  # MAIN BALANCE
  subplot(122)
  # plot(data[1], data[3], "ok",
  #       label="Weber's experimental data")
  plot(alphas*180/pi, _m_CLs, "or", label="FLOWVLM")
  xlim([-4,24])
  xticks(-4:4:24)
  ylim([-0.6, 1.4])
  yticks(-0.6:0.2:1.4)
  xlabel("Angle of attack (deg)")
  ylabel("CL")
  grid(true, color="0.8", linestyle="--")
  legend(loc="best")
  title("Main balance")
end

"""
Validation of wing and tail interaction based on experimental data in Buell's
1957 *The effects at subsonic speeds of wing fences and a tail on the
longitudinal characteristics of a 63deg swept-wing and fuselage
combination*.
"""
function wing_tail_interaction(; n=1, save_w=false, save_fdom=false)
  run_name = "wtinter"
  rhoinf = 1.177            # kg/m^3, density of air at 1 atm and 300K
  muinf = 1.846/10^5        # kg/ms, dynamic viscosity of air at 1 atm and 300K

  # Experimental setup
  Re = 7*10^6             # Reynolds number
  Ma = 0.20               # Mach number
  w_barc = 1.20*0.3048    # (m) mean chord
  magVinf = Re*muinf/(rhoinf*w_barc)
  alphas = [i for i in -4:1:22]*pi/180  # Angle of attack

  # Model geometry
  ## Wing
  w_lambda = 63.0             # Leading edge sweep
  w_ar = 3.50                 # Aspect ratio
  w_tr = 0.25                 # Taper ratio
  w_b = 3.75*0.3048           # (m) span
  w_S = 4.02*0.3048^2         # (m^2) planform area
  ## Vertical tail
  vt_lambda = 54
  vt_ar = 1.51
  vt_tr = 0.16
  vt_b = 1.27*0.3048
  vt_S = 1.07*0.3048^2
  ## Horizontal tail
  ht_alphas = vcat(["off"], [0.2, -3.9, -7.8, -11.7]*pi/180)
  ht_lambdas = [0, 60.0]
  ht_ars = [4.0, 2.5]
  ht_trs = [0.33, 0.20]
  ht_bs = [1.87, 1.58]*0.3048
  ht_Ss = [0.87, 1.0]*0.3048^2
  ht_pivots = [0.45, 0.84]    # Pivot axis, fraction of root chord
  ## Fuselage
  b_l = (72+10.5)*0.0254      # (m) long fuselage length
  b_r = 0.13*w_b/2            # Radius
  b_S = 0.13*0.3048^2
  ## Fence
  f_poss = [.29, .45, .70]    # X fence positions, factor of semi-span

  # Simulation parameters
  ht_alphas_to_calculate = [1, 4]
  ht_i = 2                    # Type of horizontal tail
  n_w = n*50                    # Number of lattices
  n_vt = n*25
  n_ht = n*25
  n_b = n*5
  n_f = n*3
  r_w = 10.0                  # Lattice expansion ratio
  r_vt = 1/10.0
  r_ht = 10.0
  r_b = 1.0
  r_f = 1.0

  # Simulation geometry
  ## Wing
  w_O = [21.92, 0,0 ]*0.0254  # (m) position
  w_Oaxis = [1.0 0 0; 0 1.0 0; 0 0 1.0]   # Orientation
  w_twist_root = -0.4
  w_twist_tip = -3.5
  wing = vlm.simpleWing(w_b, w_ar*0 + w_b/(2*w_barc/(1/w_tr+1)), w_tr, w_twist_root, w_lambda, 0.0;
                      twist_tip=w_twist_tip, n=n_w, r=r_w)
  vlm.setcoordsystem(wing, w_O, w_Oaxis)
  ## Vertical tail
  vt_barc = vt_S/vt_b
  vt_ar = vt_b/(2*vt_barc/(1/vt_tr+1))
  vt_xroot = 0.0
  vt_yroot = 0.0
  vt_zroot = 0.0
  vt_ytip = vt_b
  vt_xtip = vt_ytip*tan(vt_lambda*pi/180)
  vt_ztip = 0.0
  vt_ctip = vt_b/vt_ar
  vt_croot = vt_ctip/vt_tr
  vt_O = [b_l-vt_croot, 0, b_r]   # Places it flush with the fuselage
  vt_Oaxis = [1.0 0 0; 0 0 1.0; 0 -1.0 0] # Places it vertically
  vertical_tail = vlm.Wing(vt_xroot, vt_yroot, vt_zroot, vt_croot, 0.0)
  vlm.addchord(vertical_tail, vt_xtip, vt_ytip, vt_ztip, vt_ctip, 0.0,
                n_vt; r=r_vt)
  vlm.setcoordsystem(vertical_tail, vt_O, vt_Oaxis)
  ## Fuselage
  b_lambda = 70*pi/180
  b_O = [0.0,0,0]
  b_Oaxis = [1.0 0 0; 0 0 1.0; 0 -1.0 0]
  body = vlm.Wing(b_r*tan(b_lambda), -b_r, 0.0, b_l-b_r*tan(b_lambda), 0.0)
  vlm.addchord(body, 0.0, 0.0, 0.0, b_l, 0.0, n_b; r=r_b)
  vlm.addchord(body, b_r*tan(b_lambda), b_r, 0.0, b_l-b_r*tan(b_lambda), 0.0,
                n_b; r=1/r_b)
  vlm.setcoordsystem(body, b_O, b_Oaxis)

  # Data storage
  CLs, CDs, CMs  = Dict(), Dict(), Dict()
  for ht_alpha in [ht_alphas[i] for i in ht_alphas_to_calculate]
    CLs[ht_alpha] = []
    CDs[ht_alpha] = []
    CMs[ht_alpha] = []
  end

  # Iterates over angles of horizontal tail incidence
  for (i,ht_alpha) in enumerate([ht_alphas[k] for k in ht_alphas_to_calculate])
    println("Case tail=$ht_alpha")
    pivot_x = (67.32+10.50)*0.0254  # (m) position of pivot from fuselage nose
    if ht_alpha!="off"
      # Builds horizontal tail
      ht_barc = ht_Ss[ht_i]/ht_bs[ht_i]
      ht_ctip = 2*ht_barc/(1/ht_trs[ht_i]+1)
      ht_ar = ht_bs[ht_i]/ht_ctip
      horizontal_tail = vlm.simpleWing(ht_bs[ht_i],
                            ht_ars[ht_i]*0 + ht_ar, ht_trs[ht_i],
                            0.0, ht_lambdas[ht_i], 0.0; twist_tip=0.0,
                            n=n_ht, r=r_ht)
      # Sets the angle of incidence
      ht_croot = ht_bs[ht_i]/ht_ar/ht_trs[ht_i]  # Root chord
      ht_pivot = ht_pivots[ht_i]      # Pivot axis, fraction of root chord
      pivotroot = ht_pivot*ht_croot   # Distance of pivot point from canard nose
      nose_x = pivot_x - pivotroot*cos(ht_alpha)    # x-position of canard nose
      nose_z = pivotroot*sin(ht_alpha)              # z-position of canard nose
      ht_O = [nose_x, 0, nose_z]      # position of canard nose
      ht_Oaxis = [cos(ht_alpha) 0 -sin(ht_alpha); 0 1.0 0;
                  sin(ht_alpha) 0 cos(ht_alpha)]   # Orientation of canard
      vlm.setcoordsystem(horizontal_tail, ht_O, ht_Oaxis)
    end

    # Builds the system
    system = vlm.WingSystem()
    vlm.addwing(system, "Wing", wing)
    vlm.addwing(system, "Fuselage", body)
    vlm.addwing(system, "VerticalTail", vertical_tail)
    if ht_alpha!="off"
      vlm.addwing(system, "HorizontalTail", horizontal_tail)
    end

    file_name = run_name*"_$(i)_"

    # Iterates over angle of attack
    for (j,alpha) in enumerate(alphas)
        println("\tAOA: $alpha")
        Vinf(X,t) = magVinf*[cos(alpha), 0, sin(alpha)]
        # Solves
        vlm.solve(system, Vinf)
        vlm.calculate_field(system, "Ftot"; rhoinf=rhoinf)
        vlm.calculate_field(system, "CFtot"; S=w_S+ht_Ss[ht_i])
        vlm.calculate_field(system, "Mtot"; r_cg=[pivot_x-2.16*0.5*0.3048, 0, 0])
        vlm.calculate_field(system, "CMtot"; qinf=(1/2)*rhoinf*magVinf^2,
                              S=w_S+ht_Ss[ht_i], l=0.5*0.3048)
        if save_w
          vlm.save(system, file_name; save_horseshoes=true, num=j)
        end

        # Fluid domain
        if save_fdom
          println("Calculating FDOM")
          P_max = [b_l*5/4, w_b/2*5/2, vt_b*2]
          fdom = vlm.PP.FluidDomain(
                  [0.0, 0.0, 0.0],                        # P_min
                  P_max, # P_max
                  [10, 5, 1]*2^3,                         # NDVIS
                            )
          vlm.PP.setcoordsystem(fdom,
                  P_max.*[-2/15*0, -1/2, -1/4],
                  [1.0 0 0; 0 1 0; 0 0 1])
          V(X) = vlm.Vind(system, X) + Vinf(X,0)
          vlm.PP.calculate(fdom,
                    [Dict("field_name"=>"U",
                          "field_type"=>"vector",
                          "field_function"=>V)]
                    )
          vlm.PP.save(fdom, file_name; num=j)
        end

        info = vlm.fields_summary(system)
        push!(CLs[ht_alpha], info["CL"]);
        push!(CDs[ht_alpha], info["CD"]);
        push!(CMs[ht_alpha], info["CMtot"]);
      end
  end
  save("$run_name.jld", "alphas", alphas, "CLs", CLs, "CDs", CDs, "CMs", CMs)
  wing_tail_interaction_plot("$run_name.jld", "Wing-tail interaction"*
                              " at Ma=$Ma, Re=$Re")
end
function wing_tail_interaction_plot(jld_file_name, dscrptn)
  alphas, CLs, CDs, CMs = load("wtinter.jld", "alphas", "CLs", "CDs", "CMs")


  # Experimental data
  data_alphas = [22, 18, 13, 9.5, 5, 0.5, -4]
  data_CLs = Dict("off" => [.96, 0.8, .6, .4, .2, 0, -.2],
                  "7.8" => [1.06, 0.88, .64, .4, .2, -.08, -.3])

  # ------- PLOTS
  fig = figure("validation_wing_tail",figsize=(7*3,5*1))
  suptitle(dscrptn, fontsize="x-large")


  # LIFT
  subplot(131)
  for key in keys(data_CLs)
    plot(data_alphas, data_CLs[key], "--.", label="Experimental tail=$key")
  end
  for key in keys(CLs)
    deg = key!="off"? round(key*180/pi,1) : key
    _CLs = [CLs[key][i] * (alphas[i]<0 ? -1 : 1) for i in 1:size(alphas)[1]]
    plot(alphas*180/pi, _CLs, "o", label="FLOWVLM tail=$deg")
  end
  xlim([-5,25]);
  ylim([-0.6, 1.4]);
  xlabel("Angle of attack (deg)");
  ylabel(L"$C_L$");
  grid(true, color="0.8", linestyle="--")
  legend(loc="best");
  title("Lift");

  # DRAG
  subplot(132)
  for key in keys(CLs)
    deg = key!="off"? round(key*180/pi,1) : key
    # plot(data[1], data[2], "ok",
    #       label="Weber's experimental data")
    plot(CLs[key]./CDs[key], CDs[key], "o", label="FLOWVLM tail=$deg")
  end
  # xlim([0, 16]);
  # ylim([-0.6, 1.4]);
  xlabel(L"\frac{L}{D}");
  ylabel(L"$C_L$");
  grid(true, color="0.8", linestyle="--")
  legend(loc="best");
  title("Lift-Drag");

  # DRAG
  subplot(133)
  for key in keys(CLs)
    deg = key!="off"? round(key*180/pi,1) : key
    # plot(data[1], data[2], "ok",
    #       label="Weber's experimental data")
    plot(CMs[key], CLs[key], "o", label="FLOWVLM tail=$deg")
  end
  # xlim([-0.1, 0.05]);
  # ylim([0.05, 1.4]);
  xlabel(L"C_m");
  ylabel(L"$C_L$");
  grid(true, color="0.8", linestyle="--")
  legend(loc="best");
  title("Pitching moment");
end
##### END OF VALIDATION ########################################################

################################################################################
# VERIFICATION CASES
################################################################################
"""
Comparison with theoretical calculation on two planar wings of sweep 0 and 35.
For details, see Validation 2 in *Great OWL Publishing - Surfaces; Vortex
Lattice Module*, pp. 106. It also performs a grid dependance study to determine
dependance of aerodynamic characteristics with lattice refinement.
"""
function ver_planarWing(; save_w=false, drag_trick=false)

  save = save_w  # Saves wings
  save_fd = false # Saves fluid domains

  # Theoretical results
  theoretical = Dict( 0.0 => Dict(
                                  "CLalpha" => 0.08846, #per deg
                                  "CL"      => 0.8846,
                                  "L"       => 299.7*4.44822, # N
                                  "CD"   => 0.02491,
                                  "D"    => 8.4*4.44822 # N
                                  ),
                      35.0 => Dict(
                                  "CLalpha" => 0.07480, #per deg
                                  "CL"      => 0.7480,
                                  "L"       => 253.4*4.44822, # N
                                  "CD"   => 0.01781,
                                  "D"    => 6.0*4.44822 # N
                                  )
                      )

  to_compare = ["CD", "CL", "L", "D"]

  lambdas = [0.0, 35.0]
  b = 10*0.3048 # m
  c = 1*0.3048 # m
  alpha = 10.0 # deg
  gamma = 0.0
  S = 10*0.092903
  AR = b^2/S
  tr = 1.0
  n = 100

  magVinf = 168.8*0.3048 # m/s
  Vinf(X, t) = magVinf*[ cos(alpha*pi/180), 0.0, sin(alpha*pi/180)]
  rhoinf = 0.002378*515.379 # kg/m^3

  for lambda in lambdas
    println("###### LAMBDA $lambda #######")
    function wing_fun(n)
      return vlm.simpleWing(b, AR, tr, 0.0, lambda, gamma; n=n, r=10.0)
    end

    println("------------ GRID CONVERGENCE")
    wing, ns, vals = vlm.grid_dependance(wing_fun, Vinf; ns=[4*2^i for i in 1:6],
                              S=S, rhoinf=rhoinf,
                              verbose=true, save_w=save,
                              save_name="planar$(Int(lambda))", save_fd=save_fd,
                              fig_title="Planar wing - Sweep $lambda")
    println("------------ GRID CONVERGENCE DONE")

    vlm.calculate_field(wing, "L"; rhoinf=rhoinf)
    info = vlm.fields_summary(wing; drag_trick=drag_trick)
    println("\n\tPARAM\tTHEORETICAL\tFLOWVLM\t\tError")
    for key in to_compare
      theo = theoretical[lambda][key]
      myvlm = info[key]
      err = (theo-myvlm)/myvlm
      sigfigs = 3
      ordr = log10(theo)
      rnd = sigfigs - ordr <= 0 ? 0 : Int(round(sigfigs-ordr))
      println("\t>$(key)\t$(round(theo,rnd))\t\t$(round(myvlm,rnd))\t\t$(round(err*100,1))%")
    end
    println("###### Done\n")
  end
end

"Aerodynamic characteristics on Bertin's wing when the AOA is included
in the incident freestream compared to when the AOA is included in the local
coordinate system"
function ver_local_coord_sys(; save_w=false)

  function print_dict(dict)
    for key in keys(dict)
      println("\t$key:\t$(dict[key])")
    end
  end

  magVinf = 163*0.3048 # m/s
  rhoinf = 9.093/10^1

  alpha = 10.5
  b=98*0.0254
  lambda = 45.0
  ar = 5.0
  tr = 1.0
  gamma = 0.0
  twist = 0.0
  n = 4*2^4
  r = 1.0
  central = false

  # ------------ AOA IN INCIDENT FREESTREAM
  Vinf1(X, t) = magVinf*[ cos(alpha*pi/180), 0.0, sin(alpha*pi/180)]
  wing1 = vlm.simpleWing(b, ar, tr, twist, lambda, gamma;
                        n=n, r=r, central=central)
  vlm.solve(wing1, Vinf1)
  vlm.calculate_field(wing1, "Ftot"; rhoinf=rhoinf)
  vlm.calculate_field(wing1, "Mtot")

  info1 = vlm.fields_summary(wing1)
  if save_w
    vlm.save(wing1, "ver_coord_sys_wing1"; save_horseshoes=true)
  end

  println("######## AOA IN INCIDENT FREESTREAM ########")
  print_dict(info1)

  # ------------ AOA IN LOCAL COORDINATE SYSTEM
  # (It aligns the axis of the wing with the z-global-axis)
  function Vinf(X, t)
    return magVinf*[ 0, 0.0, 1]
  end

  wing2 = vlm.simpleWing(b, ar, tr, twist, lambda, gamma;
                          n=n, r=r, central=central)

  # Local coordinate system
  kp = -[ cos(alpha*pi/180), 0.0, -sin(alpha*pi/180)]
  ip = [ sin(alpha*pi/180), 0.0, cos(alpha*pi/180)]
  jp = [0.0,1.0,0.0]
  T = [b/1,b/2,-b/3]
  vlm.setcoordsystem(wing2, T, [ip,jp,kp])

  vlm.solve(wing2, Vinf)
  vlm.calculate_field(wing2, "Ftot"; rhoinf=rhoinf)
  vlm.calculate_field(wing2, "Mtot")

  info2 = vlm.fields_summary(wing2)
  if save_w
    vlm.save(wing2, "ver_coord_sys_wing2"; save_horseshoes=true)
  end

  println("######## AOA IN LOCAL COORDINATE SYSTEM ########")
  print_dict(info2)

  return wing1, wing2
end
##### END OF VERIFICATION ######################################################
