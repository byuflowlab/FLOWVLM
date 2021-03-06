import Dierckx

include("../src/FLOWVLM.jl")
vlm = FLOWVLM


# ------------ GLOBAL VARIABLES ------------------------------------------------
global module_path; module_path,_ = splitdir(@__FILE__);   # Path to this module
global data_path = joinpath(module_path,"data/")            # Path to data
global save_horseshoes = true

function distributedpropulsion(; save_path=module_path*"/../temps/distprop04/",
                                  paraview=true, verbose=false, delete=false,
                                  prompt=false,
                                  rotor_file="apc9x4_5.jl",
                                  runccblade=false)

  save_name = "distprop"

  println("Defining parameters...")
  ##############################################################################
  # PARAMETERS
  ##############################################################################
  # ------------ OPERATION PARAMETERS ------------------------------------
  # Propellers
  np_w = 2                      # Number of propellers on wing on each side
  np_tw = 2                     # Number of propellers on tandem wing each side
  RPM = 3000                    # (rev/min) RPM
  omega = 2*pi*RPM/60           # (rad/s) angular velocity
  pitch = 0.0                   # (deg) pitch of propellers
  CW = false                    # Rotation orientation of first propeller
  n_ccb = 8                     # Number of CCBlade elements

  # Environment
  magVinf = 21.0                 # (m/s) inflow
  Vinf(X,t) = magVinf*[1.0,0,0]
  rho = 1.225                   # (kg/m^3) air density
  Re = 5*10^5                   # This is especulation

  # ------------ GEOMETRIC PARAMETERS ------------------------------------
  # Wing
  b_w = 6.0                     # (m) span
  AR_w = 6.0                    # Aspect ratio
  tr_w = 1.0                    # Taper ratio
  twist_r_w = 7.5               # (deg) twist at root
  twist_t_w = 0.0               # (deg) twist at tip
  lambda_w = 0.0                # (deg) sweep
  gamma_w = 0.0                 # (deg) dihedral
  n_w = 12                      # Number of horseshoes per side of wing
  r_w = 2.0                     # Horseshoe expansion ratio

  # Wing winglets
  b_wl = b_w/4                  # (m) span of winglet from top to bottom
  AR_wl = 3.0                   # Aspect ratio
  tr_wl = (b_wl/AR_wl)/(b_w/AR_w/tr_w)      # Taper ratio
  twist_r_wl = 2.5              # (deg) twist at root
  twist_t_wl = 0.0              # (deg) twist at tip
  lambda_wl = 40.0              # (deg) sweep
  gamma_wl = 15.0               # (deg) dihedral
  n_wl = 12                     # Number of horseshoes per side of wing
  r_wl = 2.0                    # Horseshoe expansion ratio

  # Tandem wing
  b_tw = b_w*0.9                # (m) span
  AR_tw = 10.0                  # Aspect ratio
  tr_tw = 1.0                   # Taper ratio
  twist_r_tw = 5.0              # (deg) twist at root
  twist_t_tw = 0.0              # (deg) twist at tip
  lambda_tw = 0.0               # (deg) sweep
  gamma_tw = 0.0                # (deg) dihedral
  n_tw = n_w                    # Number of horseshoes per side of wing
  r_tw = r_w                    # Horseshoe expansion ratio

  # Fuselage
  l_f = 5.8                     # (m) length
  h_f = 2.75/2                  # (m) height
  c_pos_f1 = 0                  # (m) span position of first chord
  c_pos_f2 = 0.2*h_f            # (m) span position of second chord
  c_pos_f3 = 0.7*h_f
  c_pos_f4 = h_f
  c_f1 = 0.66*l_f               # (m) length of first chord
  c_f2 = 0.75*l_f
  c_f3 = 0.75*l_f
  c_f4 = 0.50*l_f
  x_pos_f1 = 0.025*l_f          # (m) x-position of first chord
  x_pos_f2 = 0.0*l_f
  x_pos_f3 = 0.20*l_f
  x_pos_f4 = 0.50*l_f


  # ------------ ASSEMBLY PARAMETERS ------------------------------------
  # Position of wings on fuselage
  h_pos_w = 0.90*h_f            # (m) height position of wing
  h_pos_tw = 0.15*h_f           # (m) height position of tandem wing
  l_pos_w = 0.95*l_f-b_w/AR_w/tr_w   # (m) length position of wing
  l_pos_tw = 0.05*l_f           # (m) length position of tandem wing

  # Position of propellers on main wing
  d_prop_w = b_w/2/np_w         # Distance between propellers on wing
  y_pos_prop_w = Float64[b_w/2 - i*d_prop_w for i in 0:np_w-1] # y-positions

  # Position of propellers on tandem wing
  d_prop_tw = b_tw/2/np_tw         # Distance between propellers on wing
  y_pos_prop_tw = Float64[b_tw/2 - i*d_prop_tw for i in 0:np_tw-1] # y-positions


  # ------------ SIMULATION PARAMETERS -----------------------------------
  # First stage: VTOL
  # Second stage: transition frontal flight
  # Third stage: frontal flight
  # Fourth stage: transition to VTL
  nstages = 4
  nrevs = [1, 3, 1, 1]            # Number of revolutions on each stage
  totnrevs = sum(nrevs)              # Total number of revolutons
  init_ori = 90.0                 # Initial orientation of wings
  rot_w = Float64[-45, -45, 0, 90]# (deg) main wing rotation on each stage
  rot_tw = Float64[0, -90, 0, 90] # (deg) tandem wing rotation on each stage
  steps_per_rev = 360/10          # Steps in each revolution
  nsteps = Int(ceil(steps_per_rev*totnrevs))       # Total number of steps
  dangle_prop = totnrevs*360/nsteps  # (deg) propeller rotation on each step
  dangle_w = rot_w./(nrevs*steps_per_rev)  # (deg) wing rotation on each step per stage
  dangle_tw = rot_tw./(nrevs*steps_per_rev)# (deg) wing rotation on each step per stage


  println("Generating geometry...")
  ##############################################################################
  # ASSEMBLY
  ##############################################################################
  # ------------ MAIN WING ---------------------------------------------
  # Generates base propellers (one on each rotation orientation)
  propellers = vlm.Rotor[]
  println("\tGenerating first propeller...")
  @time push!(propellers, generate_rotor(pitch; n=n_ccb, CW=!CW, Re=Re,
                          verbose=verbose, xfoil=false, rotor_file=rotor_file))
  println("\tGenerating second propeller...")
  @time push!(propellers, generate_rotor(pitch; n=n_ccb, CW=CW, Re=Re,
                          verbose=verbose, xfoil=false, rotor_file=rotor_file))

  # Generates wing
  println("\tGenerating main wing assembly...")
  wing = vlm.simpleWing(b_w, AR_w, tr_w, twist_r_w, lambda_w, gamma_w;
                                              twist_tip=twist_t_w, n=n_w, r=r_w)
  # Generates winglets
  winglet_R = vlm.simpleWing(b_wl, AR_wl, tr_wl, twist_r_wl, lambda_wl, gamma_wl;
                                          twist_tip=twist_t_wl, n=n_wl, r=r_wl)
  winglet_L = vlm.simpleWing(b_wl, AR_wl, tr_wl, twist_r_wl, lambda_wl, gamma_wl;
                                          twist_tip=twist_t_wl, n=n_wl, r=r_wl)
  O_wl_R = (b_w/2)*[tan(lambda_w*pi/180), 1, tan(gamma_w*pi/180)]
  O_wl_L = [1 0 0; 0 -1 0; 0 0 1]*O_wl_R
  Oaxis_wl_R = vlm.vtk.rotation_matrix(0.0, 0.0, 90.0)
  Oaxis_wl_L = vlm.vtk.rotation_matrix(0.0, 0.0, -90.0)
  vlm.setcoordsystem(winglet_R, O_wl_R, Oaxis_wl_R)
  vlm.setcoordsystem(winglet_L, O_wl_L, Oaxis_wl_L)

  ## Generates propellers on wing (from right to left)
  println("\t\tGenerating main wing propellers...")
  O_prop_w = [ ypos*[tan(lambda_w*pi/180), 1, tan(gamma_w*pi/180)]
                                              for ypos in y_pos_prop_w]
  props_w = vlm.Rotor[]
  for i in 1:2*np_w
    right = i<=np_w    # Indicates wich side of the wing
    this_prop = deepcopy(propellers[1+i%2]) # Alternates rotation orientation
    this_O = O_prop_w[ right ? i : np_w-(i-np_w-1)] # Chooses position
    this_O = [1 0 0; 0 (-1)^!right 0; 0 0 1]*this_O   # Places it in correct side

    vlm.setcoordsystem(this_prop, this_O, Float64[1 0 0; 0 1 0; 0 0 1]; user=true)
    push!(props_w, this_prop)
  end

  # Assembles main wing
  main_wing = vlm.WingSystem()
  vlm.addwing(main_wing, "Wing", wing)
  vlm.addwing(main_wing, "WingletR", winglet_R)
  vlm.addwing(main_wing, "WingletL", winglet_L)
  for (i, prop) in enumerate(props_w)
    vlm.addwing(main_wing, "Prop$i", prop)
  end

  # Position of main wing
  O_w = [l_pos_w, 0, h_pos_w]
  Oaxis_w = vlm.vtk.rotation_matrix(0.0, -init_ori, 0.0)
  vlm.setcoordsystem(main_wing, O_w, Oaxis_w)

  # ------------ TANDEM WING -------------------------------------------
  println("\tGenerating tandem wing assembly...")
  # Generates tandem wing
  twing = vlm.simpleWing(b_tw, AR_tw, tr_tw, twist_r_tw, lambda_tw,
                                gamma_tw; twist_tip=twist_t_tw, n=n_tw, r=r_tw)

  ## Generates propellers on tandem wing (from right to left)
  println("\t\tGenerating tandem wing propellers...")
  O_prop_tw = [ ypos*[tan(lambda_tw*pi/180), 1, tan(gamma_tw*pi/180)]
                                              for ypos in y_pos_prop_tw]
  props_tw = vlm.Rotor[]
  for i in 1:2*np_tw
    right = i<=np_tw    # Indicates wich side of the wing
    this_prop = deepcopy(propellers[1+i%2]) # Alternates rotation orientation
    this_O = O_prop_tw[ right ? i : np_tw-(i-np_tw-1)] # Chooses position
    this_O = [1 0 0; 0 (-1)^!right 0; 0 0 1]*this_O   # Places it in correct side

    vlm.setcoordsystem(this_prop, this_O, Float64[1 0 0; 0 1 0; 0 0 1]; user=true)
    push!(props_tw, this_prop)
  end

  # Assembles tandem wing
  tandem_wing = vlm.WingSystem()
  vlm.addwing(tandem_wing, "Wing", twing)
  for (i, prop) in enumerate(props_tw)
    vlm.addwing(tandem_wing, "Prop$i", prop)
  end

  # Position of tandem wing
  O_tw = [l_pos_tw, 0, h_pos_tw]
  Oaxis_tw = vlm.vtk.rotation_matrix(0.0, -init_ori, 0.0)
  vlm.setcoordsystem(tandem_wing, O_tw, Oaxis_tw)

  # ------------ FUSELAGE ----------------------------------------------
  # Generates fuselage
  fuselage = vlm.Wing(x_pos_f1, c_pos_f1, 0.0, c_f1, 0.0)
  vlm.addchord(fuselage, x_pos_f2, c_pos_f2, 0.0, c_f2, 0.0, 1)
  vlm.addchord(fuselage, x_pos_f3, c_pos_f3, 0.0, c_f3, 0.0, 1)
  vlm.addchord(fuselage, x_pos_f4, c_pos_f4, 0.0, c_f4, 0.0, 1)
  Oaxis_f = vlm.vtk.rotation_matrix(0.0, 0.0, -90.0)
  vlm.setcoordsystem(fuselage, zeros(Float64,3), Oaxis_f)


  # ------------ SYSTEM ------------------------------------------------
  # Creates system assembly
  system = vlm.WingSystem()
  vlm.addwing(system, "MainWing", main_wing)
  vlm.addwing(system, "TandemWing", tandem_wing)
  vlm.addwing(system, "Fuselage", fuselage)




  println("Starting simulation...")
  ##############################################################################
  # SIMULATION
  ##############################################################################
  # ------------ SIMULATION SETUP ----------------------------------------
  vlm.setVinf(system, Vinf)
  for prop in props_w; vlm.setRPM(prop, RPM); end;
  for prop in props_tw; vlm.setRPM(prop, RPM); end;

  # Creates save path
  if save_path!=nothing; vlm.create_path(save_path, prompt); end;

  # ------------ RUN SIMULATION ------------------------------------------
  for i in 0:nsteps-1
    if i%10==0; println("\tStep $i out of $nsteps"); end;
    cur_rev = floor(i/steps_per_rev) # Revolutions completed

    # Rotates the propellers
    rotation = i==0 ? 0.0 : dangle_prop
    for prop in props_w; vlm.rotate(prop, rotation); end;
    for prop in props_tw; vlm.rotate(prop, rotation); end;

    # Identifies operation stage
    stg_i = 0
    for j in 1:nstages
      if cur_rev<sum(nrevs[1:j])
        stg_i = j
        break
      end
    end

    # Tilts main wing
    rot_Oaxis_w = vlm.vtk.rotation_matrix(0.0, -dangle_w[stg_i], 0.0)
    new_Oaxis_w = rot_Oaxis_w*main_wing.Oaxis
    vlm.setcoordsystem(main_wing, main_wing.O, new_Oaxis_w)

    # Tilts tandem wing
    rot_Oaxis_tw = vlm.vtk.rotation_matrix(0.0, -dangle_tw[stg_i], 0.0)
    new_Oaxis_tw = rot_Oaxis_tw*tandem_wing.Oaxis
    vlm.setcoordsystem(tandem_wing, tandem_wing.O, new_Oaxis_tw)

    # # Calculates distributed loads from CCBlade
    # if runccblade
    #   vlm.calc_distributedloads(rotor1, Vinf, RPM, rho; include_comps=true)
    #   vlm.calc_distributedloads(rotor2, Vinf, RPM, rho; include_comps=true)
    #   vlm.calc_distributedloads(rotor3, Vinf, RPM, rho; include_comps=true)
    #   vlm.calc_distributedloads(rotor4, Vinf, RPM, rho; include_comps=true)
    #   vlm.calc_distributedloads(rotor5, Vinf, RPM, rho; include_comps=true)
    # end

    # Saves
    if save_path!=nothing
      vlm.save(system, save_name; save_horseshoes=save_horseshoes,
                              path=save_path, num=i)
      # # Overwrites the wing so to save horseshoes
      # vlm.save(wing, save_name*"_Wing"; save_horseshoes=true,
      #                         path=save_path, num=i)
    end
  end


  # ------------ VISUALIZATION -------------------------------------------
  # Call Paraview and deletes temp files
  if save_path!=nothing

    if paraview
      strn = ""
      strn = strn * save_name * "_MainWing_Wing_vlm...vtk;"
      strn = strn * save_name * "_MainWing_WingletR_vlm...vtk;"
      strn = strn * save_name * "_MainWing_WingletL_vlm...vtk;"
      strn = strn * save_name * "_TandemWing_Wing_vlm...vtk;"
      strn = strn * save_name * "_Fuselage_vlm...vtk;"
      num_blades = propellers[1].B
      for (sys_name, num_props) in [("MainWing", np_w*2), ("TandemWing", np_tw*2)]
        for j in 1:num_props
          for i in 1:num_blades
            strn = strn * save_name * "_" * sys_name * "_Prop$j" * "_Blade$(i)_vlm...vtk;"
          end
          for i in 1:num_blades
            strn = strn * save_name * "_" * sys_name * "_Prop$j" * "_Blade$(i)_loft...vtk;"
          end
        end
      end
      # println("Calling paraview on $save_path$strn")
      run(`paraview --data="$save_path$strn"`)
    end

    if delete
      run(`rm -rf $save_path`)
    end
  end

  return system
end



"Generates the Rotor object. `pitch` is the pitch of the blades in degrees,
`n` is the number of lattices in the VLM."
function generate_rotor(pitch::Float64; n=10, CW=true, Re=5*10^5,
                          verbose=false, xfoil=false,
                          rotor_file="apc9x4_5.jl")
  n_bem = n
  # include(data_path*"caradonna_rotor.jl")
  # include(data_path*"apc9x4_5.jl")
  include(data_path*rotor_file)

  # Splines
  spline_k = 3
  spline_bc = "extrapolate"
  spline_s = 0.001
  _spl_chord = Dierckx.Spline1D(chorddist[:, 1]*Rtip, chorddist[:, 2]*Rtip;
                                k= size(chorddist)[1]>2 ? spline_k : 1,
                                s=spline_s, bc=spline_bc)
  _spl_theta = Dierckx.Spline1D(pitchdist[:, 1]*Rtip, pitchdist[:, 2];
                                k= size(pitchdist)[1]>2 ? spline_k : 1,
                                s=spline_s, bc=spline_bc)
  _spl_LE_x = Dierckx.Spline1D(sweepdist[:, 1]*Rtip, sweepdist[:, 2]*Rtip;
                                k= size(sweepdist)[1]>2 ? spline_k : 1,
                                s=spline_s, bc=spline_bc)
  _spl_LE_z = Dierckx.Spline1D(heightdist[:, 1]*Rtip, heightdist[:, 2]*Rtip;
                                k= size(heightdist)[1]>2 ? spline_k : 1,
                                s=spline_s, bc=spline_bc)
  spl_chord(x) = Dierckx.evaluate(_spl_chord, x)
  spl_theta(x) = pitch + Dierckx.evaluate(_spl_theta, x)
  spl_LE_x(x) = Dierckx.evaluate(_spl_LE_x, x)
  spl_LE_z(x) = Dierckx.evaluate(_spl_LE_z, x)

  # Geometry for CCBlade & FLOWVLM
  r = [Rhub + i*(Rtip-Rhub)/n_bem for i in 0:n_bem] # r is discretized in n+1 sections
  chord = [spl_chord(ri) for ri in r]
  theta = [spl_theta(ri) for ri in r]
  LE_x = [spl_LE_x(ri) for ri in r]
  LE_z = [spl_LE_z(ri) for ri in r]


  airfoils = []
  for (pos, contour, file_name) in airfoil_contours
    x, y = contour[:,1], contour[:,2]
    # Calls XFOIL to calculate polars of each airfoil
    if xfoil
      polar = vlm.ap.runXFOIL(x, y, Re;
                        alphas=[i for i in -20:1.0:20],
                        verbose=verbose, Mach=0.0,
                        iter=100, alpha_ite=10)
    # Reads polars from files
    else
      println(file_name)
      polar = vlm.ap.read_polar(file_name; path=data_path*"airfoils/", x=x, y=y)
    end

    if verbose; vlm.ap.plot(polar; label="pos=$pos"); end;
    push!(airfoils, (pos, polar))
  end

  propeller = vlm.Rotor(CW, r, chord, theta, LE_x, LE_z, B, airfoils)
  # propeller = vlm.Rotor(CW, r, chord, theta, LE_x, LE_z, B)
  vlm.initialize(propeller, n)

  return propeller
end
