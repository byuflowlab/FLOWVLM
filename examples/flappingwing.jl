# Bald Eagly flapping around

include("../src/FLOWVLM.jl")
vlm = FLOWVLM

save_path = "temp_flapping/"
run_name = "flappingwing"
save_horseshoes = true
save_code = "examples/flappingwing.jl"
prompt = true

function flappingwing(; verbose=true, save_fdom=false)
  # ---------------- RUN PARAMETERS -------------------------
  flaps_per_sec = 6/2.6       # (1/s) Strikes per second
  num_strikes = 3             # Simulation strikes length
  nsteps = 20                 # Number of time steps
  A = 55                      # (deg) Strike amplitude
  offset = 35/55              # Strike offset
  amplif = 1.5                # Wingtip section amplification

  # ---------------- SIMULATION PARAMETERS ------------------
  omega = flaps_per_sec * 2*pi# (rad/s) Angular speed
  tend = num_strikes/flaps_per_sec      # Simulation end time
  dt = tend/nsteps            # Time step
  verbose_nsteps = 10         # Verbose frequency
  n = 1                       # Lattice density
  fdom_n = 1                  # Fluid domain density

  dihedrals(t) = forearm_wingtip_dihedral(t, omega, A, offset, amplif)

  # ---------------- PREPARATION ----------------------------
  # Creates save path
  if save_path!=nothing
    # Checks if folder already exists
    if isdir(save_path)
      if prompt
        inp1 = ""
        opts1 = ["y", "n"]
        while false==(inp1 in opts1)
          print("\n\nFolder $save_path already exists. Remove? (y/n) ")
          inp1 = readline()[1:end-1]
        end
        if inp1=="y"
          run(`rm $save_path -rf`)
          println("\n")
        else
          return
        end
      else
        run(`rm $save_path -rf`)
      end
    end
    run(`mkdir $save_path`)

    # Saves code
    if save_code!=""
      run(`cp $save_code $save_path`)
    end
  end


  # ---------------- RUN SIMULATION -------------------------
  # RUN
  if verbose
    time_beg = Dates.DateTime(now())
  println("*******************************************************************")
  println("START $(save_path!=nothing ? save_path*run_name : "")\t$(time_beg)")
  println("*******************************************************************")
  end
  for i in 0:nsteps
    # Verbose
    if verbose && i%verbose_nsteps==0
      println("Time step $i out of $nsteps")
    end

    # Time step
    t = dt*i
    forearm_gamma, wingtip_gamma = dihedrals(t)
    system, fdom = generate_eagle(; gamma2=forearm_gamma*pi/180,
                                  gamma3=wingtip_gamma*pi/180, num=i,
                                  save_path=nothing, save_fdom=false,
                                  generate_fdom=save_fdom,
                                  openParaview=false)

    # Saves vtks
    if save_path!=nothing
      vlm.save(system, run_name; num=i,
              save_horseshoes=save_horseshoes, path=save_path)

      # Saves fluid domain
      if save_fdom
        vlm.PP.save(fdom, run_name; path=save_path, num=i)
      end
    end
  end
  if verbose
    time_end = Dates.DateTime(now())
    hrs,mins,secs = timeformat(time_beg, time_end)
  println("*******************************************************************")
  println("END $(save_path!=nothing ? save_path*run_name : "")\t$(time_end)")
  println("*******************************************************************")
  println("ELAPSED TIME: $(hrs) hours $mins minutes $secs seconds")
  end
end

"t is the current time, omega is the angular speed, A is amplitude of
the full strike, and offset is the offset of the strike (value from
-1 to 1), amplif is the amplification factor of the wingtip section.
It returns (forearm, wingtip) dihedral of each section"
function forearm_wingtip_dihedral(t, omega, A, offset, amplif)
    forearm = A*( sin(omega*t + asin(-offset)) )/2 + A/2*offset
    wingtip = forearm * amplif
    return forearm, wingtip
end

function generate_eagle(; gamma2=4*pi/180, gamma3=15*pi/180, n=2, fdom_n=2,
                        num=nothing,
                        save_path=save_path, save_fdom=true,
                        openParaview=false, generate_fdom=true)

  # ---------------- RUN PARAMETERS -------------------------
  Vinf_mag = 20.0   # (m/s) Flight speed
  AOA = 7.5*pi/180  # Angle of attack

  # n = 2             # Lattice density factor
  # fdom_n = 2        # Fluid domain density factor

  Vinf_vec = Vinf_mag*[1,0,0]
  Vinf(X,t) = Vinf_vec


  # ---------------- DIMENSIONS -----------------------------
  b = 2.3           # (m) Wingspan
  # mean_c = 0.7      # (m) Mean chord
  mean_c = 0.5      # (m) Mean chord

  # Fraction of the span
  b_sec1 = 1/9       # Body
  b_sec2 = 5/9       # Forearm
  b_sec3 = 4/9       # Tip

  # Chords
  c1 = 2*mean_c     # Head-to-tail chord
  c2 = 0.8*mean_c   # Root chord
  c3 = 1.2*mean_c   # Joint chord
  c4 = 0.7*mean_c   # Tip chord

  # Twists
  twist1 = 0*pi/180         # Head-to-tail chord twist
  twist2 = 0*pi/180         # Root twist
  twist3 = -2*pi/180        # Joint chord twist
  twist4 = -7*pi/180        # Tip chord twist

  # Sweeps
  lambda1 = 60*pi/180       # Head-to-root sweep
  lambda2 = 0*pi/180        # Root-to-joint sweep
  lambda3 = 20*pi/180       # Joint-to-tip sweep

  # Dihedral
  gamma1 = 0*pi/180         # Head-to-root dihedral
  # gamma2 = 4*pi/180         # Root-to-joint dihedral
  # gamma3 = 15*pi/180        # Joint-to-tip dihedral

  # Lattices in each section
  n1 = n*2          # Body
  n2 = n*10         # Forearm
  n3 = n*10         # Tip
  # Lattice expansion factor
  r1, r1_central = 1.0, false          # Head-to-root
  r2, r2_central = 2.0, true           # Root-to-joint
  r3, r3_central = 4.0, true           # Joint-to-tip

  # Fluid domain
  fdom_x_len = mean_c*6     # x length
  fdom_y_len = b*1.5        # y length
  fdom_z_len = mean_c*1.25  # z length
  fdom_NDIVS = fdom_n*[20, 15, 10]       # Subdivisions
  # --------------------------------------------------------


  # ---------------- WING SECTIONS--------------------------
  # (Builds each section as an independent wing)
  # Body
  b_y_tip = b*b_sec1/2
  b_x_tip = b_y_tip*tan(lambda1)
  b_z_tip = b_y_tip*tan(gamma1)
  b_c_tip = c2
  b_twist_tip = twist2
  b_x_root = 0.0
  b_y_root = 0.0
  b_z_root = 0.0
  b_c_root = c1
  b_twist_root = twist1

  body = vlm.Wing(b_x_tip, -b_y_tip, b_z_tip, b_c_tip, b_twist_tip)
  vlm.addchord(body, b_x_root, b_y_root, b_z_root, b_c_root, b_twist_root,
                    n1; r=r1, central=r1_central)
  vlm.addchord(body, b_x_tip, b_y_tip, b_z_tip, b_c_tip, b_twist_tip,
                    n1; r=1/r1, central=r1_central)

  # Forearms
  f_y_tip = b*b_sec2/2
  f_x_tip = f_y_tip*tan(lambda2)
  f_z_tip = 0.0
  f_c_tip = c3
  f_twist_tip = twist3
  f_x_root = 0.0
  f_y_root = 0.0
  f_z_root = 0.0
  f_c_root = c2
  f_twist_root = twist2

  ## Left forearm (root at the origin)
  forearm_l = vlm.Wing(f_x_tip, -f_y_tip, f_z_tip, f_c_tip, f_twist_tip)
  vlm.addchord(forearm_l, f_x_root, f_y_root, f_z_root, f_c_root, f_twist_root,
                    n2; r=r2, central=r2_central)
  ## Right forearm (root at the origin)
  forearm_r = vlm.Wing(f_x_root, f_y_root, f_z_root, f_c_root, f_twist_root)
  vlm.addchord(forearm_r, f_x_tip, f_y_tip, f_z_tip, f_c_tip, f_twist_tip,
                    n2; r=r2_central ? r2 : 1/r2, central=r2_central)

  # Wing tips
  t_y_tip = b*b_sec3/2
  t_x_tip = t_y_tip*tan(lambda3)
  t_z_tip = 0.0
  t_c_tip = c4
  t_twist_tip = twist4
  t_x_root = 0.0
  t_y_root = 0.0
  t_z_root = 0.0
  t_c_root = c3
  t_twist_root = twist3

  ## Left wingtip (root at the origin)
  wingtip_l = vlm.Wing(t_x_tip, -t_y_tip, t_z_tip, t_c_tip, t_twist_tip)
  vlm.addchord(wingtip_l, t_x_root, t_y_root, t_z_root, t_c_root, t_twist_root,
                    n3; r=r3, central=r3_central)
  ## Left wingtip (root at the origin)
  wingtip_r = vlm.Wing(t_x_root, t_y_root, t_z_root, t_c_root, t_twist_root)
  vlm.addchord(wingtip_r, t_x_tip, t_y_tip, t_z_tip, t_c_tip, t_twist_tip,
                    n3; r=r3_central ? r3 : 1/r3, central=r3_central)
  # --------------------------------------------------------

  # ---------------- CONSTRUCTS THE SYSTEM -----------------
  # Position of each section
  ## Body
  b_O = [0.0, 0.0, 0.0]
  b_Oaxis = [1.0 0 0; 0 1.0 0; 0 0 1.0]
  ## Left forearm
  # fl_O = [b_x_tip, -b_y_tip, b_z_tip] + b_O
  fl_O = vlm.countertransform([b_x_tip, -b_y_tip, b_z_tip], inv(b_Oaxis), b_O)
  fl_Oaxis = [1.0 0 0; 0 cos(gamma2) -sin(gamma2); 0 sin(gamma2) cos(gamma2)]
  ## Right forearm
  # fr_O = [b_x_tip, b_y_tip, b_z_tip] + b_O
  fr_O = vlm.countertransform([b_x_tip, b_y_tip, b_z_tip], inv(b_Oaxis), b_O)
  fr_Oaxis = [1.0 0 0; 0 cos(gamma2) sin(gamma2); 0 -sin(gamma2) cos(gamma2)]
  ## Left wingtip
  # tl_O = [f_x_tip, -f_y_tip, f_z_tip] + fl_O
  tl_O = vlm.countertransform([f_x_tip, -f_y_tip, f_z_tip], inv(fl_Oaxis), fl_O)
  tl_Oaxis = [1.0 0 0; 0 cos(gamma3) -sin(gamma3); 0 sin(gamma3) cos(gamma3)]
  ## Right wingtip
  # tr_O = [f_x_tip, f_y_tip, f_z_tip] + fr_O
  tr_O = vlm.countertransform([f_x_tip, f_y_tip, f_z_tip], inv(fr_Oaxis), fr_O)
  tr_Oaxis = [1.0 0 0; 0 cos(gamma3) sin(gamma3); 0 -sin(gamma3) cos(gamma3)]

  vlm.setcoordsystem(body, b_O, b_Oaxis)
  vlm.setcoordsystem(forearm_l, fl_O, fl_Oaxis)
  vlm.setcoordsystem(forearm_r, fr_O, fr_Oaxis)
  vlm.setcoordsystem(wingtip_l, tl_O, tl_Oaxis)
  vlm.setcoordsystem(wingtip_r, tr_O, tr_Oaxis)

  # System
  system = vlm.WingSystem()
  vlm.addwing(system, "Body", body)
  vlm.addwing(system, "ForearmLeft", forearm_l)
  vlm.addwing(system, "ForearmRight", forearm_r)
  vlm.addwing(system, "WingTipLeft", wingtip_l)
  vlm.addwing(system, "WingTipRight", wingtip_r)
  vlm.setVinf(system, Vinf)
  sys_O = [0.0, 0.0, 0.0]
  sys_Oaxis = [cos(AOA) 0 -sin(AOA); 0 1 0; sin(AOA) 0 cos(AOA)]
  vlm.setcoordsystem(system, sys_O, sys_Oaxis)
  # --------------------------------------------------------


  # ---------------- FLUID DOMAIN --------------------------
  P_min = [0, -fdom_y_len/2, 0]
  P_max = [fdom_x_len, fdom_y_len/2, fdom_z_len]
  fdom = vlm.PP.FluidDomain(P_min, P_max, fdom_NDIVS)
  fdom_O = vlm.countertransform(b_O , inv(sys_Oaxis), sys_O)
  fdom_O += -[1, 0, 0]*fdom_x_len/14  - [0, 0, 1]*fdom_z_len*1/2
  # fdom_Oaxis = sys_Oaxis
  fdom_Oaxis = [1.0 0 0; 0 1 0; 0 0 1]
  vlm.PP.setcoordsystem(fdom, fdom_O, fdom_Oaxis)
  # --------------------------------------------------------


  # ---------------- SOLVES VLM AND FLUID DOMAIN -----------
  vlm.solve(system, Vinf; t=0.0)

  if generate_fdom
    fdom_V(X) = vlm.Vind(system, X; t=0.0) + system.Vinf(X,0.0)
    fdom_field = [Dict( "field_name" => "U",
                        "field_type" => "vector",
                        "field_function" => fdom_V
                       )]
    vlm.PP.calculate(fdom, fdom_field; verbose=true)
  end
  # --------------------------------------------------------


  # ---------------- OUTPUTS VTKS --------------------------
  if save_path!=nothing

    # Gets rid of the Gamma field
    # for wing in system.wings
    #   pop!(wing.sol, "Gamma")
    # end

    # Checks if folder already exists
    if isdir(save_path)
      if prompt
        inp1 = ""
        opts1 = ["y", "n"]
        while false==(inp1 in opts1)
          print("\n\nFolder $save_path already exists. Remove? (y/n) ")
          inp1 = readline()[1:end-1]
        end
        if inp1=="y"
          run(`rm $save_path -rf`)
          println("\n")
        else
          return
        end
      else
        run(`rm $save_path -rf`)
      end
    end
    run(`mkdir $save_path`)

    # Saves code
    if save_code!=""
      run(`cp $save_code $save_path`)
    end

    # Saves vtks
    vlm.save(system, run_name; num=num,
            save_horseshoes=save_horseshoes, path=save_path)

    # Saves fluid domain
    if save_fdom && generate_fdom
      vlm.PP.save(fdom, run_name; path=save_path, num=(num==nothing?-1:num))
    end
  end

  # --------------------------------------------------------

  # ---------------- CALLS PARAVIEW ------------------------
  if save_path!=nothing && openParaview
    fyles = ""
    for fyle in ["_ForearmLeft_vlm.vtk", "_ForearmRight_vlm.vtk",
                  "_WingTipLeft_vlm.vtk", "_WingTipRight_vlm.vtk",
                  "_Body_vlm.vtk"]
      fyles = fyles*run_name*fyle*";"
    end
    fyles = fyles*(save_fdom ? run_name*"_fdom.vtk;" : "")
    run(`paraview --data=$(save_path*fyles)`)
  end

  return system, fdom
end


function timeformat(time_beg, time_end)
  time_delta = Int(time_end-time_beg)
  hrs = Int(floor(time_delta/1000/60/60))
  mins = Int(floor((time_delta-hrs*60*60*1000)/1000/60))
  secs = Int(floor((time_delta-hrs*60*60*1000-mins*1000*60)/1000))
  return hrs,mins,secs
end
