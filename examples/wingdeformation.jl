vlm_path = "../"
include(vlm_path*"src/FLOWVLM.jl")
vlm = FLOWVLM


global module_path; module_path,_ = splitdir(@__FILE__);   # Path to this module
global save_horseshoes = true
global run_name = "wingdeform"




"""
  Example of modeling a wing that has been deformed during flutter. In this case
  we received chord positions along the deformed wing calculated from ASWING and
  we use this information to construct a FLOWVLM Wing.
"""
function wingdeformation(; save_path=module_path*"/../temps/wingdeformation00/",
                                  paraview=true, verbose=false, delete=false,
                                  prompt=false)

println("Rabbit1: vars def")
prev_t = time()
  # ------------------- GEOMETRY FROM ASWING -----------------------------------
  # Wing
  w_x=[3.70114, 3.68389, 3.63262, 3.54877, 3.43465, 3.29332, 3.12842, 2.94402,
            2.74437, 2.53369, 2.36651, 2.13647, 1.90578, 1.67749, 1.45376,
            1.23582, 1.02398, 0.849822, 0.678993, 0.510528, 0.5, 0.510528,
            0.678993, 0.849822, 1.02398, 1.23582, 1.45376, 1.67749, 1.90578,
            2.13647, 2.36651, 2.53369, 2.74437, 2.94402, 3.12842, 3.29332,
            3.43465, 3.54877, 3.63262, 3.68389]
  w_y=[-20.4679, -20.3627, -20.0503, -19.5392, -18.8436, -17.982, -16.9767,
            -15.8524, -14.6347, -13.3495, -12.3294, -10.925, -9.51598,
            -8.12092, -6.75307, -5.42009, -4.12408, -3.05835, -2.01297,
            -0.982085, 0.0, 0.982085, 2.01297, 3.05835, 4.12408, 5.42009,
            6.75307, 8.12092, 9.51598, 10.925, 12.3294, 13.3495, 14.6347,
            15.8524, 16.9767, 17.982, 18.8436, 19.5392, 20.0503, 20.3627]
  w_z=[2.2904, 2.27214, 2.2179, 2.12927, 2.00887, 1.86028, 1.68818, 1.49839,
            1.29799, 1.09526, 0.943213, 0.749407, 0.573784, 0.418929,
            0.286959, 0.17961, 0.0983574, 0.0499936, 0.0188427, 0.0034735,
            -2.79235e-28, 0.0034735, 0.0188427, 0.0499936, 0.0983574, 0.17961,
            0.286959, 0.418929, 0.573784, 0.749407, 0.943213, 1.09526, 1.29799,
            1.49839, 1.68818, 1.86028, 2.00887, 2.12927, 2.2179, 2.27214]
  w_chord=[1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4,
            1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 2.0, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4,
            1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4]
  w_twist=[-0.0766763, -0.0766224, -0.0764595, -0.0761964, -0.0758377,
            -0.0753736, -0.0747697, -0.0739611, -0.0728547, -0.0713373,
            -0.0698069, -0.0672505, -0.0647935, -0.0623463, -0.0598169,
            -0.0571208, -0.0541716, -0.051353, -0.0487555, -0.0457306, 0.0,
            -0.0457294, -0.0487555, -0.051353, -0.0541716, -0.0571208,
            -0.0598169, -0.0623463, -0.0647935, -0.0672505, -0.0698069,
            -0.0713373, -0.0728547, -0.0739611, -0.0747697, -0.0753736,
            -0.0758377, -0.0761964, -0.0764595, -0.0766224]

  # Left winglet
  l_x=[3.51592, 3.52455, 3.54638, 3.57264, 3.59713, 3.61993, 3.64442, 3.67068,
            3.6925]
  l_y=[-20.4742, -20.4739, -20.4733, -20.4726, -20.4719, -20.4712, -20.4703,
            -20.4693, -20.4683]
  l_z=[4.50915, 4.40566, 4.14428, 3.82968, 3.53639, 3.26321, 2.96992, 2.6553,
            2.39389]
  l_chord=[0.5, 0.541978, 0.648009, 0.775625, 0.894594, 1.00541, 1.12437,
            1.25199, 1.35802]
  l_twist=[-0.0464419, -0.0462372, -0.0346505, -0.0199209, -0.00619664,
            -0.00621593, -0.0199834, -0.0347656, -0.0464038]

  # Right winglet
  r_x=[3.70114, 3.6925, 3.67068, 3.64442, 3.61993, 3.59713, 3.57264, 3.54638,
            3.52455]
  r_y=[20.4679, 20.4683, 20.4693, 20.4703, 20.4712, 20.4719, 20.4726, 20.4733,
            20.4739]
  r_z=[2.2904, 2.39389, 2.6553, 2.96992, 3.26321, 3.53639, 3.82968, 4.14428,
            4.40566]
  r_chord=[1.4, 1.35802, 1.25199, 1.12437, 1.00541, 0.894594, 0.775625,
            0.648009, 0.541978]
  r_twist=[-0.0466301, -0.0464038, -0.0347656, -0.0199834, -0.00621593,
            -0.00619664, -0.0199209, -0.0346505, -0.0462372]

  # NOTE: ASWING's twist is defined in relation to the s-axis, an axis which
  #         runs normal to the airfoil cross sections.  For the undeformed
  #         state, it was at an eighty degree dihedral angle.
  wl_angle = 80.0                       # (deg) winglet plane inclination

  # NOTE: Here I join each winglet with respective wing tips if needed
  if [l_x[end], l_y[end], l_z[end]]!=[w_x[1], w_y[1], w_z[1]]
    push!(l_x, w_x[1])
    push!(l_y, w_y[1])
    push!(l_z, w_z[1])
    push!(l_chord, w_chord[1])
    push!(l_twist, w_twist[1])
  end
  if [r_x[1], r_y[1], r_z[1]]!=[w_x[end], w_y[end], w_z[end]]
    r_x = vcat(w_x[end], r_x)
    r_y = vcat(w_y[end], r_y)
    r_z = vcat(w_z[end], r_z)
    r_chord = vcat(w_chord[end], r_chord)
    r_twist = vcat(w_twist[end], r_twist)
  end

  # ------------------- FLOWVLM PARAMETERS -------------------------------------
  magVinf = 1.0                         # Freestream magnitude
  Vinf(X,t) = magVinf*[1,0,0]           # Freestream
  Sref = 1.0                            # Reference area for coefficients

  # ------------------- PRECALCULATIONS ----------------------------------------
  nchords_w = size(w_y, 1)              # Number of chords along wing
  nchords_lwl = size(l_y, 1)            # Number of chords along left winglet
  nchords_rwl = size(r_y, 1)            # Number of chords along right winglet

  # Position and inclination of winglets
  O_lwl = [l_x[end], l_y[end], l_z[end]]# Centers it about its joint
  Oaxis_lwl = vlm.vtk.rotation_matrix(0.0, 0.0, wl_angle)  # Inclination

  O_rwl = [r_x[1], r_y[1], r_z[1]]      # Centers it about its joint
  Oaxis_rwl = vlm.vtk.rotation_matrix(0.0, 0.0, -wl_angle)  # Inclination
println("\ttime: $(time()-prev_t) s\n")


println("Rabbit2: FLOWVLM geometry def")
prev_t = time()
  # ------------------- FLOWVLM WINGSYSTEM -------------------------------------
  # Builds wing
  wing = nothing                        # Wing object
  for i in 1:nchords_w
    x, y, z, c, twist, n = w_x[i], w_y[i], w_z[i], w_chord[i], w_twist[i], 1

    if i==1 # Initiate Wing object with left-furthest chord
      wing = vlm.Wing(x, y, z, c, twist)

    else  # Add chords
      vlm.addchord(wing, x, y, z, c, twist, n)
    end
  end


  # Builds left winglet
  lwinglet = nothing                    # Wing object
  for i in 1:nchords_lwl
    # NOTE: Here I define the span in the z-component and center it about joint
    x = l_x[i]-O_lwl[1]
    y = -(l_z[i]-O_lwl[3])
    z = -(l_y[i]-O_lwl[2])
    c, twist, n = l_chord[i], l_twist[i], 1

    if i==1 # Initiate Wing object with left-furthest chord
      lwinglet = vlm.Wing(x, y, z, c, twist)

    else  # Add chords
      vlm.addchord(lwinglet, x, y, z, c, twist, n)
    end
  end
  vlm.setcoordsystem(lwinglet, O_lwl, Oaxis_lwl) # Places it at joint

  # Builds right winglet
  rwinglet = nothing                    # Wing object
  for i in 1:nchords_rwl
    # NOTE: Here I define the span in the z-component and center it about joint
    x = r_x[i]-O_rwl[1]
    y = (r_z[i]-O_rwl[3])
    z = -(r_y[i]-O_rwl[2])
    c, twist, n = r_chord[i], r_twist[i], 1

    if i==1 # Initiate Wing object with left-furthest chord
      rwinglet = vlm.Wing(x, y, z, c, twist)

    else  # Add chords
      vlm.addchord(rwinglet, x, y, z, c, twist, n)
    end
  end
  vlm.setcoordsystem(rwinglet, O_rwl, Oaxis_rwl) # Places it at joint

  # Builds the WingSystem
  system = vlm.WingSystem()
  vlm.addwing(system, "Wing", wing)
  vlm.addwing(system, "LWinglet", lwinglet)
  vlm.addwing(system, "RWinglet", rwinglet)
println("\ttime: $(time()-prev_t) s\n")


println("Rabbit3: FLOWVLM solver")
prev_t = time()
  # ------------ RUN FLOWVLM ---------------------------------------------------
  # Solves the lattice
  vlm.solve(system, Vinf)
println("\ttime: $(time()-prev_t) s\n")
println("Rabbit4: Aerodynamic properties")
prev_t = time()

  # Calculates aerodynamic properties
  vlm.calculate_field(system, "CFtot"; S=Sref)

println("\ttime: $(time()-prev_t) s\n")
println("Rabbit5: Generate vtk")
prev_t = time()
  # Outputs vtk geometry
  if save_path!=nothing
    vlm.vtk.create_path(save_path, prompt)
    vlm.save(system, run_name;
              save_horseshoes=save_horseshoes, path=save_path)
  end
println("\ttime: $(time()-prev_t) s\n")

  # ------------ VISUALIZATION -------------------------------------------------
  # Calls Paraview and deletes temp files
  if save_path!=nothing

    if paraview
      strn = ""
      for name in system.wing_names
        strn = strn * run_name * "_" * name * "_vlm.vtk;"
      end
      run(`paraview --data="$save_path$strn"`)
    end

    if delete
      run(`rm -rf $save_path`)
    end

  end

end
