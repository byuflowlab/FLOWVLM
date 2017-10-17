include("../src/FLOWVLM.jl")
vlm = FLOWVLM

using PyPlot

prompt = true
run_name = "propeller"
save_horseshoes = true
plot_rel_vel = false

function def_Vinf(X,t)
  return [1.0, 0, 0]         # Incoming free stream
end

global propeller = nothing

function singlepropeller(; save_path="temp_singlepropeller/")

  # ---------------- GEOMETRIC PARAMETERS -------------------------
  CW = true                         # Clockwise rotation
                                    # Radial positions
  r = .0254*[0.7526, 0.7928, 0.8329, 0.8731, 0.9132, 0.9586, 1.0332,
       1.1128, 1.1925, 1.2722, 1.3519, 1.4316, 1.5114, 1.5911,
       1.6708, 1.7505, 1.8302, 1.9099, 1.9896, 2.0693, 2.1490, 2.2287,
       2.3084, 2.3881, 2.4678, 2.5475, 2.6273, 2.7070, 2.7867, 2.8661, 2.9410]
                                    # Chord at each radial position
  chord = .0254*[0.6270, 0.6255, 0.6231, 0.6199, 0.6165, 0.6125, 0.6054, 0.5973,
            0.5887,
            0.5794, 0.5695, 0.5590, 0.5479, 0.5362, 0.5240, 0.5111, 0.4977,
            0.4836, 0.4689, 0.4537, 0.4379, 0.4214, 0.4044, 0.3867, 0.3685,
            0.3497, 0.3303, 0.3103, 0.2897, 0.2618, 0.1920]
                                    # Twist at each radial position
  theta = [40.2273, 38.7657, 37.3913, 36.0981, 34.8803, 33.5899, 31.6400,
                     29.7730, 28.0952, 26.5833, 25.2155, 23.9736, 22.8421,
                     21.8075,
                     20.8586, 19.9855, 19.1800, 18.4347, 17.7434, 17.1005,
                     16.5013,
                     15.9417, 15.4179, 14.9266, 14.4650, 14.0306, 13.6210,
                     13.2343,
                     12.8685, 12.5233, 12.2138]
  LE_x = zeros(r)                   # x-position of leading edge at radial pos
  LE_z = zeros(r)                   # z-position of leading edge at radial pos
  B = 3                             # Number of blades

  # ---------------- VLM PARAMETERS -------------------------------
  n = 25                            # Number of lattices per blade

  # ---------------- SIMULATION PARAMETERS -------------------------
  RPM = 600                         # Revs per minutes
  omega = 2*pi*RPM/60 * (-1)^(CW==false)   # (rad/s) angular velocity
  cycles = 6                        # Number of rotations
  cycles_rampup = 3                 # Number of rampup rotations to RPM
  nsteps_per_cycle = 30             # Steps per each rotation
  nsteps = cycles*nsteps_per_cycle  # Total number of steps
  init_angle = 0*pi/180             # (rad) initial rotation
  # d_angle = 2*pi/nsteps_per_cycle * (-1)^(CW==false)    # (rad) angle step
  dt = 1/(nsteps_per_cycle*RPM/60)  # (s) time step
  t_rampup = dt * nsteps_per_cycle*cycles_rampup        # Rampup time



  # ---------------- PREPARATION ----------------------------------
  # VLM setup
  global propeller = vlm.Rotor(CW, r, chord, theta, LE_x, LE_z, B)
  vlm.initialize(propeller, n)
  # vlm.setVinf(propeller, def_Vinf)

  # Creates save path
  if save_path!=nothing; vlm.create_path(save_path, prompt); end;

  # Angular velocity function
  # omega_fun(t) = omega
  function omega_fun(t)
    if t<t_rampup
      # Using a Weibull distribution to simulate rampup
      out = omega*(1 - exp(-(1.5*t/t_rampup)^5))
    else
      out = omega
    end
    return out
  end

  # ---------------- RUN -----------------------------------------
  prev_rotation = init_angle
  for i in 0:nsteps
    t = dt*i
    rotation = prev_rotation + omega_fun(t)*dt
    Oaxis = [1 0 0;
              0 cos(rotation) sin(rotation);
              0 -sin(rotation) cos(rotation)]

    # Rotates the propeller
    vlm.setcoordsystem(propeller, [0.0,0,0], Oaxis)
    # Resets its horseshoes to force recalculation
    vlm._reset(propeller; keep_Vinf=true)
    # Solves
    vlm.solve(propeller._wingsystem, def_Vinf; t=t,
              extraVinf=local_relative_vel, wingsystem=propeller._wingsystem,
              omega_fun=omega_fun)

    # vlm.calculate_field(propeller._wingsystem)

    vlm.save(propeller, run_name; save_horseshoes=save_horseshoes,
              path=save_path, num=i, t=dt)

    prev_rotation = rotation
  end

  if plot_rel_vel
    plot_relative_vel(; save_path=save_path,
                      ylims=[-maximum(r)*omega, maximum(r)*omega])
  end

  strn = ""
  for i in 1:B
    strn = strn * run_name * "_Blade$(i)_vlm...vtk;"
  end
  run(`paraview --data="$save_path$strn"`)
end

"Local relative velocity due to rotation"
global plots = Dict()
function local_relative_vel(i, t; wingsystem=nothing, omega_fun=nothing, wing=nothing)
  CP = vlm.getControlPoint(wing==nothing ? wingsystem : wing, i)  # i-th control point
  omega = omega_fun(t)                     # Angular velocity
  if wing==nothing
    _wing, _ = vlm._fetch_wing(wingsystem, i) # Wing of the i-th lattice
  else
    _wing = wing
  end
  r = norm(vlm.transform(CP, _wing.invOaxis, _wing.O))    # Radius
  local_V = Float64[omega*r, 0, 0]         # Velocity in local ref system
                                           # Velocity in global ref system
  global_V = vlm.countertransform(local_V, _wing.invOaxis, zeros(Float64, 3))
  if false==(t in keys(plots))
    plots[t] = Dict()
  end
  if false==(_wing in keys(plots[t]))
    this_dict = Dict(("x"=>[], "y"=>[], "label"=>"Blade$(length(plots[t])+1)"))
    plots[t][_wing] = this_dict
  end
  push!(plots[t][_wing]["x"], r)
  push!(plots[t][_wing]["y"], global_V)
  return global_V
end

function plot_relative_vel(; save_path=nothing, ylims=nothing)
  if save_path!=nothing; vlm.create_path(save_path*"plot1", false); end;
  j = 0
  for t in sort(collect(keys(plots)))
    wings = plots[t]
    fig = figure("local_flapping_vel", figsize=(7*3, 5*2))
    suptitle("Time = $(round(t,3)) (s)")
    fig_labels = ["x", "y", "z"]

    for i in 1:3
      subplot(230+i)
      title("Rotation $(fig_labels[i])-velocity")
      for (key,val) in wings
        plot(val["x"], [vel[i] for vel in val["y"]], "o", label=val["label"])
      end
      xlabel("r-distance")
      ylabel("local relative $(fig_labels[i])-velocity")
      legend(loc="best")
      grid(true, color="0.8", linestyle="--")
      if ylims!=nothing; ylim(ylims); end;
    end

    subplot(230+4)
    title("Rotation velocity magnitude")
    for (key,val) in wings
      plot(val["x"], [norm(vel) for vel in val["y"]], "o", label=val["label"])
    end
    xlabel("r-distance")
    ylabel("local relative velocity")
    legend(loc="best")
    grid(true, color="0.8", linestyle="--")
    if ylims!=nothing; ylim(ylims); end;

    if save_path!=nothing
      num = (j<10 ? "000" : j<100 ? "00" : j<1000 ? "0" : "")*"$j"
      savefig(save_path*"plot1/"*run_name*num*".png")
    end
    j += 1
    clf()
  end
end
