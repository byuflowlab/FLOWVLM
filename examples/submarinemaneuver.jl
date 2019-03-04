# Required packages
import FLOWVLM
import MyVPM
vlm = FLOWVLM
vpm = MyVPM

using PyPlot


module_path = splitdir(@__FILE__)[1]      # Path to this module




function submarinemaneuver(;
                            # SIMULATION OPTIONS
                            nsteps=500,
                            # OUTPUT OPTIONS
                            run_name="submarine",
                            save_path=joinpath(module_path, "../temps/submarine00"),
                            prompt=!true,
                            paraview=true,
                            verbose=true, v_lvl=0,
                            disp_plot=true,
                            save_horseshoes=!true,
                            )

    # ---------- SIMULATION PARAMETERS -----------------------------------------
    RPM = 600                   # Propeller RPM
    Vs = 10.0                   # (m/s) ship forward speed
    Dp = 1.0                    # (m) propeller diameter
    Ls = 17*Dp                  # (m) ship length
    Ds = 1.5*Dp                 # (m) ship diameter
    N::Int64 = 10               # Lattice refinement parameter

    tend = 7*Ls/Vs              # (s) simulation end time
    # nsteps = 500               # Number of time steps in simulation
    dt = tend/nsteps            # (s) time step size
    nsteps_save = 1             # Steps in between vtk outputs

    solver = "direct"           # VPM solver method

    Vinf(x,t) = [1e-14, 0, 0]        # Freestream velocity

    if verbose; println("\t"^v_lvl*"Revs in simulation: $(RPM/60*tend)"); end;


    # ---------- MANEUVER DEFINITION -------------------------------------------

    # Maneuvers: I  - Straight forward travel
    #            II - Foward travel lifting up
    #            III- Foward travel lifting up on pitch
    #            IV - Turn left
    #            V  - Straight forward travel
    # Ship starts facing -x, aligned with x-axis, and sail aligned with +z

    nman = 5                                     # Number of maneuvers
    Lmans = Ls*[1, 0.5, 1.5, 2, 2]               # Distance traveled in each maneuver
    tmans = [sum(Lmans[1:i]) for i in 1:nman]/Vs # Time interval in each maneuver

    VI = Vs*[-1, 0, 0]                           # Nominal velocity at stage
    VII = Vs*[-1, 0, 0.25]
    VV = Vs*[0, -1, 0]

    OaxisI = Float64[1 0 0; 0 0 1; 0 -1 0]       # Nominal orientation at stage
    OaxisV = Float64[0 1 0; 0 0 1; 1 0 0]

    spI = zeros(3)                               # Nominal sailplane angle at stage
    spII = Float64[0,-15,0]
    stI = zeros(3)                               # Nominal sternplane angle at stage
    stII = Float64[0,-10,0]
    stIII = Float64[0,-60,0]
    rdIII = Float64[0,0,90]                      # Nominal rudder angle at stage
    rdIV = Float64[0,-15,90]
    rdV = Float64[0,90,90]


    "Velocity of center of gravity in time"
    function velocity(t)
        if t<tmans[1]
            return VI
        elseif t<tmans[2]
            return VII
        elseif t<tmans[3]
            V = VII + [0, 0, 2*exp(-( (t-mean(tmans[2:3]))/((tmans[3]-tmans[2])/4) )^2)]
            return Vs*V/norm(V)
        elseif t<tmans[4]
            V = VII + (t-tmans[3])/(tmans[4]-tmans[3])*(VV - VII)
            return Vs*V/norm(V)
        else
            return VV
        end
    end

    "Submarine orientation in time"
    function orientation(t)
        if t<tmans[2]
            return OaxisI
        elseif t<tmans[3]
            a = 3*pi/32*exp(-( (t-mean(tmans[2:3]))/((tmans[3]-tmans[2])/4) )^2)
            return Float64[cos(a) 0 -sin(a); sin(a) 0 cos(a); 0 -1 0]
        elseif t<tmans[4]
            Oaxis = OaxisI + (t-tmans[3])/(tmans[4]-tmans[3])*(OaxisV - OaxisI)
            return hcat([Oaxis[i, 1:3]/norm(Oaxis[i, 1:3]) for i in 1:3]...)'
        else
            return OaxisV
        end
    end

    "Control surface deflections in time"
    function control_surface_angle(t)
        if t<tmans[1]  # Control surface => [yaw, pitch, roll] (degrees)
            return Dict("Sailplane"     => spI,
                        "RudderCS"      => nothing,
                        "SternplaneCS"  => stI)
        elseif t<tmans[2]
            auxsp = 0.2*(tmans[2]-tmans[1])
            sp = t-tmans[1]<auxsp ? spI + ((t-tmans[1])/auxsp)*(spII - spI) : spII
            st = t-tmans[1]<auxsp ? stI + ((t-tmans[1])/auxsp)*(stII - stI) : stII
            return Dict("Sailplane"     => sp,
                        "RudderCS"      => nothing,
                        "SternplaneCS"  => st)
        elseif t<tmans[3]
            auxsp = 0.8*(tmans[3]-tmans[2])
            sp = t-tmans[2]<auxsp ? spII : spII + ((t-auxsp)/(tmans[3]-auxsp))*(spI - spII)

            auxst1 = 0.2*(tmans[3]-tmans[2])
            auxst2 = 0.8*(tmans[3]-tmans[2])
            if t-tmans[2] < auxst1
                st = stII + ((t-tmans[2])/auxst1)*(stIII - stII)
            elseif t-tmans[2] < auxst2
                st = stIII + ((t-tmans[2])/auxst2)*(stII - stIII)
            else
                st = stII + ((t-auxst2)/(tmans[3]-auxst2))*(stI - stII)
            end
            return Dict("Sailplane"     => sp,
                        "RudderCS"      => nothing,
                        "SternplaneCS"  => st)
        elseif t<tmans[4]
            auxrd1 = 0.2*(tmans[4]-tmans[3])
            auxrd2 = 0.9*(tmans[4]-tmans[3])
            if t-tmans[3]<auxrd1
                rd = rdIII + ((t-tmans[3])/auxrd1)*(rdIV - rdIII)
            elseif t-tmans[3]<auxrd2
                rd = nothing
            else
                rd = rdIV + ((t-auxrd2)/(tmans[4]-auxrd2))*(rdV - rdIV)
            end
            return Dict("Sailplane"     => nothing,
                        "RudderCS"      => rd,
                        "SternplaneCS"  => nothing)
        else
            return Dict("Sailplane"     => nothing,
                        "RudderCS"      => rdV,
                        "SternplaneCS"  => nothing)
        end
    end

    if disp_plot
        for i in 1:nman
            plot(tmans[i]*ones(2), [-Vs, Vs], ":k", alpha=0.75)
        end
        ts = linspace(0, tend, nsteps)
        vs = velocity.(ts)
        plot(ts, [V[1] for V in vs], "-", alpha=0.75, label=L"$x$-component")
        plot(ts, [V[2] for V in vs], "-", alpha=0.75, label=L"$y$-component")
        plot(ts, [V[3] for V in vs], "-", alpha=0.75, label=L"$z$-component")
        xlabel("Time (s)")
        ylabel("Velocity (m/s)")
        title("Maneuvers")
        legend(loc="best", frameon=!false)
        grid(true, color="0.8", linestyle="--")
    end



    # ---------- SIMULATION SETUP ----------------------------------------------

    # Generate vortex lattice
    submarine = generate_submarine(Dp, Ls, Ds, N;
                                run_name=run_name,
                                save_path=nothing,
                                prompt=prompt,
                                paraview=paraview,
                                verbose=verbose, v_lvl=v_lvl+1,
                                save_horseshoes=save_horseshoes,
                                )
    liftsurfs = vlm.get_wing(submarine, "LF")

    # Initiate particle field
    pfield = vpm.ParticleField(Int(1e6), Vinf, nothing, solver)
    vpm.addparticle(pfield, [zeros(3)..., ones(3)..., 1.0, 1.0])





    vlm.setVinf(submarine, Vinf)


    # --------------- SETS UP RUNTIME ROUTINE ----------------------------------

    # Runtime function
    function runtime_function(PFIELD, T, DT)

        prev_liftsurfs = deepcopy(liftsurfs)

        # Deflect control surfaces
        cs_angle = control_surface_angle(T)
        for (wing, label) in [
                            (vlm.get_wing(liftsurfs, "Sailplane"), "Sailplane"),
                            (vlm.get_wing(liftsurfs, "RudderCS"), "RudderCS"),
                            (vlm.get_wing(liftsurfs, "SternplaneCS"), "SternplaneCS"),
                            ]
            if cs_angle[label]!=nothing
                yaw, pitch, roll = cs_angle[label]
                Oaxis = vlm.vtk.rotation_matrix(yaw, pitch, roll)

                vlm.setcoordsystem(wing, wing.O, Oaxis)
            end
        end

        # Move the submarine
        O_sub = submarine.O + velocity(T)*DT
        Oaxis_sub = orientation(T)

        vlm.setcoordsystem(submarine, O_sub, Oaxis_sub)

        # Recalculate horseshoes with local translational velocities
        for i in eachindex(liftsurfs.wings)
            cur = vlm.get_wing(liftsurfs, i)
            prev = vlm.get_wing(prev_liftsurfs, i)
            vlm._calculateHSs(cur; t=T, extraVinf=Vlocal, dt=DT,
                                                    prev_sys=prev, cur_sys=cur)
        end

        # Solves the VLM of lifting surfaces
        vlm.solve(liftsurfs, Vinf; t=T, extraVinf=Vlocal, dt=DT,
                                    prev_sys=prev_liftsurfs, cur_sys=liftsurfs)


        # Saves VLM
        if save_path!=nothing && PFIELD.nt%nsteps_save==0
            vlm.save(submarine, run_name; save_horseshoes=save_horseshoes,
                        path=save_path, num=PFIELD.nt, t=DT)
        end

        return false
    end


    # --------------- RUN THE VPM ----------------------------------------------
    vpm.run_vpm!(pfield, dt, nsteps; save_path=save_path, run_name=run_name,
                 runtime_function=runtime_function, solver_method=solver,
                 nsteps_save=nsteps_save,
                 save_code=joinpath(module_path,"submarinemaneuver.jl"),
                 prompt=prompt)



  # ---------- VISUALIZATION -------------------------------------------------
  if save_path != nothing
      if paraview

          strn = ""
          for (i,subsystem) in enumerate(submarine.wings)
              sysname = submarine.wing_names[i]
              for wname in subsystem.wing_names
                  strn *= run_name*"_"*sysname*"_"*wname*"_vlm...vtk;"
              end
          end

          run(`paraview --data=$save_path/$strn`)
      end
  end

  #   cur_rev = floor(PFIELD.nt/steps_per_rev) # Revolutions completed
  #
  #   # Rotates the propellers
  #   rotation = PFIELD.nt==0 ? 0.0 : dangle_prop
  #   for prop in props_w; vlm.rotate(prop, rotation); end;
  #   for prop in props_tw; vlm.rotate(prop, rotation); end;

  #
  #   # Calculates circulation from CCBlade
  #   for prop in props_w; vlm.solvefromCCBlade(prop, Vinf, RPM, rho; t=T); end;
  #   for prop in props_tw; vlm.solvefromCCBlade(prop, Vinf, RPM, rho; t=T); end;
  #
  #   # Adds particles
  #   if PFIELD.nt%n_steps_shedding==0
  #     for prop in props_w; VLM2VPM(prop, PFIELD, DT; t=T); end;
  #     for prop in props_tw; VLM2VPM(prop, PFIELD, DT; t=T); end;
  #   end
  #


end

"Returns the local translational velocity of a control point"
function Vlocal(i, t; prev_sys=nothing, cur_sys=nothing, dt=nothing, wing=nothing)
    prev_X = vlm.getControlPoint(prev_sys, i)
    cur_X = vlm.getControlPoint(cur_sys, i)
    return -(cur_X-prev_X)/dt
end

function generate_submarine(Dp, Ls, Ds, N;
                            # OUTPUT OPTIONS
                            run_name="submarine",
                            save_path=joinpath(module_path, "../temps/submarine00"),
                            prompt=!true,
                            paraview=true,
                            verbose=true, v_lvl=0,
                            save_horseshoes=!true,
                            )

    # ---------- GEOMETRIC PARAMETERS ------------------------------------------

    # Hull
    n_h = 1*N                               # Number of lattices
    b_h = Ds                                # Span (diameter)
    # ar_h = b_h/(0.9*Ls)                     # Aspect ratio (span/tip chord)
    pos_h = collect(linspace(0, 1, n_h))    # Position of chords along semi-span
    # clen_h = (2.0*sqrt.( 1 .- pos_h.^2 ) + 1.0)/(2.0 + 1.0) # Chord distribution
    # sweep_h = 45*sqrt.( 1 .- (1-pos_h[1:end-1]).^2 )      # Sweep distribution
    clen_h = (0.5*sqrt.( 1 .- pos_h.^2 ) + 0.25)/(0.5 + 0.25) # Chord distribution
    ar_h = (b_h/Ls) * clen_h[1]
    sweep_h = -35*sqrt.( 1 .- (1-pos_h[1:end-1]).^2 )      # Sweep distribution
    twist_h = zeros(size(pos_h))            # Twist distribution
    dihed_h = zeros(size(pos_h[1:end-1]))   # Dihedral distribution
    chordalign_h = 0.25                     # Chord alignment

    hull = vlm.complexWing(b_h, ar_h, n_h, pos_h, clen_h,
                            twist_h, sweep_h, dihed_h; chordalign=chordalign_h)

    # Sail
    n_sl = 1*N                              # Number of lattices
    b_sl = 2*0.6*Ds                         # Full span
    ar_sl = 1.0                             # Aspect ratio (span/tip chord)
    pos_sl = [0.0, 1.0]                     # Position of chords along semi-span
    clen_sl = [1.0, 0.75]                   # Chord distribution
    twist_sl = zeros(pos_sl)                # Twist distribution
    sweep_sl = zeros(pos_sl[1:end-1])       # Sweep distribution
    dihed_sl = zeros(pos_sl[1:end-1])       # Dihedral distribution
    chordalign_sl = 0.25                    # Chord alignment

    x_sl = 0.25*Ls                          # x-position of nose
    y_sl = Ds/2                             # y-position of nose
    z_sl = 0.0                              # z-position of nose

    sail = vlm.complexWing(b_sl, ar_sl, n_sl, pos_sl, clen_sl, twist_sl,
                                sweep_sl, dihed_sl; chordalign=chordalign_sl,
                                symmetric=false)

    # Sailplane
    n_sp = 2*N                              # Number of lattices
    b_sp = 0.85*Ds                          # Full span
    ar_sp = 5.0                             # Aspect ratio (span/tip chord)
    pos_sp = [0.0, 1.0]                     # Position of chords along semi-span
    clen_sp = [1.5, 1.0]                    # Chord distribution
    twist_sp = zeros(pos_sp)                # Twist distribution
    sweep_sp = 0*ones(pos_sp[1:end-1])      # Sweep distribution
    dihed_sp = 0*ones(pos_sp[1:end-1])      # Dihedral distribution
    chordalign_sp = 1.0                     # Chord alignment

    x_sp = x_sl + 0.25*(b_sl/2)             # x-position of nose
    y_sp = y_sl + (b_sl/2)*7/12             # y-position of nose
    z_sp = 0.0                              # z-position of nose

    sailplane = vlm.complexWing(b_sp, ar_sp, n_sp, pos_sp, clen_sp, twist_sp,
                                sweep_sp, dihed_sp; chordalign=chordalign_sp)

    # Rudder fixed surface
    f_rd = 0.95                             # Fraction position of rudder
    n_rdf = 2*N                             # Number of lattices
    b_rdf = 1.5*Ds                          # Full span
    ar_rdf = 5.0                            # Aspect ratio (span/tip chord)
    pos_rdf = [0.0, 1.0]                    # Position of chords along semi-span
    clen_rdf = [1.9, 1.0]                   # Chord distribution
    twist_rdf = zeros(pos_rdf)              # Twist distribution
    sweep_rdf = 0*ones(pos_rdf[1:end-1])    # Sweep distribution
    dihed_rdf = 0*ones(pos_rdf[1:end-1])    # Dihedral distribution
    chordalign_rdf = 1.0                    # Chord alignment

    x_rdf = f_rd*Ls - clen_rdf[1]*b_rdf/ar_rdf         # x-position of nose
    y_rdf = 0.0                             # y-position of nose
    z_rdf = 0.0                             # z-position of nose

    rudder_fix = vlm.complexWing(b_rdf, ar_rdf, n_rdf, pos_rdf, clen_rdf, twist_rdf,
                                sweep_rdf, dihed_rdf; chordalign=chordalign_rdf)


    # Rudder control surface
    n_rdc = n_rdf                           # Number of lattices
    b_rdc = b_rdf                           # Full span
    ar_rdc = 2.5*ar_rdf                     # Aspect ratio (span/tip chord)
    pos_rdc = pos_rdf                       # Position of chords along semi-span
    clen_rdc = clen_rdf                     # Chord distribution
    twist_rdc = twist_rdf                   # Twist distribution
    sweep_rdc = sweep_rdf                   # Sweep distribution
    dihed_rdc = dihed_rdf                   # Dihedral distribution
    chordalign_rdc = 0.0                    # Chord alignment

    x_rdc = x_rdf + clen_rdf[1]*b_rdf/ar_rdf         # x-position of nose
    y_rdc = y_rdf                           # y-position of nose
    z_rdc = z_rdf                           # z-position of nose

    rudder_cs = vlm.complexWing(b_rdc, ar_rdc, n_rdc, pos_rdc, clen_rdc, twist_rdc,
                                sweep_rdc, dihed_rdc; chordalign=chordalign_rdc)


    # Sternplane fixed surface
    n_stf = n_rdf                           # Number of lattices
    b_stf = b_rdf                           # Full span
    ar_stf = 1.2*ar_rdf                     # Aspect ratio (span/tip chord)
    pos_stf = pos_rdf                       # Position of chords along semi-span
    clen_stf = [2.2, 1.0]                     # Chord distribution
    twist_stf = twist_rdf                   # Twist distribution
    sweep_stf = sweep_rdf                   # Sweep distribution
    dihed_stf = dihed_rdf                   # Dihedral distribution
    chordalign_stf = chordalign_rdf         # Chord alignment

    x_stf = f_rd*Ls - clen_stf[1]*b_stf/ar_stf     # x-position of nose
    y_stf = 0.0                             # y-position of nose
    z_stf = 0.0                             # z-position of nose

    sternplane_fix = vlm.complexWing(b_stf, ar_stf, n_stf, pos_stf, clen_stf, twist_stf,
                                sweep_stf, dihed_stf; chordalign=chordalign_stf)


    # Sternplane control surface
    n_stc = n_stf                           # Number of lattices
    b_stc = b_stf                           # Full span
    ar_stc = 1.1*ar_stf                     # Aspect ratio (span/tip chord)
    pos_stc = pos_stf                       # Position of chords along semi-span
    clen_stc = [1.0, 1.0]                     # Chord distribution
    twist_stc = twist_stf                   # Twist distribution
    sweep_stc = sweep_stf                   # Sweep distribution
    dihed_stc = dihed_stf                   # Dihedral distribution
    chordalign_stc = 0.0                    # Chord alignment

    x_stc = x_stf + clen_stf[1]*b_stf/ar_stf         # x-position of nose
    y_stc = y_stf                           # y-position of nose
    z_stc = z_stf                           # z-position of nose

    sternplane_cs = vlm.complexWing(b_stc, ar_stc, n_stc, pos_stc, clen_stc, twist_stc,
                                sweep_stc, dihed_stc; chordalign=chordalign_stc)


    # ---------- CONSTRUCT SYSTEM ----------------------------------------------
    # Sail
    O_sl = [x_sl, y_sl, z_sl]               # Position of nose
    Oaxis_sl = eye(3)                       # Orientation

    vlm.setcoordsystem(sail, O_sl, Oaxis_sl)

    # Sailplane
    O_sp = [x_sp, y_sp, z_sp]               # Position of nose
    Oaxis_sp = Float64[1 0 0; 0 0 -1; 0 1 0]# Orientation

    vlm.setcoordsystem(sailplane, O_sp, Oaxis_sp)

    # Rudder fixed surface
    O_rdf = [x_rdf, y_rdf, z_rdf]           # Position of nose
    Oaxis_rdf = eye(3)                      # Orientation

    vlm.setcoordsystem(rudder_fix, O_rdf, Oaxis_rdf)

    # Rudder control surface
    O_rdc = [x_rdc, y_rdc, z_rdc]           # Position of nose
    Oaxis_rdc = eye(3)                      # Orientation

    vlm.setcoordsystem(rudder_cs, O_rdc, Oaxis_rdc)

    # Sternplane fixed surface
    O_stf = [x_stf, y_stf, z_stf]           # Position of nose
    Oaxis_stf = Float64[1 0 0; 0 0 -1; 0 1 0]# Orientation

    vlm.setcoordsystem(sternplane_fix, O_stf, Oaxis_stf)

    # Sternplane control surface
    O_stc = [x_stc, y_stc, z_stc]           # Position of nose
    Oaxis_stc = Float64[1 0 0; 0 0 -1; 0 1 0]# Orientation

    vlm.setcoordsystem(sternplane_cs, O_stc, Oaxis_stc)



    # System of lifting surfaces
    lfsystem = vlm.WingSystem()
    vlm.addwing(lfsystem, "Sail", sail)
    vlm.addwing(lfsystem, "Sailplane", sailplane)
    vlm.addwing(lfsystem, "RudderFixed", rudder_fix)
    vlm.addwing(lfsystem, "RudderCS", rudder_cs)
    vlm.addwing(lfsystem, "SternplaneFixed", sternplane_fix)
    vlm.addwing(lfsystem, "SternplaneCS", sternplane_cs)

    # System of non-lifting surfaces
    nlfsystem = vlm.WingSystem()
    vlm.addwing(nlfsystem, "Hull", hull)

    # Submarine system
    submarine = vlm.WingSystem()
    vlm.addwing(submarine, "NLF", nlfsystem)
    vlm.addwing(submarine, "LF", lfsystem)

    # Orientation of the system
    O_sub = zeros(3)                            # Origin
    Oaxis_sub = Float64[1 0 0; 0 0 1; 0 -1 0]   # Orientation

    vlm.setcoordsystem(submarine, O_sub, Oaxis_sub)




    # ---------- VISUALIZATION -------------------------------------------------
    if save_path != nothing

        Vinf(x,t) = [1e-14, 0, 0]
        vlm.setVinf(submarine, Vinf)

        vlm.vtk.create_path(save_path, prompt)
        vlm.save(submarine, run_name; save_horseshoes=save_horseshoes, path=save_path)

        if paraview

            strn = ""
            for (i,subsystem) in enumerate(submarine.wings)
                sysname = submarine.wing_names[i]
                for wname in subsystem.wing_names
                    strn *= run_name*"_"*sysname*"_"*wname*"_vlm.vtk;"
                end
            end

            run(`paraview --data=$save_path/$strn`)
        end
    end

    return submarine
end
