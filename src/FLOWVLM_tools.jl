
################################################################################
# CALCULATION WING PROPERTIES
################################################################################
"""
Returns the area of the m-th panel, with `self` either Wing or WingSystem.
For derivation, see notebook entry 20170619.
"""
function get_A(self, m::Int64)
  Ap, A, B, Bp, CP, infDA, infDB, Gamma = getHorseshoe(self, m)
  CPAB = (A+B)/2
  lcp = norm(CPAB-CP)/(pm-pn)
  h = sqrt( (A[2]-B[2])^2 + (A[3]-B[3])^2  )
  Apanel = lcp*h
  return Apanel
end

"""
Returns the center of gravity of the m-th panel, with `self` either Wing or
WingSystem.
"""
function get_r(self, m::Int64)
  Ap, A, B, Bp, CP, infDA, infDB, Gamma = getHorseshoe(self, m)
  CPAB = (A+B)/2
  r = (CPAB+CP)/2
  return r
end

"""
Returns the center of gravity of the wing, with `self` either Wing or
WingSystem.
"""
function get_CG(self)
  sum_A = 0.0
  sum_rA = zeros(3)
  # Iterates over each panel
  for i in 1:get_m(self)
    r = get_r(self, i)
    A = get_A(self, i)
    sum_A += A
    sum_rA += r*A
  end
  r_cg = sum_rA/sum_A
  return r_cg
end

"""
  Returns the mean chord length of a Wing or WingSystem
"""
function get_barc(self)
  sum_lA = 0.0
  sum_A = 0.0

  # Iterates over each panel
  for i in 1:get_m(self)
    Ap, A, B, Bp, CP, infDA, infDB, Gamma = getHorseshoe(self, i)
    CPAB = (A+B)/2

    # Mean chord of this panel
    lcp = norm(CPAB-CP)/(pm-pn)
    # Area of this panel
    Apanel = get_A(self, i)

    sum_A += Apanel
    sum_lA += lcp*Apanel
  end
  barc = sum_lA/sum_A
  return barc
end

"""
Receives a Wing or WingsSystem and return the planform area, which is the
area of the projection on the xy-plane of its local coordinate system.
"""
function planform_area(wing)
  area = 0.0

  # Calculates the area of each panel
  ip = wing.Oaxis[1,:]
  jp = wing.Oaxis[2,:]
  for i in 1:get_m(wing)
    # Obtains the points
    Ap, A, B, Bp, CP, infDA, infDB, Gamma = getHorseshoe(wing, i)
    # Projects the points to the xy-plane
    _A = dot(A,ip)*ip + dot(A,jp)*jp
    _B = dot(B,ip)*ip + dot(B,jp)*jp
    _CP = dot(CP,ip)*ip + dot(CP,jp)*jp
    # Calculates the area of the trapezoid
    CPab = (_A+_B)/2
    lcp = norm(CPab-_CP)/(pm-pn) # Chord at the control point
    h = sqrt( (_A[2]-_B[2])^2 + (_A[3]-_B[3])^2 ) # Trapezoid height
    Apnl = lcp*h # Area of the panel

    area += Apnl
  end

  return area
end


################################################################################
# UTILITIES
################################################################################
"Example of a system of lifting surfaces (WingSystem)"
function exampleWing(; n=1, solveVLM=true, aileron_angleL = 15*pi/180,
                      aileron_angleR = -15*pi/180, num=nothing,
                      openParaview=true)

  # SIMULATION PARAMETERS
  magVinf = 55.0
  AOA = 4.0*pi/180
  Vinf(X,t) = magVinf*[cos(AOA), 0, sin(AOA)]
  rhoinf = 9.093/10^1
  qinf = (1/2)*rhoinf*magVinf^2
  # aileron_angleL = 15*pi/180      # Inclination of left aileron
  # aileron_angleR = -15*pi/180     # Inclination of right aileron

  # VLM PARAMETERS
  ## Number of lattices
  n_w = Int(ceil(n*40))           # Semi-span of the wing
  n_wl = Int(ceil(n_w/2))         # Winglets
  n_f = Int(ceil(n*10))           # Fuselage
  n_vt = Int(ceil(n*20))          # Vertical tail
  n_ht = Int(ceil(n*20))          # Semi-span of horizontal tail
  ## Panel expansion ratio
  r_w = 5.0
  r_wl = 1.0
  r_f = 1.0
  r_vt = 4.0
  r_ht = 5.0


  # BUILDS EACH LIFTING SURFACE
  ## --------- Wing --------------------
  ### Parameters
  w_b = 40.0                      # Span
  w_lambda = 30.0                 # Leading edge sweep
  w_gamma = 5.0                   # Dihedral
  w_twist_root = 0.0              # Twist at the root
  w_twist_tip = -2.5              # Twist at the tip
  w_c_root = 10.0                 # Chord at the root
  w_c_tip = 5.0                   # Chord at the tip
  w_sec_div = 0.60                # Spanwise position of the aileron
  w_aileron_c = 0.15*w_c_tip      # Chord of the aileron
  n_w1 = Int(ceil(n_w*w_sec_div*3/4))
  n_w2 = n_w-n_w1
  n_a = n_w2
  r_w1, r_w2, r_a = r_w, r_w, r_w
  ### Dimensions inner section
  w1_y_tip = (w_b*w_sec_div)/2                    # y-position of the tip chord
  w1_x_tip = w1_y_tip*tan(w_lambda*pi/180)        # x-position of the tip chord
  w1_z_tip = w1_y_tip*tan(w_gamma*pi/180)         # z-position of the tip chord
  w1_c_tip = w_c_root + w_sec_div*(w_c_tip-w_c_root)  # Chord length at the tip
  w1_c_angle_tip = w_sec_div*w_twist_tip          # Twist at the tip
  w1_x_root = 0.0
  w1_y_root = 0.0
  w1_z_root = 0.0
  w1_c_root = w_c_root
  w1_c_angle_root = w_twist_root
  ### Dimensions outer section
  w2_x_root = w1_x_tip
  w2_y_root = w1_y_tip
  w2_z_root = w1_z_tip
  w2_c_root = w1_c_tip-w_aileron_c
  w2_c_angle_root = w1_c_angle_tip
  w2_y_tip = w_b/2
  w2_x_tip = w2_y_tip*tan(w_lambda*pi/180)
  w2_z_tip = w2_y_tip*tan(w_gamma*pi/180)
  w2_c_tip = w_c_tip-w_aileron_c
  w2_c_angle_tip = w_twist_tip
  ### Dimensions aileron
  a_OL = [w2_x_root+w2_c_root*cos(w2_c_angle_root*pi/180),# Origin of left aileron
          -w2_y_root,
          w2_z_root-w2_c_root*sin(w2_c_angle_root*pi/180)]
  a_OR = a_OL.*[1, -1, 1]                                # Origin of right aileron
  a_x_root = 0.0
  a_y_root = 0.0
  a_z_root = 0.0
  a_c_root = w_aileron_c
  a_c_angle_root = w2_c_angle_root
  a_x_tip = w2_x_tip+w2_c_tip*cos(w2_c_angle_tip*pi/180) - a_OL[1]
  a_y_tip = w_b*(1-w_sec_div)/2
  a_z_tip = w2_z_tip-w2_c_tip*sin(w2_c_angle_tip*pi/180) - a_OL[3]
  a_c_tip = w_aileron_c
  a_c_angle_tip = w2_c_angle_tip

  ### Sets the angle of inclination of the aileron through its axis orientation
  aux1 = [a_x_tip, -a_y_tip, a_z_tip]
  a_OaxisL = axis_rotation(aux1/norm(aux1), -aileron_angleL*180/pi)
  aux1 = [a_x_tip, a_y_tip, a_z_tip]
  a_OaxisR = axis_rotation(aux1/norm(aux1), aileron_angleR*180/pi)

  ### Creates lifting surfaces of the wing
  wing_inner = Wing(w1_x_tip, -w1_y_tip, w1_z_tip, w1_c_tip, w1_c_angle_tip)
  addchord(wing_inner, w1_x_root, w1_y_root, w1_z_root, w1_c_root,
            w1_c_angle_root, n_w1; r=r_w1)
  addchord(wing_inner, w1_x_tip, w1_y_tip, w1_z_tip, w1_c_tip,
            w1_c_angle_tip, n_w1; r=1/r_w1)

  wing_outerL = Wing(w2_x_tip, -w2_y_tip, w2_z_tip, w2_c_tip, w2_c_angle_tip)
  addchord(wing_outerL, w2_x_root, -w2_y_root, w2_z_root, w2_c_root,
            w2_c_angle_root, n_w2; r=r_w2, central=true)

  wing_outerR = Wing(w2_x_root, w2_y_root, w2_z_root, w2_c_root,w2_c_angle_root)
  addchord(wing_outerR, w2_x_tip, w2_y_tip, w2_z_tip, w2_c_tip,
            w2_c_angle_tip, n_w2; r=r_w2, central=true)

  aileronL = Wing(a_x_tip, -a_y_tip, a_z_tip, a_c_tip, a_c_angle_tip)
  addchord(aileronL, a_x_root, a_y_root, a_z_root, a_c_root,
            a_c_angle_root, n_a; r=r_a, central=true)
  setcoordsystem(aileronL, a_OL, a_OaxisL)

  aileronR = Wing(a_x_root, a_y_root, a_z_root, a_c_root, a_c_angle_root)
  addchord(aileronR, a_x_tip, a_y_tip, a_z_tip, a_c_tip,
            a_c_angle_tip, n_a; r=r_a, central=true)
  setcoordsystem(aileronR, a_OR, a_OaxisR)


  ## --------- C-WINGLETS --------------------
  ### Parameters
  wl1_lambda = 87.0               # First section swept-back angle
  wl2_lambda = 25.0               # Second section swept-back
  wl1_gamma = 45.0                # First section dihedral
  wl2_gamma = 180.0               # Second section dihedral
  wl1_l = w_b/8                   # First section length
  wl2_l = wl1_l*1.25              # Second section length
  wl1_c_tip = w_c_tip*0.5         # First section tip chord
  wl2_c_tip = wl1_c_tip*0.75      # Second section tip chord
  n_wl1 = Int(ceil(n_wl*wl1_l/(wl1_l+wl2_l)))
  n_wl2 = n_wl - n_wl1
  r_wl1, r_wl2 = r_wl, r_wl
  ### Dimensions first section
  wl1_x_root = w2_x_tip
  wl1_y_root = w2_y_tip
  wl1_z_root = w2_z_tip
  wl1_c_root = w_c_tip
  wl1_c_angle_root = w2_c_angle_tip
  wl1_x_tip = wl1_x_root + wl1_l*cos(wl1_gamma*pi/180)*sin(wl1_lambda*pi/180)
  wl1_y_tip = wl1_y_root + wl1_l*cos(wl1_gamma*pi/180)*cos(wl1_lambda*pi/180)
  wl1_z_tip = wl1_z_root + wl1_l*sin(wl1_gamma*pi/180)
  wl1_c_tip = wl1_c_tip
  wl1_c_angle_tip = 0.0
  ### Dimensions second section
  wl2_x_root = wl1_x_tip
  wl2_y_root = wl1_y_tip
  wl2_z_root = wl1_z_tip
  wl2_c_root = wl1_c_tip
  wl2_c_angle_root = wl1_c_angle_tip
  wl2_x_tip = wl2_x_root + wl2_l*abs(cos(wl2_gamma*pi/180))*sin(wl2_lambda*pi/180)
  wl2_y_tip = wl2_y_root + wl2_l*cos(wl2_gamma*pi/180)*cos(wl2_lambda*pi/180)
  wl2_z_tip = wl2_z_root + wl2_l*sin(wl2_gamma*pi/180)
  wl2_c_tip = wl2_c_tip
  wl2_c_angle_tip = 0.0
  ### Creates the lifting surfaces
  wingletL = Wing(wl2_x_tip, -wl2_y_tip, wl2_z_tip, wl2_c_tip, wl2_c_angle_tip)
  addchord(wingletL, wl1_x_tip, -wl1_y_tip, wl1_z_tip, wl1_c_tip,
    wl1_c_angle_tip, n_wl2; r=r_wl2)
  addchord(wingletL, wl1_x_root, -wl1_y_root, wl1_z_root, wl1_c_root,
            wl1_c_angle_root, n_wl1; r=r_wl1)
  wingletR = Wing(wl1_x_root, wl1_y_root, wl1_z_root, wl1_c_root,
            wl1_c_angle_root)
  addchord(wingletR, wl1_x_tip, wl1_y_tip, wl1_z_tip, wl1_c_tip,
    wl1_c_angle_tip, n_wl1; r=r_wl1)
  addchord(wingletR, wl2_x_tip, wl2_y_tip, wl2_z_tip, wl2_c_tip,
    wl2_c_angle_tip, n_wl2; r=r_wl2)

  ## --------- FUSELAGE --------------------
  ### Parameters
  f_l = w_b*13/8                   # Length
  f_h = f_l/15                    # Height
  f_w_lpos = 0.2                  # Longitudinal position of wing
  f_w_vpos = 0.6                  # Vertical position of wing
  f_nose_angle = 30*pi/180        # Nose angle
  ### Dimensions
  f_x_root = 0.0
  f_y_root = 0.0
  f_z_root = 0.0
  f_c_root = f_l
  f_c_angle_root = 0.0
  f_y_tip = f_h/2
  f_x_tip = f_y_tip/tan(f_nose_angle)
  f_z_tip = 0.0
  f_c_tip = f_l - f_x_tip
  f_c_angle_tip = 0.0
  ### Creates the lifting surfaces
  fuselage = Wing(f_x_tip, -f_y_tip, f_z_tip, f_c_tip, f_c_angle_tip)
  addchord(fuselage, f_x_root, f_y_root, f_z_root, f_c_root,
          f_c_angle_root, n_f; r=r_f)
  addchord(fuselage, f_x_tip, f_y_tip, f_z_tip, f_c_tip,
          f_c_angle_tip, n_f; r=1/r_f)
  setcoordsystem(fuselage, [0.0,0,0], [1.0 0 0; 0 0 1; 0 -1 0])

  ## --------- VERTICAL TAIL --------------------
  ### Parameters
  vt_b = f_h*7/2
  vt_lambda = 45*pi/180
  vt_c_root = vt_b/2*0.85
  vt_c_tip = vt_c_root/2
  ### Dimensions
  vt_x_root = 0.0
  vt_y_root = 0.0
  vt_z_root = 0.0
  vt_c_angle_root = 0.0
  vt_y_tip = vt_b/2
  vt_x_tip = vt_y_tip*tan(vt_lambda)
  vt_z_tip = 0.0
  vt_c_angle_tip = 0.0
  ### Creates the lifting surfaces
  vtail = Wing(vt_x_root, vt_y_root, vt_z_root, vt_c_root, vt_c_angle_root)
  addchord(vtail, vt_x_tip, vt_y_tip, vt_z_tip, vt_c_tip,
          vt_c_angle_tip, n_vt; r=1/r_vt)
  setcoordsystem(vtail, [f_l-vt_c_root,0,f_h/2], [1.0 0 0; 0 0 1; 0 -1 0])

  ## --------- HORIZONTAL TAIL --------------------
  ### Parameters
  ht_b = w_b*5/16
  ht_ar = 6.0
  ht_tr = 0.75
  ht_twist_root= 0.0
  ht_twist_tip = -2.5
  ht_lambda = 15.0
  ht_gamma = 0.0
  ht_vt_pos = 0.6            # Vertical position of horizontal tail on vertical
  ht_c_root = ht_b/ht_ar/ht_tr
  ### Creates the lifting surface
  htail = simpleWing(ht_b, ht_ar, ht_tr,
                      ht_twist_root, ht_lambda, ht_gamma;
                      twist_tip=ht_twist_tip,
                      n=n_ht, r=r_ht, central=false, refinement=[])
  c_vt_at_ht = vt_c_root + (vt_c_tip-vt_c_root)*ht_vt_pos
  x_vt_at_ht = vt_x_root + (vt_x_tip-vt_x_root)*ht_vt_pos
  setcoordsystem(htail,
    [f_l-vt_c_root+x_vt_at_ht+c_vt_at_ht-ht_c_root, 0 , f_h/2+vt_b/2*ht_vt_pos],
         [1.0 0 0; 0 1 0; 0 0 1])


  # CREATES THE WINGSYSTEM
  ## Wing
  wing = WingSystem()
  addwing(wing, "Inner", wing_inner)
  addwing(wing, "OuterL", wing_outerL)
  addwing(wing, "OuterR", wing_outerR)
  addwing(wing, "AileronL", aileronL)
  addwing(wing, "AileronR", aileronR)
  addwing(wing, "WingletL", wingletL)
  addwing(wing, "WingletR", wingletR)
  ## Body
  body = WingSystem()
  addwing(body, "Fuselage", fuselage)
  addwing(body, "VerticalTail", vtail)
  addwing(body, "HorizontalTail", htail)
  b_O = [-f_l*f_w_lpos, 0, -f_h*(f_w_vpos-0.5)]
  setcoordsystem(body, b_O, [1.0 0 0; 0 1 0; 0 0 1])
  ## Entire system
  system = WingSystem()
  addwing(system, "Wing", wing)
  addwing(system, "Body", body)

  # SOLVES
  if solveVLM
    solve(system, Vinf)
    calculate_field(system, "Ftot"; rhoinf=rhoinf)
    calculate_field(system, "CFtot")
    calculate_field(system, "Mtot")
    calculate_field(system, "CMtot"; qinf=qinf)
  else
    setVinf(system, Vinf)
  end

  # GENERATES FLUID DOMAIN
  P_max = [f_l*5/4, w_b*5/4, (f_h+vt_b/2)*5/4]
  fdom = PP.FluidDomain([0.0, 0.0, 0.0], P_max, [10, 5, 1]*2^3)
  PP.setcoordsystem(fdom, b_O + P_max.*[-1/60, -1/2, -1/6],
                      [1.0 0 0; 0 1 0; 0 0 1])
  fdom_V(X) = Vind(system, X) + Vinf(X,0)


  # GENERATES VTKS
  # run(`rm *.vtk -f`)
  save(system, "test"; save_horseshoes=false, num=num)
  PP.save(fdom, "test"; num=(num==nothing?-1:num))

  # OPENS PARAVIEW
  if openParaview
    fls = ""
    for fl in ["test_Wing_Inner_vlm.vtk", "test_Wing_OuterR_vlm.vtk",
                "test_Wing_OuterL_vlm.vtk",
                "test_Wing_AileronL_vlm.vtk", "test_Wing_AileronR_vlm.vtk",
                "test_Wing_WingletL_vlm.vtk", "test_Wing_WingletR_vlm.vtk",
                "test_Body_Fuselage_vlm.vtk", "test_Body_VerticalTail_vlm.vtk",
                "test_Body_HorizontalTail_vlm.vtk",
                "test_fdom.vtk"]
      fls = fls*fl*";"
    end
    run(`paraview --data=$fls`)
  end


  # # ------- Moves the system's coordinate system
  # displacement = [10, 20, 30]+0.0
  # roll = 30*pi/180
  # pitch = 4*pi/180
  # yaw = 10*pi/180
  #
  # # Transformation matrices
  # Mr = [1 0 0; 0 cos(roll) -sin(roll); 0 sin(roll) cos(roll)]
  # Mp = [cos(pitch) 0 -sin(pitch); 0 1 0; sin(pitch) 0 cos(pitch)]
  # My = [cos(yaw) -sin(yaw) 0; sin(yaw) cos(yaw) 0; 0 0 1]
  # M = Mp*Mr*My
  #
  # setcoordsystem(system, displacement, M)

  return system, fdom, fdom_V
end


"""
    `simpleWing(b, ar, tr, twist, lambda, gamma; twist_tip=twist, n=20, r=2.0)`
  Generates a single-section wing.

  # Arguments
  *   b       :   (float) span.
  *   ar      :   (float) aspect ratio defined as the span over the tip's chord.
  *   tr      :   (float) taper ratio (tip chord / root chord).
  *   twist   :   (float) twist of the root in degrees.
  *   lambda  :   (float) sweep in degrees.
  *   gamma   :   (float) dihedral in degrees.
  *   (OPTIONALS)
  *   twist_tip : (float) twist of the tip if different than root.
  *   n       :   (int)   number of horseshoes per side of the wing.
  *   r       :   (float) horseshoes' expansion ratio
"""
function simpleWing(b::Float64, ar::Float64, tr::Float64,
                    twist::Float64, lambda::Float64, gamma::Float64;
                    twist_tip=nothing,
                    n::Int64=20, r::Float64=2.0, central=false, refinement=[])
  cr = 1/tr
  c_tip = b/ar
  c_root = cr*c_tip
  twist_t = twist_tip==nothing ? twist : twist_tip

  y_tip = b/2
  x_tip = y_tip*tan(lambda*pi/180)
  z_tip = y_tip*tan(gamma*pi/180)

  # Inverts the complex refinement for the opposite wing
  _ref = []
  for i in size(refinement)[1]:-1:1
    push!(_ref, [refinement[i][1], refinement[i][2], 1/refinement[i][3]])
  end

  wing = Wing(x_tip, -y_tip, z_tip, c_tip, twist_t)
  addchord(wing, 0.0, 0.0, 0.0, c_root, twist, n;
              r=r, central=central, refinement=refinement)
  addchord(wing, x_tip, y_tip, z_tip, c_tip, twist_t, n;
              r=central!=false ? r : 1/r,
              central=typeof(central)!=Bool ? 1-central : central,
              refinement=_ref)

  return wing
end

"Saves the wing domain in VTK legacy format"
function save(self::Wing, filename::String;
                  save_horseshoes::Bool=true,
                  path::String="", comment::String="",
                  num=nothing)
  aux = num!=nothing ? ".$num" : ""
  ext = "_vlm"*aux*".vtk"

  if path !=""
    _path = string(path, (path[end]!="/" ? "/" : ""))
  else
    _path = ""
  end

  # HEADER
  n = get_m(self)+1 # Number of chords
  f = open(string(_path, filename, ext), "w")
  header = "# vtk DataFile Version 4.0" # File version and identifier
  header = string(header, "\n", "[MyVLM domain - lattices] ", comment) # Title
  header = string(header, "\n", "ASCII") # File format
  header = string(header, "\n", "DATASET UNSTRUCTURED_GRID")
  write(f, header)

  # POINTS
  nle = n
  nte = n
  ncp = n-1
  if save_horseshoes
    nhs = n-1
  else
    nhs = 0
  end
  write(f, string("\n", "POINTS ", nle+nte+ncp+nhs*6, " float"))
  ## Leading edge
  for i in 1:nle
    LE = getLE(self, i)
    line1 = string(LE[1], " ", LE[2], " ", LE[3])
    write(f, string("\n", line1))
  end
  ## Trailing edge
  for i in 1:nte
    TE = getTE(self, i)
    line2 = string(TE[1], " ", TE[2], " ", TE[3])
    write(f, string("\n", line2))
  end
  ## Control points (m)
  for i in 1:ncp
    CP = getControlPoint(self, i)
    line = string(CP[1], " ", CP[2], " ", CP[3])
    write(f, string("\n", line))
  end
  ## Horseshoes
  x_vor_end = maximum(self._xtwingdcr)*1.25
  for i in 1:nhs
    Ap, A, B, Bp, CP, infDA, infDB, Gamma  = getHorseshoe(self, i)

    # Calculates how long to extend the semi-infinite vortices
    # factor = (x_vor_end - Ap[1])/infDA[1]
    factor = (x_vor_end - self._xtwingdcr[i])/dot(infDA, self.Oaxis[1,:])
    if factor==Inf
      warn("Infinite vortex sheet avoided in visualization")
      factor = (x_vor_end - self._xtwingdcr[i]) / (x_vor_end - self._xtwingdcr[1])/100
    end
    Apinf = Ap + factor*infDA
    # factor = (x_vor_end - Bp[1])/infDB[1]
    factor = (x_vor_end - self._xtwingdcr[i+1])/dot(infDB, self.Oaxis[1,:])
    if factor==Inf
      warn("Infinite vortex sheet avoided in visualization")
      factor = (x_vor_end - self._xtwingdcr[i]) / (x_vor_end - self._xtwingdcr[1])/100
    end
    Bpinf = Bp + factor*infDB

    for point in [Apinf, Ap, A, B, Bp, Bpinf]
      line = string(point[1], " ", point[2], " ", point[3])
      write(f, string("\n", line))
    end
  end

  # CELLS
  write(f, string("\n\n", "CELLS ", n-1+ncp+nhs,
                  " ", (n-1)*5 + ncp*2 + nhs*7))
  ## Lattices
  for i in 0:n-2
    line = string(4, " ", i, " ", i+nle, " ", i+nle+1, " ", i+1)
    write(f, string("\n", line))
  end
  ## Control points
  for i in nle+nte:nle+nte+ncp-1
    line = string(1, " ", i)
    write(f, string("\n", line))
  end
  ## Horseshoes
  aux1 = nle+nte+ncp+nhs*6
  count = 1
  line = string(6, " ")
  for i in nle+nte+ncp:aux1-1
    line = string(line, " ", i)
    if count%6==0
      write(f, string("\n", line))
      line = string(6, " ")
    end
    count += 1
  end


  # CELL TYPES
  write(f, string("\n\n", "CELL_TYPES ", n-1+ncp+nhs))
  ## Lattices
  for i in 0:n-2
    write(f, string("\n", 9))
  end
  ## Control points
  for i in 1:ncp
    write(f, string("\n", 1))
  end
  ## Horseshoes
  for i in 1:nhs
    write(f, string("\n", 4))
  end


  # FIELDS
  initiated = false
  for field_name in keys(self.sol)
    if false==(field_name in keys(FIELDS))
      error(string("CRITICAL ERROR: field ", field_name, " not found in FIELDS"))
    end

    if initiated==false
      write(f, string("\n\n", "CELL_DATA ", n-1+ncp+nhs))
      initiated = true
    end

    if FIELDS[field_name][2]=="vector"
      write(f, string("\n\n", "VECTORS ", field_name," float"))
      for i in 1:n-1+ncp+nhs
            vect = self.sol[field_name][(i-1)%ncp+1]
            line = string(vect[1], " ", vect[2], " ", vect[3])
            write(f, string("\n", line))
      end
    elseif FIELDS[field_name][2]=="scalar"
      write(f, string("\n\n", "SCALARS ", field_name," float"))
      write(f, string("\n", "LOOKUP_TABLE default"))
      for i in 1:n-1+ncp+nhs
            sclr = self.sol[field_name][(i-1)%ncp+1]
            try
              line = string(isnan(sclr) ? -1 : sclr)
            catch
              error(string(sclr))
            end

            write(f, string("\n", line))
      end
    else
      error(string("CRITICAL ERROR: field type ", FIELDS[field_name][2],
            " hasn't been implemented!"))
    end
  end

  close(f)
end

function save(self::WingSystem, filename::String;
                    save_horseshoes::Bool=true,
                    path::String="", comment::String="",
                    num=nothing)
  for (i, wing) in enumerate(self.wings)
    save(wing, "$(filename)_$(self.wing_names[i])",
                save_horseshoes=save_horseshoes,
                path=path, comment=comment, num=num)
  end
end

"Receives a function generate_wing(n) and plots aerodynamic characteristics
at different lattice resolutions"
function grid_dependance(wing_function, Vinf; ns=[4*2^i for i in 1:6],
                          S="automatic", rhoinf=9.093/10^1,
                          verbose=true,
                          save_w=false, save_name="griddep", save_fd=false,
                          fig_title="")

  to_plot = ["CD", "CL", "Gamma"]
  to_print = ["CD", "CL", "Gamma"]
  vals = Dict([ (key, []) for key in to_print])
  wing = nothing
  for (i,n) in enumerate(ns)
    verbose ? println("Starting n=$n ...") : nothing

    wing = wing_function(n)

    solve(wing, Vinf)
    calculate_field(wing, "CFtot"; rhoinf=rhoinf)

    info = fields_summary(wing)

    str = "Iteration #$i"
    for key in keys(info)
      if key in to_plot || key in to_print
        push!(vals[key], info[key])
        str *="\n\t $(key): $(info[key])"
      end
    end
    verbose ?  println("\t Done!") : nothing

    if save_w
      save(wing, save_name; save_horseshoes=false, num=i)
      if save_fd
        # # fd = PP.FluidDomain(wing, [-0.25, -0.75, -0.125], [2.0, 0.75, 0.25], [8, 8, 5]*4);
        # fd = PP.FluidDomain(wing, [-0.25, -1.75, -0.125], [2.5, 1.75, 0.5], [8, 8, 5]*4);
        # PP.solve(fd, "U"; Vinf=Vinf);
        # # PP.solve(fd, "Uind"; Vinf=Vinf);
        # PP.solve(fd, "Cp"; Vinf=Vinf);
        # PP.savedomain(fd, save_name*"_fd")
      end
    end
  end


  n_vals = length(to_plot)
  fig = figure(fig_title, figsize=(7*n_vals,5))
  suptitle(fig_title, fontsize="x-large")
  for (i,key) in enumerate(to_plot)
    subplot(100 + 10*n_vals + i)
    title("$key dependance to grid")
    y = vals[key]
    plot(ns, y, "o")
    xlim([minimum(ns),maximum(ns)*1.25]);
    xlabel("Grid");
    y_min = minimum(y)
    y_max = maximum(y)
    y_low = y_min - (y_max-y_min)*0.25
    y_up = y_max + (y_max-y_min)*0.25
    ylim([y_low, y_up]);
    ylabel(key);
  end

  println("VALUES AT n=$(ns[end])")
  for key in to_print
    println("\t$key : $(vals[key][end])")
  end
  return wing, ns, vals
end


################################################################################
# CALCULATIONS
################################################################################
# """
#   `V(self::Wing or WingSystem, P::Float64[])`
# Calculates the induced velocity at point P (Gamma field must have been solved)
# """
# function V(self, P; ign_col::Bool=false)
#   if false==("Gamma" in keys(self.sol))
#     error("ERROR: Gamma field not found!")
#   end
#
#   m = get_m(self)
#
#   V_tot = zeros(3)
#   # Iterates over every horseshoe in the wing
#   for i in 1:m
#     this_HS = getHorseshoe(self, i)
#     V_this = VLMSolver.V(this_HS, P; ign_col=ign_col)
#     V_tot += V_this
#   end
#   return V_tot
# end



################################################################################
# ALGEBRA
################################################################################
"""
Rotates and translates the vector V.

Receives the i', j', k' unit vectors of an euclidean system with origin T, and
returns V'. (In this version, the unit vectors have been organized as a matrix
M)
"""
function transform(V::typeof(Float64[]),
                    M::Array{Float64,2}, T::typeof(Float64[]))
  return M*(V-T)
end

function transform(Vs::Array{Array{Float64,1},1},
                    M::Array{Float64,2}, T::typeof(Float64[]))
  out = Array{Float64,1}[]
  for V in Vs
    push!(out, transform(V, M, T))
  end
  return out
end

"""
Rotates and translates back a vector V' that had been rotated and translated
into the system (i', j', k') with origin T, and returns the original V.
To ease repetitive computation, instead of giving the unit vectors, give the
inverse of their matrix.
"""
function countertransform(Vp::typeof(Float64[]),
                          invM::Array{Float64,2}, T::typeof(Float64[]))
  return invM*Vp + T
end

function countertransform(Vps::Array{Array{Float64,1},1},
                          invM::Array{Float64,2}, T::typeof(Float64[]))
  out = Array{Float64,1}[]
  for Vp in Vps
    push!(out, countertransform(Vp, invM, T))
  end
  return out
end

"Checks that the unit vectors given as the matrix M=[i;j;k] define a coordinate
system"
function check_coord_sys(M::Array{Float64,2}; raise_error::Bool=true)
  # Checks normalization
  for i in 1:size(M)[1]
    if abs(norm(M[i,:])-1) > 0.00000001
      println(M)
      if raise_error
        error("Not unitary axis: $(M[i,:])")
      else
        return false
      end
    end
  end

  # Checks ortogonality
  for i in size(M)[1]
    xi = M[i, :]
    xip1 = M[(i%size(M)[1])+1, :]
    proj = abs(dot(xi, xip1))
    if proj>0.00000001
      if raise_error
        error("Non-ortogonal system $M")
      else
        return false
      end
    end
  end
  return true
end

function check_coord_sys(M::Array{Array{Float64,1},1}; raise_error::Bool=true)
  dims = 3
  newM = zeros(dims,dims)
  for i in 1:dims
    newM[i, :] = M[i]
  end
  return check_coord_sys(newM; raise_error=raise_error)
end

"returns the transformation matrix of rotation around an arbitrary axis of unit
vector `r`"
function axis_rotation(r, angle_deg)
  ux, uy, uz = r
  C = cos(angle_deg*pi/180)
  S = sin(angle_deg*pi/180)
  t = 1 - C
  M = [t*ux^2+C t*ux*uy-S*uz t*ux*uz+S*uy;
        t*ux*uy+S*uz t*uy^2+C t*uy*uz-S*ux;
        t*ux*uz-S*uy t*uy*uz+S*ux t*uz^2+C]
  return M
end

##### END OF ALGEBRA ###########################################################
