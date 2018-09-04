
# APC Thin Electric 10x7 propeller as described in McCrink, M. H., & Gregory, J.
# W. (2017), *Blade Element Momentum Modeling of Low-Reynolds Electric
# Propulsion Systems*

scale = 5.25                     # Scaling factor

Rtip = 10*0.0254/2*scale      # (m) Radius of blade tip
Rhub = 0.05*Rtip              # (m) Radius of hub
# Rhub = 1.2/100              # (m) Radius of hub
B = 2                         # Number of blades

data_sec = CSV.read(joinpath(data_path, "apc10x7_sections.csv"))
data_chord = CSV.read(joinpath(data_path, "apc10x7_chord.csv"),
                                                header=["r","c"], datarow=1)
data_twist = CSV.read(joinpath(data_path, "apc10x7_twist.csv"),
                                                header=["r","theta"], datarow=1)

# r/R y/R (y-distance of LE from the middle point of hub)
sweepdist = [data_sec[i,j] for i=size(data_sec,1):-1:1, j=[2,5]]*0.0254/(Rtip/scale)

# r/R z/R  (height of leading edge from top face of hub)
heightdist = [data_sec[i,j] for i=size(data_sec,1):-1:1, j=[2,6]]*0.0254/(Rtip/scale)

# r/R c/R
chorddist = [data_chord[i,j] for i=1:size(data_chord,1), j=[1,2]]/(Rtip/scale)

# r/R twist
pitchdist = [data_twist[i,j]/(j==1 ? (Rtip/scale) : 1) for i=1:size(data_twist,1), j=[1,2]]

n4412_file = "n4412-1500000.csv"
x,y = vlm.ap.readcontour(data_path*"airfoils/naca4412.dat")
airfoiln4412 = hcat(x,y)

# Airfoils along the blade as
# airfoil_contours=[ (pos1, contour1), (pos2, contour2), ...] with contour=(x,y)
# and pos the position from root to tip between 0 and 1. pos1 must equal 0
# (root airfoil) and the last must be 1 (tip airfoil)
airfoil_contours = [ (0, airfoiln4412, n4412_file),
                      (0.125, airfoiln4412, n4412_file),
                      (0.25, airfoiln4412, n4412_file),
                      (0.375, airfoiln4412, n4412_file),
                      (0.50, airfoiln4412, n4412_file),
                      (0.625, airfoiln4412, n4412_file),
                      (0.75, airfoiln4412, n4412_file),
                      (0.875, airfoiln4412, n4412_file),
                      (1, airfoiln4412, n4412_file)]
