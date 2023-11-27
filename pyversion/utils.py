import getopt

def usage():
    print(" -h or --help : This info")
    print(" -v or --verbose")
    print(" -n or --nodes Nodes_per_dimension (int) ")
    print(" -s or --separation Grid_separation (float)")
    print(" -m or --mass Mass_of_node (float)")
    print(" -f or --fcon Force_constant (float)")
    print(" -i or --interact Node_interaction_level (int)")
    print(" -g or --gravity Gravity (float)")
    print(" -b or --ballsize Radius_of_ball (float)")
    print(" -o or --offset offset_of_falling_cloth (float)")
    print(" -t or --timestep timestep (float)")
    print(" -u or --update Timesteps_per_display_update (int)")
    return

def read_arg(argv):
    verbose, dt, N, mass, fcon, separation, ballsize, gravity, offset, interact, update = (
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    )
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "hvn:s:m:f:i:g:t:o:u:g:b:d:",
            [
                "help",
                "verbose",
                "nodes=",
                "separation=",
                "mass=",
                "fcon=",
                "interact=",
                "gravity=",
                "ballsize=",
                "offset=",
                "timestep=",
                "update=",
            ],
        )
    except getopt.GetoptError:
        print("using default parameters")
        opts = {}
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-v", "--verbose"):
            verbose = 1
        elif o in ("-n", "--nodes"):
            N = int(a)
        elif o in ("-s", "--separation"):
            separation = float(a)
        elif o in ("-m", "--mass"):
            mass = float(a)
        elif o in ("-f", "--fcon"):
            fcon = float(a)
        elif o in ("-i", "--interact"):
            interact = int(a)
        elif o in ("-g", "--gravity"):
            gravity = float(a)
        elif o in ("-b", "--ballsize"):
            ballsize = float(a)
        elif o in ("-o", "--offset"):
            offset = float(a)
        elif o in ("-t", "--timestep"):
            dt = float(a)
        elif o in ("-u", "--update"):
            update = int(a)
        else:
            assert False, "unhandled option"

    print("The cloth ")
    print("  Nodes per dimension ", N)
    print("  Grid Separation     ", separation)
    print("  Mass of node        ", mass)
    print("  Force constant      ", fcon)
    print("  Node interaction    ", interact)
    print("The Environment")
    print("  Gravity             ", gravity)
    print("  Ballsize            ", ballsize)
    print("  Offset              ", offset)
    print("The Simulation")
    print("  Timestep            ", dt)
    print("  Updates per display ", update)
    print("  Verbose             ", verbose)
    return verbose, dt, N, mass, fcon, separation, ballsize, gravity, offset, interact, update