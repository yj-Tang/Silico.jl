using Pkg: Pkg
Pkg.activate("$(@__DIR__)")

using Silico: Visualizer, set_floor!, set_light!, set_background!, Diagonal, get_3d_polytope_drop, Mehrotra, simulate!, build_mechanism!, set_mechanism!, visualize!, normalize, planar_sliding