# Add packages. This only has to be done once under the environment of Disorderness.jl
# After that, before every time running the project you just simply:
# 	Pkg.activate(rootdir)
# 	Pkg.instantiate(); # install the packages in the same state that is given by that manifest.
# and you are freely to `using` the packages in dependency.
using Pkg;
Pkg.add(["Plots","PlutoUI","DifferentialEquations","Formatting","LaTeXStrings"])
