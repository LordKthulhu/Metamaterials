#Installing relevant packages
using Pkg
Pkg.add(["Plots","Printf","DelimitedFiles","Glob","ProgressMeter","LaTeXStrings","Reinforce","Shell","PyCall","PyPlot","Dates","JLD","Colors","Conda"])
Pkg.build("PyCall")
using Conda
Conda.add("matplotlib")
