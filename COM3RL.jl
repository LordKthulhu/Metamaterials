using Plots
using Printf
using DelimitedFiles
using Glob
using ProgressMeter
using LaTeXStrings

include("Ressources/structures.jl")
include("Ressources/simulationTools.jl")
include("Ressources/COM3RWTools.jl")

const H = 3
const unitSize = 1
const dEpsilon = 2e-4
const area = 12*(H-1)*unitSize*1.0

include("Ressources/geometryTools.jl")

ENV["GKSwstype"] = "100"

simulations = []
