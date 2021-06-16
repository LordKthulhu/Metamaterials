module MetamatModule

using Plots
using Printf
using DelimitedFiles
using Glob
using ProgressMeter
using LaTeXStrings
using Shell
using Statistics
using Reinforce: AbstractEnvironment


import ProgressMeter: next!
import Reinforce: reset!, actions, finished, step!, state, maxsteps
import Base: copy

export
    next!,
    Point,
    Element,
    Material, randonMaterial,
    Skeleton, randomSkeleton, alteredSkeleton, copy,
    Model, modelFromSkeleton,
    Simulation, emptySimulation, runSimulation,
    makeNodeWeights, makeLinkWeights, makeBasePoints, makeBaseElements,
    randomGeometry, plotGeometry,
    pBar, execute, parseArguments,
    runSteps, energy,
    MetamatEnv, maxsteps, reset!, step!, actions, finished,
    MetamatEnv2,
    fullScaleModelFromModel, plotModel

include("structures.jl")
include("simulationTools.jl")
include("COM3RWTools.jl")
include("geometryTools.jl")
include("MetamatRLEnv.jl")
include("MetamatRLEnv2.jl")
end
