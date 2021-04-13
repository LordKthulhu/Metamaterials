using Plots
using Printf
using DelimitedFiles
using Statistics
using Glob
using Dates
using ProgressMeter

include("Ressources/structures.jl")
include("Ressources/simulationTools.jl")
include("Ressources/COM3RWTools.jl")

const H,iterations,randomMat = parseArguments()
const unitSize = 1
const dEpsilon = 2e-4
const area = 12*(H-1)*unitSize*1.0

include("Ressources/geometryTools.jl")

ENV["GKSwstype"] = "100"
barPlots = [[] for i in 1:iterations]

simulations = [ emptySimulation(iter,dEpsilon) for iter in 1:iterations ]
simulationsPlain = [ emptySimulation("plain$iter",dEpsilon) for iter in 1:iterations ]

@printf("Starting %d simulations on %d threads.\n",iterations,Threads.nthreads())

################################################################################
#####                            SIMULATIONS                               #####
################################################################################

progress = pBar(iterations,"Computing progess... ")

Threads.@threads for iter = 1:iterations

    while simulations[iter].exit == 1

        skeleton = randomSkeleton(H)
        plotGeometry(barPlots,iter,skeleton) # Adding data for geometry plot

        if randomMat
            tensile = round(3+2*rand(),digits=2)
            compressive = round(30+20*rand(),digits=2)
            material = Material(true,compressive,tensile)
            materialPlain = Material(false,compressive,tensile)
        else
            material = Material(true,45.0,4.8)
            materialPlain = Material(false,45.0,4.8)
        end

        model = modelFromSkeleton(skeleton, material, unitSize, nodeWeights, linkWeights)
        modelPlain = modelFromSkeleton(skeleton, materialPlain, unitSize, nodeWeights, linkWeights)

        simulations[iter].model = model
        simulationsPlain[iter].model = modelPlain

        runSimulation(simulations[iter])
        runSimulation(simulationsPlain[iter])
    end

    filename = simulations[iter].filename
    run(`rm $filename/$filename-MECH.crk $filename/$filename-MECH.fld
            $filename/$filename-MECH.int $filename/$filename-MECH.tmp $filename-restart.aux`)
    io = open(filename*"/"*filename*"-results.csv","a")
    writedlm(io,transpose(simulations[iter].strain),",")
    writedlm(io,transpose(simulations[iter].stress),",")
    close(io)
    filename = simulationsPlain[iter].filename
    run(`rm $filename/$filename-MECH.crk $filename/$filename-MECH.fld
            $filename/$filename-MECH.int $filename/$filename-MECH.tmp $filename-restart.aux`)
    io = open(filename*"/"*filename*"-results.csv","a")
    writedlm(io,transpose(simulations[iter].strain),",")
    writedlm(io,transpose(simulations[iter].stress),",")
    close(io)
    next!(progress)
end

################################################################################
#####                           RESULTS OUTPUT                             #####
################################################################################

# Geometry plots

progress = pBar(2*iterations,"Plotting...          ",dt=0.5)

for i=1:iterations
    plt = plot(showaxis=false,size=(400,400))
    for line in barPlots[i]
        plot!(plt, line[1], line[2], lw =3, lc = :black, label = false)
    end
    png(plt,"metamat"*string(i)*"/metamat"*string(i)*"-barplot.png")
    next!(progress)
    sleep(0.1)
end

# Stress-strain plots

for i=1:iterations
    plotfile = "metamat$i/metamat$i-plot.png"
    plt = plot(simulations[i].strain,simulations[i].stress, lw = 3, ylabel = "Overall Stress (kg/cm2)",lc="green", label="PVA-ECC")
    plot!(plt,simulationsPlain[i].strain,simulationsPlain[i].stress, lw = 3,lc="red", label="Plain Concrete")
    png(plt,plotfile)
    next!(progress)
    sleep(0.1)
end

# CSV output

weights = [ simulation.model.weight for simulation in simulations ]
maxStrains = [ maximum(simulation.strain) for simulation in simulations ]
maxStresses = [ maximum(simulation.stress) for simulation in simulations ]
energyAbsorptions = [ energy(simulation.strain,simulation.stress) for simulation in simulations ]
maxStrainsPlain = [ maximum(simulation.strain) for simulation in simulationsPlain ]
maxStressesPlain = [ maximum(simulation.stress) for simulation in simulationsPlain ]
energyAbsorptionsPlain = [ energy(simulation.strain,simulation.stress) for simulation in simulationsPlain ]

io = open("results.csv","a")
writedlm(io,transpose(weights),",")
writedlm(io,transpose(maxStrains),",")
writedlm(io,transpose(maxStresses),",")
writedlm(io,transpose(energyAbsorptions),",")
writedlm(io,transpose(maxStrainsPlain),",")
writedlm(io,transpose(maxStressesPlain),",")
writedlm(io,transpose(energyAbsorptionsPlain),",")
close(io)

# General plots over all simulations

strainsPlt = plot(weights,maxStrains, seriestype = :scatter, xlabel = "Area (cm2)",ylabel = "Failure strain",label="PVA-ECC", color = :blue)
plot!(strainsPlt, weights,maxStrainsPlain, seriestype = :scatter, label="Plain Concrete", color = :red)
png(strainsPlt,"strains.png")

sleep(0.1)
stressPlt = plot(weights,maxStresses, seriestype = :scatter, xlabel = "Area (cm2)",ylabel = "Failure stress",label="PVA-ECC", color = :blue)
plot!(stressPlt, weights,maxStressesPlain, seriestype = :scatter, label="Plain Concrete", color = :red)
png(stressPlt,"stresses.png")

sleep(0.1)
energyPlt = plot(weights,energyAbsorptions, seriestype = :scatter, xlabel = "Area (cm2)",ylabel = "Absorbed energy at failure",label="PVA-ECC", color = :blue)
plot!(energyPlt, weights, energyAbsorptionsPlain, seriestype = :scatter, label="Plain Concrete", color = :red)
png(energyPlt,"energies.png")
