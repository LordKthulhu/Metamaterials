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
include("Ressources/geometryTools.jl")

ENV["GKSwstype"] = "100"



function pBar(len,text; dt=1)
    Progress(len, dt = dt, desc = text , barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',), barlen=20)
end

#### MAIN ####

const H,iterations,randomMat = parseArguments()

@printf("Starting %d simulations on %d threads.\n",iterations,Threads.nthreads())

global maxStrains = zeros(iterations)
global MaxStresses = zeros(iterations)
global energyAbsorptions = zeros(iterations)
global weights = zeros(iterations)

global simulations = zeros(iterations)

global plt2 = [plot() for i in 1:iterations]

global barPlots = [[] for i in 1:iterations]

const unitSize = 1
const dEpsilon = 2e-4
const area = 12*(H-1)*unitSize*1.0

nodeWeights = [4*ones(1,H);8*ones(H-2,H);4*ones(1,H)]

nodeWeights[:,1] = 4*ones(H); nodeWeights[:,end] = 4*ones(H)
nodeWeights[1,1] = 2; nodeWeights[1,end] = 2; nodeWeights[end,1] = 2; nodeWeights[end,end] = 2
nodeWeights = unitSize^2 .* nodeWeights

linkWeights = zeros(H,H,4)
linkWeights[:,:,1] = (20*sqrt(2)-2)*unitSize^2 .* ones(H,H); linkWeights[:,:,3] = (20*sqrt(2)-2)*unitSize^2 .* ones(H,H)
linkWeights[:,:,2] = (20-4*sqrt(2))*unitSize^2 .* ones(H,H); linkWeights[:,:,4] = (20-4*sqrt(2))*unitSize^2 .* ones(H,H)

basePointsData = readdlm("Ressources/points.csv",',')
basePoints = [ Point(Int(basePointsData[i,1]),[basePointsData[i,2],0.0,basePointsData[i,3]]) for i in 1:length(basePointsData[:,1])]
baseElementsData = readdlm("Ressources/elements.csv",',',Int)

for index in findall(x -> x > 1000, baseElementsData)
    if baseElementsData[index] < 2000
        baseElementsData[index] = baseElementsData[index]%1000 + (H-1)*145
    elseif baseElementsData[index] < 3000
        baseElementsData[index] = baseElementsData[index]%2000 + H*145
    elseif baseElementsData[index] < 4000
        baseElementsData[index] = baseElementsData[index]%3000 + (H+1)*145
    else
        baseElementsData[index] = baseElementsData[index]%4000 + 145
    end
end

baseElements = [ Element(i,baseElementsData[i,2:end]) for i in baseElementsData[:,1] ]

simulations = [ emptySimulation(iter,dEpsilon) for iter in 1:iterations ]

################################################################################
#####                             MAIN LOOP                                #####
################################################################################

progress = pBar(iterations,"Computing progess... ")

Threads.@threads for iter = 1:iterations

    while simulations[iter].exit == 1

        skeleton = randomSkeleton(H)
        plotGeometry(barPlots,iter,skeleton) # Adding data for geometry plot

        if randomMat
            tensile = round(3+2*rand(),digits=2)
            compressive = round(30+20*rand(),digits=2)
            material = Material(compressive,tensile)
        else
            material = Material(45.0,4.8)
        end

        model = modelFromSkeleton(skeleton, material, unitSize, nodeWeights, linkWeights)
        simulations[iter].model = model

        runSimulation(simulations[iter])
    end

    ###### If all COM3 runs terminate, go on to results processing ######

    filename = simulations[iter].filename

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
    sleep(0.5)
end

# Stress-strain plots

for i=1:iterations
    plotfile = "metamat$i/metamat$i-plot.png"
    plt = plot(simulations[i].strain,simulations[i].stress, lw = 3, ylabel = "Overall Stress (kg/cm2)",lc="green", label="PVA-ECC")
    png(plt,plotfile)
    next!(progress)
    sleep(0.5)
end

# CSV output

weights = [ simulation.model.weight for simulation in simulations ]
maxStrains = [ maximum(simulation.strain) for simulation in simulations ]
maxStresses = [ maximum(simulation.stress) for simulation in simulations ]
energyAbsorptions = [ energy(simulation.strain,simulation.stress) for simulation in simulations ]

display(maxStresses)

io = open("results.csv","a")
writedlm(io,transpose(weights),",")
writedlm(io,transpose(maxStrains),",")
writedlm(io,transpose(MaxStresses),",")
writedlm(io,transpose(energyAbsorptions),",")
close(io)

# General plots over all simulations

strainsPlt = plot(weights,maxStrains, seriestype = :scatter, xlabel = "Area (cm2)",ylabel = "Failure strain",label="PVA-ECC", color = :blue)
png(strainsPlt,"strains.png")

sleep(0.5)
stressPlt = plot(weights,MaxStresses, seriestype = :scatter, xlabel = "Area (cm2)",ylabel = "Failure stress",label="PVA-ECC", color = :blue)
png(stressPlt,"stresses.png")

sleep(0.5)
energyPlt = plot(weights,energyAbsorptions, seriestype = :scatter, xlabel = "Area (cm2)",ylabel = "Absorbed energy at failure",label="PVA-ECC", color = :blue)
png(energyPlt,"energies.png")
