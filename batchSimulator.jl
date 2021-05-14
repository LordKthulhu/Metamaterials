push!(LOAD_PATH,pwd()*"/Ressources")
using MetamatModule
using DelimitedFiles
using Plots
using Printf
using LaTeXStrings
using Shell
using Dates
using Glob

const H,iterations,randomMat,parameters = parseArguments()
const unitSize = 1
const dEpsilon = 2e-4
const nodeWeights = makeNodeWeights(H,unitSize)
const linkWeights = makeLinkWeights(H,unitSize)
const basePoints = makeBasePoints()
const baseElements = makeBaseElements(H)

if Sys.iswindows()
    Shell.run("cleanup.sh")
else
    run(`./cleanup.sh`)
end

ENV["GKSwstype"] = "100"
barPlots = [[] for i in 1:iterations]
skeletonLinks = [ falses(H,H,4) for i in 1:iterations ]

parameterValues = []
repeatedSimulation = 1
materials = []

if parameters != "none"
    parameterValues = readdlm(parameters,',')
    repeatedSimulation = size(parameterValues,1)
    if randomMat
        error("Cannot use random material and parametric study simultaneously.")
    end
end

simulations = [ emptySimulation("$iter-$j",dEpsilon) for iter in 1:iterations,  j in 1:repeatedSimulation ]

@printf("Starting %dx%d simulations on %d threads.\n",iterations,repeatedSimulation,Threads.nthreads())

################################################################################
#####                            SIMULATIONS                               #####
################################################################################

progress = pBar(iterations,"Computing progess... ")

Threads.@threads for iter = 1:iterations

    while (s->s.exit).(simulations[iter,:]) != zeros(repeatedSimulation)

        materials = []
        if parameters != "none"
            materials = [ Material(true,parameterValues[i,1],parameterValues[i,2]) for i in 1:size(parameterValues,1) ]
        elseif randomMat
            tensile = round(3+2*rand(),digits=2)
            compressive = round(30+20*rand(),digits=2)
            materials = [ Material(true,compressive,tensile) ]
        else
            materials = [ Material(true,45.0,4.8) ]
        end

        skeleton = randomSkeleton(H)
        barPlots[iter] = plotGeometry(skeleton) # Adding data for geometry plot

        for i in 1:repeatedSimulation
            model = modelFromSkeleton(skeleton, materials[i], unitSize, basePoints, baseElements, nodeWeights, linkWeights)
            simulations[iter,i].model = model
            simulations[iter,i].step = 0; simulations[iter,i].strain = []; simulations[iter,i].stress = []
            runSimulation(simulations[iter,i])
        end
        skeletonLinks[iter] = skeleton.links
    end

    for simulation in simulations[iter,:]
        filename = simulation.filename
        run(`rm $filename/$filename-MECH.crk $filename/$filename-MECH.fld $filename/$filename-MECH.int $filename/$filename-MECH.tmp $filename-restart.aux`)
        io = open(filename*"/"*filename*"-results.csv","a")
        writedlm(io,transpose(simulation.strain),",")
        writedlm(io,transpose(simulation.stress),",")
        close(io)
    end
    next!(progress)
end

################################################################################
#####                           RESULTS OUTPUT                             #####
################################################################################

# Number of unique simulations

@printf("Number of unique simulations : %d out of %d simulations total.\n",length(unique(skeletonLinks)),iterations)

# Geometry plots

progress = pBar(2*iterations,"Plotting...          ",dt=0.5)

for i=1:iterations
    plt = plot(showaxis=false,size=(400,400),title="Metamaterial $i")
    for line in barPlots[i]
        plot!(plt, line[1], line[2], lw =3, lc = :black, label = false)
    end
    png(plt,"metamat"*string(i)*"-1/metamat"*string(i)*"-barplot.png")
    next!(progress)
    #sleep(0.01)
end

# Stress-strain plots

if randomMat
    labels = [ "Random material properties" ]
elseif repeatedSimulation > 1
    labels = [ L"\sigma_c = %$(parameterValues[i,1]) MPa, \sigma_t = %$(parameterValues[i,2]) MPa" for i in 1:repeatedSimulation ]
else
    labels = [ L"\sigma_c = 45 MPa, \sigma_t = 4.8 MPa" for i in 1:repeatedSimulation ]
end

for iter=1:iterations
    plotfile = "metamat$iter-1/metamat$iter-plot.png"
    plt = plot(lw = 3, xlabel = "Strain", ylabel = "Overall Stress (kg/cm2)", palette = :tab10, ylims = (0,Inf))
    for i in 1:repeatedSimulation
        plot!(plt,simulations[iter,i].strain,simulations[iter,i].stress, label=labels[i])
    end
    png(plt,plotfile)
    next!(progress)
    #sleep(0.01)
end

# CSV output

weights = [ simulation.model.weight for simulation in simulations ]
maxStrains = [ maximum(simulation.strain) for simulation in simulations ]
maxStresses = [ maximum(simulation.stress) for simulation in simulations ]
energyAbsorptions = [ energy(simulation.strain,simulation.stress) for simulation in simulations ]

io = open("results.csv","a")
writedlm(io,transpose(weights),",")
writedlm(io,transpose(maxStrains),",")
writedlm(io,transpose(maxStresses),",")
writedlm(io,transpose(energyAbsorptions),",")
close(io)

# General plots over all simulations

strainsPlt = plot(xlabel = "Area (cm2)",ylabel = "Failure strain", palette = :tab10)
for i in 1:repeatedSimulation
    plot!(strainsPlt, weights[:,i], maxStrains[:,i], seriestype = :scatter, label=labels[i])
end
png(strainsPlt,"strains.png")

sleep(0.1)
stressPlt = plot(xlabel = "Area (cm2)",ylabel = "Failure stress", palette = :tab10)
for i in 1:repeatedSimulation
     plot!(stressPlt,weights[:,i],maxStresses[:,i], seriestype = :scatter, label=labels[i])
end
png(stressPlt,"stresses.png")

sleep(0.1)
energyPlt = plot(xlabel = "Area (cm2)",ylabel = "Absorbed energy at failure", palette = :tab10)
for i in 1:repeatedSimulation
    plot!(energyPlt,weights[:,i],energyAbsorptions[:,i], seriestype = :scatter,label=labels[i])
end
png(energyPlt,"energies.png")

currentTime = Dates.format(now(),"dd-mm-yyyy_HH:MM:SS")
run(`mkdir Batch_$currentTime`)
folders = glob("metamat*/")
plots = glob("*.png")
run(`mv $folders $plots results.csv Batch_$currentTime/`)
