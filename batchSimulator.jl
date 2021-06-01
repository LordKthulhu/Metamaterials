push!(LOAD_PATH,pwd()*"/Ressources")
using MetamatModule
using DelimitedFiles
using Plots
using Printf
using LaTeXStrings
using Shell
using Dates
using Glob
using PyCall
using PyPlot
using JLD

const H,iterations,randomMat,parameters = parseArguments()
const unitSize = 1
const dEpsilon = 4e-4
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
            materials = [ Material(parameterValues[i,1],parameterValues[i,2],parameterValues[i,3],parameterValues[i,4],parameterValues[i,5],parameterValues[i,6]) for i in 1:size(parameterValues,1) ]
        elseif randomMat
            compressive = round(30+20*rand(),digits=2)
            tensile = round(3+2*rand(),digits=2)
            tensilePeak = round(0.002+0.005*rand(),digits=4)
            peakStrain = round(0.01+0.05*rand(),digits=3)
            crackStrainRatio = round(0.05+0.045*rand(),digits=3)
            materials = [ Material(compressive,tensile,tensilePeak,peakStrain,crackStrainRatio) ]
        else
            materials = [ Material(45.0,4.8) ]
        end

        skeleton = randomSkeleton(H)
        barPlots[iter] = plotGeometry(skeleton) # Adding data for geometry plot

        for i in 1:repeatedSimulation
            model = modelFromSkeleton(skeleton, materials[i], unitSize, basePoints, baseElements, nodeWeights, linkWeights)
            #model = fullScaleModelFromModel(fullScaleModelFromModel(model))
            simulations[iter,i].model = model
            simulations[iter,i].step = 0; simulations[iter,i].strain = []; simulations[iter,i].stress = []
            runSimulation(simulations[iter,i])
        end
        skeletonLinks[iter] = skeleton.links
    end

    for simulation in simulations[iter,:]
        filename = simulation.filename
        run(`rm $filename/$filename-MECH.int $filename/$filename-MECH.tmp $filename-restart.aux`)
        io = open(filename*"/"*filename*"-results.csv","a")
        writedlm(io,transpose(simulation.strain),",")
        writedlm(io,transpose(simulation.stress),",")
        close(io)
    end
    next!(progress)
end

GC.gc(true)

################################################################################
#####                           RESULTS OUTPUT                             #####
################################################################################

isStrainHardening = [ false for i in 1:iterations ]
for iter in 1:iterations
    stress = simulations[iter,1].stress
    maxi = 0; current = stress[1]
    index = 1
    while index <= length(stress)-1 && current > 0.9 * maxi
        current > maxi ? maxi = current : nothing
        index += 1
        current = stress[index]
    end
    mini = current
    while index <= length(stress)-1
        current < mini ? mini = current : nothing
        if current > 1.15 * mini
            isStrainHardening[iter] = true
            break
        end
        index += 1
        current = stress[index]
    end
end

io=open("StrainHard.csv","w")
writedlm(io,isStrainHardening,',')
close(io)

# Number of unique simulations

@printf("Number of unique simulations : %d out of %d simulations total.\n",length(unique(skeletonLinks)),iterations)

# Geometry plots

progress = pBar(2*iterations,"Plotting...          ",dt=0.5)

for i=1:iterations
    # fig = PyPlot.figure()
    # PyPlot.title("Metamaterial $i")
    # points = simulations[i,1].model.points
    # for element in simulations[i,1].model.elements
    #     P = [ points[findall(p->p.n==element.points[i],points)[1]] for i in 1:4 ]
    #     X = (p->p.x).(P); Y = (p->p.z).(P)
    #     PyPlot.fill(X,Y, facecolor = "gray", edgecolor = "black")
    #     PyPlot.fill(24*(H-1).-X, Y, facecolor = "gray", edgecolor = "black")
    #     PyPlot.fill(X, 24*(H-1).-Y, facecolor = "gray", edgecolor = "black")
    #     PyPlot.fill(24*(H-1).-X, 24*(H-1).-Y, facecolor = "gray", edgecolor = "black")
    # end
    # PyPlot.savefig("metamat"*string(i)*"-1/metamat"*string(i)*"-barplot.png",dpi=200)

    plt = Plots.plot(showaxis=false,size=(400,400),title="Metamaterial $i")
    for line in barPlots[i]
        plot!(plt, line[1], line[2], lw =3, lc = :black, label = false)
    end
    # for polygon in  barPlots[i]
    #     plot!(plt,polygon,label=false,color=:grey)
    # end
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
    plt = Plots.plot(lw = 3, xlabel = "Strain", ylabel = "Overall Stress (kg/cm2)", palette = :tab10, ylims = (0,Inf))
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
materials = [ simulation.model.material for simulation in simulations ]

io = open("results.csv","a")
writedlm(io,transpose(weights),",")
writedlm(io,transpose(maxStrains),",")
writedlm(io,transpose(maxStresses),",")
writedlm(io,transpose(energyAbsorptions),",")
writedlm(io,transpose((m->m.isECC).(materials)),",")
writedlm(io,transpose((m->m.compressive).(materials)),",")
writedlm(io,transpose((m->m.tensile).(materials)),",")
writedlm(io,transpose((m->m.tensilePeak).(materials)),",")
writedlm(io,transpose((m->m.peakStrain).(materials)),",")
writedlm(io,transpose((m->m.crackStrainRatio).(materials)),",")
close(io)

# General plots over all simulations

strainsPlt = Plots.plot(xlabel = "Area (cm2)",ylabel = "Failure strain", palette = :tab10, legend =:topleft)
for i in 1:repeatedSimulation
    plot!(strainsPlt, weights[:,i], maxStrains[:,i], seriestype = :scatter, label=labels[i])
end
png(strainsPlt,"strains.png")

sleep(0.1)
stressPlt = Plots.plot(xlabel = "Area (cm2)",ylabel = "Failure stress", palette = :tab10, legend =:topleft)
for i in 1:repeatedSimulation
     plot!(stressPlt,weights[:,i],maxStresses[:,i], seriestype = :scatter, label=labels[i])
end
png(stressPlt,"stresses.png")

sleep(0.1)
energyPlt = Plots.plot(xlabel = "Area (cm2)",ylabel = "Absorbed energy at failure", palette = :tab10, legend =:topleft)
for i in 1:repeatedSimulation
    plot!(energyPlt,weights[:,i],energyAbsorptions[:,i], seriestype = :scatter,label=labels[i])
end
png(energyPlt,"energies.png")

save("skeletons.jld","skeletons",skeletonLinks)
save("materials.jld","materials",materials)

currentTime = Dates.format(now(),"dd-mm-yyyy_HH:MM:SS")
run(`mkdir Batch_$currentTime`)
folders = glob("metamat*/")
plots = glob("*.png")
run(`mv $folders $plots results.csv skeletons.jld materials.jld Batch_$currentTime/`)
