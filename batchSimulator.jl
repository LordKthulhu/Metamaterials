using Plots
using Printf
using DelimitedFiles
using Statistics
using Glob
using Dates
using ProgressMeter

include("Ressources/consoleTools.jl")
include("Ressources/COM3RWTools.jl")
include("Ressources/geometryTools.jl")

ENV["GKSwstype"] = "100"

function runSteps(strain, sigmaOverall, filename, dEpsilonMax,loadPoints,maxOverall,exit,sigmaOverallPlain)

    nSteps = length(strain)

    run(`cp $(filename)-restart.aux $(filename).dat`)
    datFile = open(filename * ".dat","a")

    loadRange = [ strain[end]+i*dEpsilonMax for i=1:5 ]

    for i=1:5
        time = 0.001*(nSteps+i)
        timeStr = @sprintf("%.4f",time)
        write(datFile, "STEP " * " "^(5-length(string(i+nSteps)))*string(i+nSteps) * " "^(10-length(timeStr))*timeStr * "     0.000     0.000     0.000                                             0.000\n")
        for load in loadPoints
            write(datFile, loadLine(load[1], 12*(H-1)*unitSize*dEpsilonMax * load[2]))
        end
    end
    write(datFile, "END\n")
    close(datFile)

    run(`cp $(filename)-plain-restart.aux $(filename)-plain.dat`)
    datFile = open(filename * "-plain.dat","a")

    loadRange = [ strain[end]+i*dEpsilonMax for i=1:5 ]

    for i=1:5
        time = 0.001*(nSteps+i)
        timeStr = @sprintf("%.4f",time)
        write(datFile, "STEP " * " "^(5-length(string(i+nSteps)))*string(i+nSteps) * " "^(10-length(timeStr))*timeStr * "     0.000     0.000     0.000                                             0.000\n")
        for load in loadPoints
            write(datFile, loadLine(load[1], 12*(H-1)*unitSize*dEpsilonMax * load[2]))
        end
    end
    write(datFile, "END\n")
    close(datFile)

    run(`mv $filename.dat $filename`)
    exit = execute(`./runCOM3.sh $filename`)

    run(`mv $filename-plain.dat $filename-plain`)
    exit = execute(`./runCOM3.sh $filename-plain`)

    forces = forceSteps(filename)

    forcesPlain = forceSteps(filename*"-plain")

    append!(strain, loadRange[1:length(forces)])
    append!(sigmaOverall, forces./area)
    append!(sigmaOverallPlain, forcesPlain./area)

    max(maxOverall,maximum(sigmaOverall[end-5:end]))
end

function energy(strain,sigmaOverall)
    energy = 0
    for i=2:length(strain)
        energy += 1/2*(sigmaOverall[i]+sigmaOverall[i-1])*(strain[i]-strain[i-1])
    end
    energy
end

#### MAIN ####

const H,iterations,randomMat = parseArguments()

@printf("Starting %d simulations on %d threads.\n",iterations,Threads.nthreads())

global maxStrains = zeros(iterations)
global maxStrainsPlain = zeros(iterations)
global maxSigmas = zeros(iterations)
global maxSigmasPlain = zeros(iterations)
global energyAbsorptions = zeros(iterations)
global energyAbsorptionsPlain = zeros(iterations)
global weights = zeros(iterations)

global plt2 = [plot() for i in 1:iterations]

global barPlots = [[] for i in 1:iterations]
progress = Progress(iterations, dt = 1, desc = "Computing progress... " , barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',), barlen=20)

const unitSize = 1
const dEpsilonMax = 2e-4
const area = 12*(H-1)*unitSize*1.0

nodeWeights = [4*ones(1,H);8*ones(H-2,H);4*ones(1,H)]

nodeWeights[:,1] = 4*ones(H); nodeWeights[:,end] = 4*ones(H)
nodeWeights[1,1] = 2; nodeWeights[1,end] = 2; nodeWeights[end,1] = 2; nodeWeights[end,end] = 2
nodeWeights = unitSize^2 .* nodeWeights

linkWeights = zeros(H,H,4)
linkWeights[:,:,1] = (20*sqrt(2)-2)*unitSize^2 .* ones(H,H); linkWeights[:,:,3] = (20*sqrt(2)-2)*unitSize^2 .* ones(H,H)
linkWeights[:,:,2] = (20-4*sqrt(2))*unitSize^2 .* ones(H,H); linkWeights[:,:,4] = (20-4*sqrt(2))*unitSize^2 .* ones(H,H)

basePoints = readdlm("Ressources/points.csv",',')
baseElements = readdlm("Ressources/elements.csv",',',Int)

for index in findall(x -> x > 1000, baseElements)
    if baseElements[index] < 2000
        baseElements[index] = baseElements[index]%1000 + (H-1)*145
    elseif baseElements[index] < 3000
        baseElements[index] = baseElements[index]%2000 + H*145
    elseif baseElements[index] < 4000
        baseElements[index] = baseElements[index]%3000 + (H+1)*145
    else
        baseElements[index] = baseElements[index]%4000 + 145
    end
end

################################################################################
#####                             MAIN LOOP                                #####
################################################################################

Threads.@threads for iter = 1:iterations
    filename = "metamat" * string(iter)

    strain = []
    sigmaIntern = []
    sigmaOverall = []
    sigmaOverallPlain = []
    exit = 1

    while exit == 1

        nodes = falses(H, H)
        links = falses(H, H, 4)
        potLinks = trues(H, H, 4)
        loadPoints = []

        potLinks[:,1,1] = falses(H)
        potLinks[:,H,3:4] = falses(H,2)
        potLinks[H,:,:] = falses(H,4)

        randomGeometry(H,nodes,links,potLinks) # Geometry generation

        plotGeometry(barPlots,iter,links,H) # Adding data for geometry plot

        ########################################################################
        #####                   COM3 FILES GENERATION                      #####
        ########################################################################

        # Main files

        run(`mkdir $filename $filename-plain`)
        datFile = open("$filename.dat","w")
        write(datFile, "Metamaterial Project\n1202000001000200000000002     0.700     0.000     0.000         0\nNODE\n")
        writeNodes(nodes, links, datFile, basePoints, loadPoints, H, unitSize)
        write(datFile, "ELEM\n")
        close(datFile)
        run(`cp $filename.dat $filename-plain.dat`)
        datFile = open("$filename.dat","a")
        datFilePlain = open("$filename-plain.dat","a")
        writeElements(nodes, links, datFile, baseElements, H)
        writeElementsPlain(nodes, links, datFilePlain, baseElements, H)
        write(datFile, "LOAD\n\n")
        write(datFilePlain, "LOAD\n\n")
        close(datFile)
        close(datFilePlain)

        # File with restart option

        auxFile = open("$filename-restart.aux","w")
        lines = readlines("$filename.dat",keep=true)
        for i=1:length(lines)
            if i == 2
                write(auxFile,"1212000001000200000000002     0.700     0.000     0.000         0\n")
            else
                write(auxFile,lines[i])
            end
        end
        close(auxFile)

        auxFilePlain = open("$filename-plain-restart.aux","w")
        lines = readlines("$filename-plain.dat",keep=true)
        for i=1:length(lines)
            if i == 2
                write(auxFilePlain,"1212000001000200000000002     0.700     0.000     0.000         0\n")
            else
                write(auxFilePlain,lines[i])
            end
        end
        close(auxFilePlain)

        ########################################################################
        #####                        SIMULATIONS                           #####
        ########################################################################

        strain = []
        sigmaIntern = []
        sigmaOverall = []
        time = 0
        step = 0
        loadRange = [ i*dEpsilonMax for i=1:10 ]

        addSteps(loadPoints, step, dEpsilonMax, unitSize, filename)
        addSteps(loadPoints, step, dEpsilonMax, unitSize, filename*"-plain")

        step += 10

        run(`mv $filename.dat $filename`)
        exit = execute(`./runCOM3.sh $filename`)

        run(`mv $filename-plain.dat $filename-plain`)
        execute(`./runCOM3.sh $filename-plain`)

        append!(strain, loadRange)
        append!(sigmaOverall, forceSteps(filename)./area)
        append!(sigmaOverallPlain, forceSteps(filename*"-plain")./area)

        maxOverall = maximum(sigmaOverall)

        while sigmaOverall[end]/maxOverall > 0.5 && exit == 0
            maxOverall = runSteps(strain, sigmaOverall, filename, dEpsilonMax, loadPoints,maxOverall,exit,sigmaOverallPlain)
        end

        mechFiles = glob("$(filename)/MECHIFI*")
        run(`rm $mechFiles`)
        mechFiles = glob("$(filename)-plain/MECHIFI*")
        run(`rm $mechFiles`)
        weights[iter] = sum(nodes .* nodeWeights) + sum(links .* linkWeights)
    end

    ###### If all COM3 runs terminate, go on to results processing ######

    run(`rm $filename/$filename-MECH.crk $filename/$filename-MECH.fld
            $filename/$filename-MECH.int $filename/$filename-MECH.tmp`)

    strain = strain[1:end-1]
    maxSigma,index = findmax(sigmaOverall)
    maxSigmaPlain,indexPlain = findmax(sigmaOverallPlain)

    maxSigmas[iter] = maxSigma
    maxStrains[iter] = strain[index]
    maxSigmasPlain[iter] = maxSigmaPlain
    maxStrainsPlain[iter] = strain[indexPlain]

    energyAbsorptions[iter] = energy(strain[1:index],sigmaOverall[1:index])
    energyAbsorptionsPlain[iter] = energy(strain[1:indexPlain],sigmaOverallPlain[1:indexPlain])


    io = open(filename*"/"*filename*"-results.csv","a")

    writedlm(io,transpose(strain),",")
    writedlm(io,transpose(sigmaOverall),",")
    close(io)

    plt2[iter] = plot(strain,sigmaOverall, lw = 3, ylabel = "Overall Stress (kg/cm2)",lc="green", label="PVA-ECC")
    plot!(plt2[iter],strain,sigmaOverallPlain, lw = 3, lc="red", label="Plain Concrete")

    run(`rm $filename-restart.aux $filename-plain-restart.aux`)
    next!(progress)
end

################################################################################
#####                           RESULTS OUTPUT                             #####
################################################################################

# Geometry plots

for i=1:iterations
    plt = plot(showaxis=false,size=(400,400))
    for line in barPlots[i]
        plot!(plt, line[1], line[2], lw =3, lc = :black, label = false)
    end
    png(plt,"metamat"*string(i)*"/metamat"*string(i)*"-barplot.png")
    sleep(0.5)
end

# Stress-strain plots

for i=1:iterations
    plotfile = "metamat$i/metamat$i-plot.png"
    png(plt2[i],plotfile)
    sleep(0.5)
end

# CSV output

io = open("results.csv","w")
writedlm(io,transpose(weights),",")
writedlm(io,transpose(maxStrains),",")
writedlm(io,transpose(maxSigmas),",")
writedlm(io,transpose(energyAbsorptions),",")
writedlm(io,transpose(maxStrainsPlain),",")
writedlm(io,transpose(maxSigmasPlain),",")
writedlm(io,transpose(energyAbsorptionsPlain),",")
close(io)

# General plots over all simulations

strainsPlt = plot(weights,maxStrains, seriestype = :scatter, xlabel = "Area (cm2)",ylabel = "Failure strain",label="PVA-ECC", color = :blue)
plot!(strainsPlt, weights,maxStrainsPlain, seriestype = :scatter, label="Plain Concrete", color = :red)
png(strainsPlt,"strains.png")

sleep(0.5)
stressPlt = plot(weights,maxSigmas, seriestype = :scatter, xlabel = "Area (cm2)",ylabel = "Failure stress",label="PVA-ECC", color = :blue)
plot!(stressPlt, weights,maxSigmasPlain, seriestype = :scatter, label="Plain Concrete", color = :red)
png(stressPlt,"stresses.png")

sleep(0.5)
energyPlt = plot(weights,energyAbsorptions, seriestype = :scatter, xlabel = "Area (cm2)",ylabel = "Absorbed energy at failure",label="PVA-ECC", color = :blue)
plot!(energyPlt, weights, energyAbsorptionsPlain, seriestype = :scatter, label="Plain Concrete", color = :red)
png(energyPlt,"energies.png")
