using Plots
using Printf
using DelimitedFiles
using Statistics
using Glob
using Dates
using ProgressMeter

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
global maxStrainsPlain = zeros(iterations)
global MaxStresses = zeros(iterations)
global MaxStressesPlain = zeros(iterations)
global energyAbsorptions = zeros(iterations)
global energyAbsorptionsPlain = zeros(iterations)
global weights = zeros(iterations)

global plt2 = [plot() for i in 1:iterations]

global barPlots = [[] for i in 1:iterations]

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

progress = pBar(iterations,"Computing progess... ")

Threads.@threads for iter = 1:iterations
    filename = "metamat" * string(iter)

    strain = []
    strainPlain = []
    stress = []
    stressPlain = []
    exit = 1
    exitPlain = 1

    while exit == 1 || exitPlain == 1

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
        strainPlain = []
        stress = []
        stressPlain = []
        time = 0
        startStep = 0
        startStepPlain = 0

        maxStress, startStep, exit = runSteps(strain, stress, startStep, filename, dEpsilonMax, loadPoints)
        maxStressPlain, startStepPlain, exitPlain = runSteps(strainPlain, stressPlain, startStepPlain, filename * "-plain", dEpsilonMax, loadPoints)


        while stress[end]/maxStress > 0.5 && exit == 0
            maxStress, startStep, exit = runSteps(strain, stress, startStep, filename, dEpsilonMax, loadPoints)
            if exitPlain == 0
                maxStressPlain, startStepPlain, exitPlain = runSteps(strainPlain, stressPlain, startStepPlain, filename * "-plain", dEpsilonMax, loadPoints)
            end
        end

        mechFiles = glob("$(filename)/MECHIFI*")
        run(`rm $mechFiles`)
        mechFiles = glob("$(filename)-plain/MECHIFI*")
        run(`rm $mechFiles`)
        weights[iter] = sum(nodes .* nodeWeights) + sum(links .* linkWeights)

        if exit == 1 || exitPlain == 1
            println("Simulation $iter failed. Starting over.")
            run(`rm -r $filename $filename-plain`)
        end
    end

    ###### If all COM3 runs terminate, go on to results processing ######

    run(`rm $filename/$filename-MECH.crk $filename/$filename-MECH.fld
            $filename/$filename-MECH.int $filename/$filename-MECH.tmp`)

    #strain = strain[1:end-1]
    maxStress,index = findmax(stress)
    maxStressPlain,indexPlain = findmax(stressPlain)

    MaxStresses[iter] = maxStress
    maxStrains[iter] = strain[index]
    MaxStressesPlain[iter] = maxStressPlain
    maxStrainsPlain[iter] = strainPlain[indexPlain]

    energyAbsorptions[iter] = energy(strain[1:index],stress[1:index])
    energyAbsorptionsPlain[iter] = energy(strainPlain[1:indexPlain],stressPlain[1:indexPlain])


    io = open(filename*"/"*filename*"-results.csv","a")

    writedlm(io,transpose(strain),",")
    writedlm(io,transpose(stress),",")
    close(io)

    plt2[iter] = plot(strain,stress, lw = 3, ylabel = "Overall Stress (kg/cm2)",lc="green", label="PVA-ECC")
    plot!(plt2[iter],strainPlain,stressPlain, lw = 3, lc="red", label="Plain Concrete")

    run(`rm $filename-restart.aux $filename-plain-restart.aux`)
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
    png(plt2[i],plotfile)
    next!(progress)
    sleep(0.5)
end

# CSV output

io = open("results.csv","a")
writedlm(io,transpose(weights),",")
writedlm(io,transpose(maxStrains),",")
writedlm(io,transpose(MaxStresses),",")
writedlm(io,transpose(energyAbsorptions),",")
writedlm(io,transpose(maxStrainsPlain),",")
writedlm(io,transpose(MaxStressesPlain),",")
writedlm(io,transpose(energyAbsorptionsPlain),",")
close(io)

# General plots over all simulations

strainsPlt = plot(weights,maxStrains, seriestype = :scatter, xlabel = "Area (cm2)",ylabel = "Failure strain",label="PVA-ECC", color = :blue)
plot!(strainsPlt, weights,maxStrainsPlain, seriestype = :scatter, label="Plain Concrete", color = :red)
png(strainsPlt,"strains.png")

sleep(0.5)
stressPlt = plot(weights,MaxStresses, seriestype = :scatter, xlabel = "Area (cm2)",ylabel = "Failure stress",label="PVA-ECC", color = :blue)
plot!(stressPlt, weights,MaxStressesPlain, seriestype = :scatter, label="Plain Concrete", color = :red)
png(stressPlt,"stresses.png")

sleep(0.5)
energyPlt = plot(weights,energyAbsorptions, seriestype = :scatter, xlabel = "Area (cm2)",ylabel = "Absorbed energy at failure",label="PVA-ECC", color = :blue)
plot!(energyPlt, weights, energyAbsorptionsPlain, seriestype = :scatter, label="Plain Concrete", color = :red)
png(energyPlt,"energies.png")
