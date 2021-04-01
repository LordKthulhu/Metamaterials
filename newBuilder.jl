using Plots
using Printf
using DelimitedFiles
using Statistics
using Glob
using Dates
using ProgressMeter

function execute(cmd::Cmd)
  out = Pipe()
  err = Pipe()

  process = run(pipeline(ignorestatus(cmd), stdout=out, stderr=err))
  close(out.in)
  close(err.in)
  code = process.exitcode
  code == 1 ? display("Error occured") : nothing
  code
end

function forceSteps(filename)
    force = 0
    forceList = []
    for line in readlines(filename * "/" * filename * "-MECH.nod",keep=true)[3:end]
        if line[2:5] == "STEP"
            push!(forceList,force)
        elseif line[2:5] == "NODE"
            force = 0
        else
            force += abs( parse( Float64, split(line)[7] ))
        end
    end
    push!(forceList,force)
    forceList
end

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

function nodeLine(n,x,y,z,loadPoints)
    sn = string(n)
    sx = @sprintf("%.3f",x); sy = @sprintf("%.3f",y); sz = @sprintf("%.3f",z)
    restraint = collect("010")

    if z == 0
        push!(loadPoints, [n,-1])
        restraint[3] = '1'
    elseif z == 12*(H-1)*unitSize
        push!(loadPoints, [n,1])
        restraint[3] = '1'
    elseif x == 12*(H-1)*unitSize
        restraint[1] = '1'
    end
    "NODE " * " "^(5-length(sn)) * sn * " "^(10-length(sx)) * sx * " "^(10-length(sy)) * sy * " "^(10-length(sz)) * sz * "                    " * join(restraint) * "\n"
end

function elementLine(n,P)

    sn = string(n)
    P = string.(P)

    "PLAT " *" "^(5-length(sn)) * sn * " "^(5-length(P[1])) * P[1] * " "^(5-length(P[2])) * P[2] * " "^(5-length(P[3])) * P[3] *
    " "^(5-length(P[4])) * P[4] * "    0    0    0    0    0    0    0    0    0\n" *
    "    0    0    0    0    0    0    0    3       0.0       1.0    3000.0 0.0025000       0.0                                                                                                         0.005     0.046       0.8\n" *
    "    450.00    48.000      0.17    0\n" *
    " 2100000.0    4000.0   0.00000 2100000.0    4000.0   0.00000       0.0       0.0       0.0\n"
end

function elementLinePlain(n,P)

    sn = string(n)
    P = string.(P)

    "PLAT " *" "^(5-length(sn)) * sn * " "^(5-length(P[1])) * P[1] * " "^(5-length(P[2])) * P[2] * " "^(5-length(P[3])) * P[3] *
    " "^(5-length(P[4])) * P[4] * "    0    0    0    0    0    0    0    0    0\n" *
    "    0    0    0    0    0    0    0    0       0.0       1.0    3000.0 0.0025000       0.0\n" *
    "    450.00    48.000      0.17    0\n" *
    " 2100000.0    4000.0   0.00000 2100000.0    4000.0   0.00000       0.0       0.0       0.0\n"
end

function loadLine(n,F)
    sn = string(n)
    sf = @sprintf("%.5f",F)
    "LOAD " * " "^(5-length(sn)) * sn * "   0.00000   0.00000" * " "^(10-length(sf)) * sf * "\n"
end


#### MAIN ####

iterations=parse(Int,ARGS[1])

@printf("Starting %d simulations on %d threads.\n",iterations,Threads.nthreads())

global maxStrains = zeros(iterations)
global maxStrainsPlain = zeros(iterations)
global maxSigmas = zeros(iterations)
global maxSigmasPlain = zeros(iterations)
global energyAbsorptions = zeros(iterations)
global energyAbsorptionsPlain = zeros(iterations)
global weights = zeros(iterations)

global barPlots = [ plot(showaxis=false,size=(400,400)) for i=1:iterations ]
progress = Progress(iterations, dt = 1, desc = "Computing progress... " , barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',), barlen=20)

const H = 3
const p = 0.7
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

basePoints = readdlm("points.csv",',')
baseElements = readdlm("elements.csv",',',Int)

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

        potLinks[:,1,1] = falses(H)
        potLinks[:,H,3:4] = falses(H,2)
        potLinks[H,:,:] = falses(H,4)

        loadPoints = []


        ##### Geometry generation #####


        for i in 1:H-1
            i == 1 ? nodes[1,:] = rand(H) .< p : nothing
            nodes[1,1] = true
            for j in 1:H
                if nodes[i,j]
                    potLinksCount = length(findall(potLinks[i,j,:]))
                    for l in 1:rand(1:potLinksCount)
                        index = rand(findall(potLinks[i,j,:]))
                        links[i,j,index] = true; potLinks[i,j,index] = false
                        if index == 1
                            nodes[i+1,j-1] = true
                        elseif index == 2
                            nodes[i+1,j] = true
                        elseif index == 3
                            nodes[i+1,j+1] = true
                            potLinks[i,j+1,1] = false
                        else
                            nodes[i,j+1] = true
                        end
                    end
                end
            end
        end

        for j in 1:H-1
            if nodes[H,j] && rand(Bool)
                links[H,j,4] = true
                nodes[H,j+1] = true
            end
        end


        ##### Geometry plot #####


        barPlots[iter] = plot(showaxis = false, size = (400,400))

        for i in 1:H, j in 1:H
            for index in findall(links[i,j,:])
                if index == 1
                    plot!(barPlots[iter], [j, j-1], [i,i+1], lc=:black, lw = 3, label = false)
                    plot!(barPlots[iter], [j, j-1], [2*(H)-i,2*(H)-(i+1)], lc=:black, lw = 3, label = false)
                    plot!(barPlots[iter], [2*(H)-j, 2*(H)-(j-1)], [i,(i+1)], lc=:black, lw = 3, label = false)
                    plot!(barPlots[iter], [2*(H)-j, 2*(H)-(j-1)], [2*(H)-i,2*(H)-(i+1)], lc=:black, lw = 3, label = false)
                elseif index == 2
                    plot!(barPlots[iter], [j, j], [i,i+1], lc=:black, lw = 3, label = false)
                    plot!(barPlots[iter], [j, j], [2*(H)-i,2*(H)-(i+1)], lc=:black, lw = 3, label = false)
                    plot!(barPlots[iter], [2*(H)-j, 2*(H)-j], [i,i+1], lc=:black, lw = 3, label = false)
                    plot!(barPlots[iter], [2*(H)-j, 2*(H)-j], [2*(H)-i,2*(H)-(i+1)], lc=:black, lw = 3, label = false)
                elseif index == 3
                    plot!(barPlots[iter], [j, j+1], [i,i+1], lc=:black, lw = 3, label = false)
                    plot!(barPlots[iter], [j, j+1], [2*(H)-i,2*(H)-(i+1)], lc=:black, lw = 3, label = false)
                    plot!(barPlots[iter], [2*(H)-j, 2*(H)-(j+1)], [i,i+1], lc=:black, lw = 3, label = false)
                    plot!(barPlots[iter], [2*(H)-j, 2*(H)-(j+1)], [2*(H)-i,2*(H)-(i+1)], lc=:black, lw = 3, label = false)
                else
                    plot!(barPlots[iter], [j, j+1], [i,i], lc=:black, lw = 3, label = false)
                    plot!(barPlots[iter], [j, j+1], [2*(H)-i,2*(H)-i], lc=:black, lw = 3, label = false)
                    plot!(barPlots[iter], [2*(H)-j, 2*(H)-(j+1)], [i,i], lc=:black, lw = 3, label = false)
                    plot!(barPlots[iter], [2*(H)-j, 2*(H)-(j+1)], [2*(H)-i,2*(H)-i], lc=:black, lw = 3, label = false)
                end
            end
        end

        ##### Meshing geometry #####

        run(`mkdir $filename $filename-plain`)
        datFile = open("$filename.dat","w")
        write(datFile, "Metamaterial Project\n1202000001000000000000002     0.700     0.000     0.000         0\nNODE\n")

        for index in findall(nodes)
            i = index[1]; j = index[2]
            nodeNumber = j + H*(i-1)
            startingPoint = 145*(nodeNumber-1)
            for pointNumber in 1:37
                write(datFile,nodeLine( startingPoint + pointNumber, unitSize * (12*(j-1) + basePoints[pointNumber,2]), 0, unitSize * (12*(i-1) + basePoints[pointNumber,3]), loadPoints ))
            end
            for link in findall(links[i, j,:])
                if link == 1
                    for pointNumber in 38:70
                        write(datFile,nodeLine( startingPoint + pointNumber, unitSize * (12*(j-1) + basePoints[pointNumber,2]), 0, unitSize * (12*(i-1) + basePoints[pointNumber,3]), loadPoints ))
                    end
                elseif link == 2
                    for pointNumber in 71:91
                        write(datFile,nodeLine( startingPoint + pointNumber, unitSize * (12*(j-1) + basePoints[pointNumber,2]), 0, unitSize * (12*(i-1) + basePoints[pointNumber,3]), loadPoints ))
                    end
                elseif link == 3
                    for pointNumber in 92:124
                        write(datFile,nodeLine( startingPoint + pointNumber, unitSize * (12*(j-1) + basePoints[pointNumber,2]), 0, unitSize * (12*(i-1) + basePoints[pointNumber,3]), loadPoints ))
                    end
                else
                    for pointNumber in 125:145
                        write(datFile,nodeLine( startingPoint + pointNumber, unitSize * (12*(j-1) + basePoints[pointNumber,2]), 0, unitSize * (12*(i-1) + basePoints[pointNumber,3]), loadPoints ))
                    end
                end
            end
        end
        write(datFile, "ELEM\n")

        close(datFile)
        run(`cp $filename.dat $filename-plain.dat`)
        datFile = open("$filename.dat","a")
        datFilePlain = open("$filename-plain.dat","a")

        for index in findall(nodes)
            nodeNumber = index[2] + H*(index[1]-1)
            startingPoint = 145*(nodeNumber-1)
            startingElement = 112*(nodeNumber-1)
            for elementNumber in 1:32
                write(datFile,elementLine( startingElement + elementNumber, startingPoint .+ baseElements[elementNumber,2:end] ))
                write(datFilePlain,elementLinePlain( startingElement + elementNumber, startingPoint .+ baseElements[elementNumber,2:end] ))
            end
            for link in findall(links[index[1], index[2],:])
                if link == 1
                    for elementNumber in 33:56
                        write(datFile,elementLine( startingElement + elementNumber, startingPoint .+ baseElements[elementNumber,2:end] ))
                        write(datFilePlain,elementLinePlain( startingElement + elementNumber, startingPoint .+ baseElements[elementNumber,2:end] ))

                    end
                elseif link == 2
                    for elementNumber in 57:72
                        write(datFile,elementLine( startingElement + elementNumber, startingPoint .+ baseElements[elementNumber,2:end] ))
                        write(datFilePlain,elementLinePlain( startingElement + elementNumber, startingPoint .+ baseElements[elementNumber,2:end] ))
                    end
                elseif link == 3
                    for elementNumber in 73:96
                        write(datFile,elementLine( startingElement + elementNumber, startingPoint .+ baseElements[elementNumber,2:end] ))
                        write(datFilePlain,elementLinePlain( startingElement + elementNumber, startingPoint .+ baseElements[elementNumber,2:end] ))
                    end
                else
                    for elementNumber in 97:112
                        write(datFile,elementLine( startingElement + elementNumber, startingPoint .+ baseElements[elementNumber,2:end] ))
                        write(datFilePlain,elementLinePlain( startingElement + elementNumber, startingPoint .+ baseElements[elementNumber,2:end] ))
                    end
                end
            end
        end

        write(datFile, "LOAD\n\n")
        write(datFilePlain, "LOAD\n\n")
        close(datFile)
        close(datFilePlain)

        #### File with restart option

        auxFile = open("$filename-restart.aux","w")
        lines = readlines("$filename.dat",keep=true)
        for i=1:length(lines)
            if i == 2
                write(auxFile,"1212000001000000000000002     0.700     0.000     0.000         0\n")
            else
                write(auxFile,lines[i])
            end
        end
        close(auxFile)

        auxFilePlain = open("$filename-plain-restart.aux","w")
        lines = readlines("$filename-plain.dat",keep=true)
        for i=1:length(lines)
            if i == 2
                write(auxFilePlain,"1212000001000000000000002     0.700     0.000     0.000         0\n")
            else
                write(auxFilePlain,lines[i])
            end
        end
        close(auxFilePlain)

        ##### Simulation #####

        strain = []
        sigmaIntern = []
        sigmaOverall = []

        datFile = open("$filename.dat","a")
        datFilePlain = open("$filename-plain.dat","a")

        time = 0
        nSteps = 0

        loadRange = [ i*dEpsilonMax for i=1:10 ]

        for i=1:10
            time += 0.001
            timeStr = @sprintf("%.4f",time)
            write(datFile, "STEP " * " "^(5-length(string(i)))*string(i) * " "^(10-length(timeStr))*timeStr * "     0.000     0.000     0.000                                             0.000\n")
            write(datFilePlain, "STEP " * " "^(5-length(string(i)))*string(i) * " "^(10-length(timeStr))*timeStr * "     0.000     0.000     0.000                                             0.000\n")
            for load in loadPoints
                write(datFile, loadLine(load[1], 12*(H-1)*unitSize*dEpsilonMax * load[2]))
                write(datFilePlain, loadLine(load[1], 12*(H-1)*unitSize*dEpsilonMax * load[2]))
            end
        end
        write(datFile, "END\n")
        close(datFile)
        write(datFilePlain, "END\n")
        close(datFilePlain)

        nSteps += 10

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

    plt2 = plot(strain,sigmaOverall, lw = 3, ylabel = "Overall Stress (kg/cm2)",lc="green", label="PVA-ECC")
    plot!(plt2,strain,sigmaOverallPlain, lw = 3, lc="red", label="Plain Concrete")
    plotfile = filename*"/"*filename*"-plot.png"
    png(plt2,plotfile)

    run(`rm $filename-restart.aux $filename-plain-restart.aux`)
    next!(progress)
end

for i=1:iterations
png(barPlots[i],"metamat"*string(i)*"/metamat"*string(i)*"-barplot.png")
end

io = open("results.csv","a")

writedlm(io,transpose(weights),",")
writedlm(io,transpose(maxStrains),",")
writedlm(io,transpose(maxSigmas),",")
writedlm(io,transpose(energyAbsorptions),",")
writedlm(io,transpose(maxStrainsPlain),",")
writedlm(io,transpose(maxSigmasPlain),",")
writedlm(io,transpose(energyAbsorptionsPlain),",")

close(io)

strainsPlt = plot(weights,maxStrains, seriestype = :scatter, xlabel = "Area (cm2)",ylabel = "Failure strain",label="PVA-ECC", color = :blue)
plot!(strainsPlt, weights,maxStrainsPlain, seriestype = :scatter, label="Plain Concrete", color = :red)
png(strainsPlt,"strains.png")

stressPlt = plot(weights,maxSigmas, seriestype = :scatter, xlabel = "Area (cm2)",ylabel = "Failure stress",label="PVA-ECC", color = :blue)
plot!(stressPlt, weights,maxSigmasPlain, seriestype = :scatter, label="Plain Concrete", color = :red)
png(stressPlt,"stresses.png")

energyPlt = plot(weights,energyAbsorptions, seriestype = :scatter, xlabel = "Area (cm2)",ylabel = "Absorbed energy at failure",label="PVA-ECC", color = :blue)
plot!(energyPlt, weights, energyAbsorptionsPlain, seriestype = :scatter, label="Plain Concrete", color = :red)
png(energyPlt,"energies.png")
