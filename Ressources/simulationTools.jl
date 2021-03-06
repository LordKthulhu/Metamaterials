# Tools for parsing arguments and executing bash scripts

function pBar(len,text; dt=1)
    Progress(len, dt = dt, desc = text , barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',), barlen=20)
end

function execute(cmd::String)
    if Sys.iswindows()
        executeWindows(cmd)
    else
        executeLinux(cmd)
    end
end


function executeLinux(cmd::String)
    err = Pipe()
    process = run(pipeline(ignorestatus(`./runCOM3.sh $cmd`),stderr=err))
    close(err.in)
    close(err.out)
    code = process.exitcode
    code == 1 ? display("Error occured") : nothing
    code
end

function executeWindows(cmd::String)
    Shell.run("runCOM3.sh $cmd")
    0
end

function parseArguments()
    currentArg = 1
    while currentArg <= length(ARGS)
        if ARGS[currentArg] == "-batchsize"
            iterations = parse(Int,ARGS[currentArg+1])
            currentArg += 2
        elseif ARGS[currentArg] == "-material-random"
            randomMat = true
            currentArg += 1
        elseif ARGS[currentArg] == "-patternsize"
            H = parse(Int,ARGS[currentArg+1])
            currentArg += 2
        elseif ARGS[currentArg] == "-parametric"
            parameters = ARGS[currentArg+1]
            currentArg += 2
        elseif ARGS[currentArg] == "-angle"
            angles = readdlm(ARGS[currentArg+1],',')
            currentArg += 2
        elseif ARGS[currentArg] == "-skeletons"
            skeletons = ARGS[currentArg+1]
            currentArg += 2
        elseif ARGS[currentArg] == "-output"
            if ARGS[currentArg+1] in ["none","all"]
                output = ARGS[currentArg+1]
                currentArg += 2
            else
                throw(ArgumentError("Unsupported argument : $(ARGS[currentArg+1]). For available arguments, type -help."))
            end
        elseif ARGS[currentArg] == "-help"
            println("Supported arguments are:\n
               -batchsize            Specify the number of desired simulations\n
               -material-random      Base-material will be randomly picked\n
               -patternsize          Size of the pattern grid (default is 3)\n
               -parametric           Please specify a parameters file for parametric study\n
               -angle                Please specify a parameter file for angles\n
               -output               output files written by COM3 : none or all (default is none)\n")
            exit(0)
        else
            throw(ArgumentError("Unsupported argument : $(ARGS[currentArg]). For available arguments, type -help."))
        end
    end
    (@isdefined H) ? nothing : H = 3
    (@isdefined iterations) ? nothing : iterations = 1
    (@isdefined randomMat) ? nothing : randomMat = false
    (@isdefined parameters) ? nothing : parameters = "none"
    (@isdefined angles) ? nothing : angles = [0]
    (@isdefined output) ? nothing : output = "none"
    (@isdefined skeletons) ? nothing : skeletons = "none"
    (H,iterations,randomMat,parameters,angles,output,skeletons)
end

function runSteps(simulation::Simulation,angle,direction)

    if simulation.step == 0
        restart = false; nSteps = 10
    else
        restart = true; nSteps = 5
    end

    addSteps(simulation.model.loadPoints, simulation.step, simulation.dEpsilon,
    simulation.model.unitSize, simulation.model.size, simulation.filename, restart = restart, nSteps = nSteps,angle = angle, direction=direction)
    filename = simulation.filename
    run(`mv $filename.dat $filename`)
    simulation.exit = execute(filename)
    forces = forceSteps(filename)
    loadRange = [ (length(simulation.strain)+i)*simulation.dEpsilon*(1+simulation.step÷20) for i=1:nSteps ]  #loadRange = [ (length(simulation.strain)+i)*simulation.dEpsilon for i=1:nSteps ]
    if restart
        a = "$filename/$filename-MECH.nod"
        b = "$filename/aux.nod"
        open(io -> run(pipeline(`cat $a`, stdout=io)), b, "a")
        a = "$filename/$filename-MECH.crk"
        b = "$filename/aux.crk"
        open(io -> run(pipeline(`cat $a`, stdout=io)), b, "a")
        a = "$filename/$filename-MECH.fld"
        b = "$filename/aux.fld"
        open(io -> run(pipeline(`cat $a`, stdout=io)), b, "a")
    else
        run(`mv $filename/$filename-MECH.nod $filename/aux.nod`)
        run(`mv $filename/$filename-MECH.crk $filename/aux.crk`)
        run(`mv $filename/$filename-MECH.fld $filename/aux.fld`)
    end

    append!(simulation.strain, loadRange[1:length(forces)])
    area = 12*(simulation.model.size-1)*simulation.model.unitSize*1.0
    append!(simulation.stress, forces./area)
    simulation.step += nSteps
end

function energy(strain,stress)
    energy = 0
    for i=2:length(strain)
        energy += 1/2*(stress[i]+stress[i-1])*(strain[i]-strain[i-1])
    end
    energy
end
