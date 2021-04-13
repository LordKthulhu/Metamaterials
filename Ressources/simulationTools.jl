# Tools for parsing arguments and executing bash scripts

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
        elseif ARGS[currentArg] == "-help"
            println("Supported arguments are:\n
               -batchsize            Specify the number of desired simulations\n
               -material-random      base-material will be randomly picked\n
               -patternsize          size of the pattern grid (default is 3)\n")
            exit(0)
        else
            throw(ArgumentError("Unsupported argument : $(ARGS[currentArg]). For available arguments, type -help."))
        end
    end
    (@isdefined H) ? nothing : H = 3
    (@isdefined iterations) ? nothing : iterations = 1
    (@isdefined randomMat) ? nothing : randomMat = false
    (H,iterations,randomMat)
end

function runSteps(simulation)

    if simulation.step == 0
        restart = false; nSteps = 10
    else
        restart = true; nSteps = 5
    end

    addSteps(simulation.model.loadPoints, simulation.step, simulation.dEpsilon,
    simulation.model.unitSize, simulation.filename, restart = restart, nSteps = nSteps)
    filename = simulation.filename
    run(`mv $filename.dat $filename`)
    simulation.exit = execute(`./runCOM3.sh $filename`)

    forces = forceSteps(filename)
    loadRange = [ (length(simulation.strain)+i)*simulation.dEpsilon for i=1:nSteps ]

    append!(simulation.strain, loadRange[1:length(forces)])
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
