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
