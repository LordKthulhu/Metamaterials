using DelimitedFiles
using Glob
using Plots

const possibleArgs = ["weight","strain","stress","energy","isECC","compressive","tensile","tensilePeak","failureStrain","crackStrainRatio"]
const labels = [ "Area (cm2)", "Max Strain", "Max Stress (kgf/cm2)", "Energy", "isECC", "Material compressive strength (MPa)",
                "Material tensile strength (MPa)", "Material strain at peak in tension", "Material strain at failure", "First cracking stress/tensile strength"]

function parseArguments()
    currentArg = 1
    while currentArg <= length(ARGS)
        if ARGS[currentArg] in ["-X","-x"]
            arg = ARGS[currentArg+1]
            if arg in possibleArgs
                X = findall(a->a==arg,possibleArgs)[1]-1
            else
                throw(ArgumentError("Unsupported argument : $(ARGS[currentArg]). For available arguments, type -help."))
            end
            currentArg += 2
        elseif ARGS[currentArg] in ["-Y","-y"]
            arg = ARGS[currentArg+1]
            if arg in possibleArgs
                Y = findall(a->a==arg,possibleArgs)[1]-1
            else
                throw(ArgumentError("Unsupported argument : $(ARGS[currentArg]). For available arguments, type -help."))
            end
            currentArg += 2
        elseif ARGS[currentArg] == "-help"
            println("Supported arguments are:\n
               -X           Value for X axis : $(prod([arg*", " for arg in possibleArgs])[1:end-2]).\n
               -Y           Value for Y axis : $(prod([arg*", " for arg in possibleArgs])[1:end-2]).\n")
            exit(0)
        else
            throw(ArgumentError("Unsupported argument : $(ARGS[currentArg]). For available arguments, type -help."))
        end
    end
    (@isdefined X) ? nothing : X = 0
    (@isdefined Y) ? nothing : Y = 2
    (X,Y)
end

X,Y = parseArguments()

data = readdlm("results.csv",',')

sets = size(data,1)÷length(possibleArgs)

plt = plot(xlabel = labels[X+1], ylabel = labels[Y+1], palette=:tab10, size=(1000,1000))
for i in 1:sets
    plot!(plt,data[i+X*sets,:],data[i+Y*sets,:], seriestype=:scatter, label=false, smooth = true)
end

png(plt,"$(ARGS[2])VS$(ARGS[4]).png")
