using PyCall
using PyPlot
using DelimitedFiles
using Glob

const possibleArgs = ["weight","strain","stress","energy","isECC","compressive","tensile","tensilePeak","failureStrain","crackStrainRatio"]
labels = [ "Area (cm2)", "Max Strain", "Max Stress (kgf/cm2)", "Energy", "isECC", "Material compressive strength (MPa)",
                "Material tensile strength (MPa)", "Material strain at peak in tension", "Material strain at failure", "First cracking stress/tensile strength"]


function hover(event)
    visLeft = [ a.get_visible() for a in annotLeft ]
    visRight = [ a.get_visible() for a in annotRight ]
    if event.inaxes == ax[1]
        cont,ind = sc[1].contains(event)
        if cont
            annotLeft[ind["ind"][1]+1].set_visible(true)
            annotRight[ind["ind"][1]+1].set_visible(true)
            fig.canvas.draw_idle()
        else
            for index in findall(visRight)
                annotRight[index].set_visible(false)
                fig.canvas.draw_idle()
            end
            for index in findall(visLeft)
                annotLeft[index].set_visible(false)
                fig.canvas.draw_idle()
            end
        end
    end
end

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
        elseif ARGS[currentArg] == "-normalize"
            normalize = true
            currentArg += 1
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
    (@isdefined normalize) ? nothing : normalize = false
    (X,Y,normalize)
end

X,Y,normalize = parseArguments()

data = readdlm("results.csv",',')

if normalize
    labels = "Normalized " .* labels
    for i in reverse(1:size(data,1))
        data[i,:] = data[i,:]./data[1,:]
    end
end


const N = length(glob("*-1/"))

barModels = [ PyPlot.imread("metamat$i-1/metamat$i-barplot.png",format="png") for i=1:N ]
barImageBoxes = [ matplotlib.offsetbox.OffsetImage(barModels[i], zoom=0.5) for i=1:N ]

stressPlots = [ PyPlot.imread("metamat$i-1/metamat$i-plot.png",format="png") for i=1:N ]
stressImageBoxes = [ matplotlib.offsetbox.OffsetImage(stressPlots[i], zoom=0.9) for i=1:N ]

fig,ax = PyPlot.subplots(1,2)

sets = size(data,1)÷length(possibleArgs)

trendlines = []

for i in 1:sets
    x=data[i+X*sets,:]; y=data[i+Y*sets,:]
    push!(trendlines,(x'*x)\(x'*y))
end

sc = [ ax[1].scatter(data[i+X*sets,:],data[i+Y*sets,:]) for i in 1:sets ]
for i in 1:sets
    ax[1].plot(data[i+X*sets,:],data[i+X*sets,:].*trendlines[i])
end
ax[1].set_xlabel(labels[X+1])
ax[1].set_ylabel(labels[Y+1])

ax[2].set_axis_off()

for i=1:N
    barImageBoxes[i].image.axes = ax[2]
    stressImageBoxes[i].image.axes = ax[2]
end

annotLeft = []
annotRight = []
for i=1:N
    ab1 = matplotlib.offsetbox.AnnotationBbox(barImageBoxes[i], [0.5,0.15])
    push!(annotLeft,ax[2].add_artist(ab1))
    ab2 = matplotlib.offsetbox.AnnotationBbox(stressImageBoxes[i], [0.5,0.7])
    push!(annotRight,ax[2].add_artist(ab2))
end
for a in annotLeft
    a.set_visible(false)
end
for a in annotRight
    a.set_visible(false)
end

fig.canvas.mpl_connect("motion_notify_event", hover)

PyPlot.show()
