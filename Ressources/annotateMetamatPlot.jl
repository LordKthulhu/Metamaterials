using PyCall
using PyPlot
using DelimitedFiles
using Glob

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

ylabels = [ "Max Strain", "Max Stress (kgf/cm2)", "Energy"]

mode = 2

if length(ARGS)>0
    if ARGS[1] == "strain"
        mode = 1
    elseif ARGS[1] == "stress"
        mode = 2
    elseif ARGS[1] == "energy"
        mode = 3
    else
        throw(ArgumentError("Unsupported argument : $(ARGS[1]). Possible arguments are strain, stress, or energy.\n"))
    end
end

data = readdlm("results.csv",',')

const N = length(glob("*-1/"))

barModels = [ PyPlot.imread("metamat$i-1/metamat$i-barplot.png",format="png") for i=1:N ]
barImageBoxes = [ matplotlib.offsetbox.OffsetImage(barModels[i], zoom=0.5) for i=1:N ]

stressPlots = [ PyPlot.imread("metamat$i-1/metamat$i-plot.png",format="png") for i=1:N ]
stressImageBoxes = [ matplotlib.offsetbox.OffsetImage(stressPlots[i], zoom=0.9) for i=1:N ]

fig,ax = PyPlot.subplots(1,2)

sets = size(data,1)÷4

trendlines = []

for i in 1:sets
    x=data[i,:]; y=data[i+mode*sets,:]
    push!(trendlines,(x'*x)\(x'*y))
end

sc = [ ax[1].scatter(data[i,:],data[i+mode*sets,:]) for i in 1:sets ]
for i in 1:sets
    ax[1].plot(data[i,:],data[i,:].*trendlines[i])
end
ax[1].set_xlabel("Area (cm2)")
ax[1].set_ylabel(ylabels[mode])

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
