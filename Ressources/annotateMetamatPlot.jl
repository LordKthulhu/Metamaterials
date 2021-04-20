using PyCall
using PyPlot
using DelimitedFiles
using Glob

function hover(event)
    visLeft = [ a.get_visible() for a in annotLeft ]
    visRight = [ a.get_visible() for a in annotRight ]
    if event.inaxes == ax[1]
        cont,ind = sc.contains(event)
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

data = readdlm("results.csv",',')

const N = length(glob("*-1/"))

barModels = [ PyPlot.imread("metamat$i-1/metamat$i-barplot.png",format="png") for i=1:N ]
barImageBoxes = [ matplotlib.offsetbox.OffsetImage(barModels[i], zoom=0.5) for i=1:N ]

stressPlots = [ PyPlot.imread("metamat$i-1/metamat$i-plot.png",format="png") for i=1:N ]
stressImageBoxes = [ matplotlib.offsetbox.OffsetImage(stressPlots[i], zoom=0.9) for i=1:N ]

fig,ax = PyPlot.subplots(1,2)

sc = ax[1].scatter(data[1,:],data[1+size(data,1)÷4,:])
ax[1].set_xlabel("Area (cm2)")
ax[1].set_ylabel("Max Stress (kgf/cm2)")

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
