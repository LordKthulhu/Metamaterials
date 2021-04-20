using PyCall
using PyPlot
using DelimitedFiles

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


barModels = [ PyPlot.imread("metamat$i/metamat$i-barplot.png",format="png") for i=1:299]
barImageBoxes = [ matplotlib.offsetbox.OffsetImage(barModels[i], zoom=0.5) for i=1:299 ]

stressPlots = [ PyPlot.imread("metamat$i/metamat$i-plot.png",format="png") for i=1:299 ]
stressImageBoxes = [ matplotlib.offsetbox.OffsetImage(stressPlots[i], zoom=0.5) for i=1:299 ]


fig,ax = PyPlot.subplots(2)

sc = ax[1].scatter( (2500/1352) .* data[1,1:299],data[3,1:299])
ax[1].set_xlabel("Density (kg/m3)")
ax[1].set_ylabel("Max Stress (kgf/cm2)")


ax[2].set_axis_off()

for i=1:299
    barImageBoxes[i].image.axes = ax[2]
    stressImageBoxes[i].image.axes = ax[2]
end

annotLeft = []
annotRight = []
for i=1:299
    ab1 = matplotlib.offsetbox.AnnotationBbox(barImageBoxes[i], [0.25,0.5])
    push!(annotLeft,ax[2].add_artist(ab1))
    ab2 = matplotlib.offsetbox.AnnotationBbox(stressImageBoxes[i], [0.75,0.5])
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
