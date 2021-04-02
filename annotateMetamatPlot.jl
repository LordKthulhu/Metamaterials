using PyCall
using PyPlot
using DelimitedFiles

function hover(event)
    visLeft = [ a.get_visible() for a in annotLeft ]
    visRight = [ a.get_visible() for a in annotRight ]
    if event.inaxes == ax[2]
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


barModels = [ PyPlot.imread("metamat$i/metamat$i-barplot.png",format="png") for i=1:500 ]
barImageBoxes = [ matplotlib.offsetbox.OffsetImage(barModels[i], zoom=0.3) for i=1:500 ]

stressPlots = [ PyPlot.imread("metamat$i/metamat$i-plot.png",format="png") for i=1:500 ]
stressImageBoxes = [ matplotlib.offsetbox.OffsetImage(stressPlots[i], zoom=0.3) for i=1:500 ]


fig,ax = PyPlot.subplots(1,3)

sc = ax[2].scatter(data[1,:],data[2,:])

ax[1].set_axis_off()
ax[3].set_axis_off()

for i=1:500
    barImageBoxes[i].image.axes = ax[1]
    stressImageBoxes[i].image.axes = ax[3]
end

annotLeft = []
annotRight = []
for i=1:500
    ab1 = matplotlib.offsetbox.AnnotationBbox(imagebox, [0.5,0.5])
    push!(annot,ax[2].add_artist(ab1))
    ab2 = matplotlib.offsetbox.AnnotationBbox(imagebox, [0.5,0.5])
    push!(annot,ax[2].add_artist(ab2))

for a in annotLeft
    a.set_visible(false)
end
for a in annotRight
    a.set_visible(false)
end


fig.canvas.mpl_connect("motion_notify_event", hover)

PyPlot.show()
