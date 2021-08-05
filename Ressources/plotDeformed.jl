using Glob
using Plots
using LinearAlgebra
using Colors
using LaTeXStrings

H=3

const defFactor = 1

const cmap2 =  reverse(colormap("Reds",101))
const cmap1 = reverse(colormap("Blues",101))

function colorElem(stress,sigmin,sigmax)
    if stress > 0
        cmap2[1+Int(round(100*(sigmax-stress)/sigmax))]
    else
        cmap1[1+Int(round(100*(sigmin-stress)/sigmin))]
    end
end

datFile = glob("*.dat")[1]
lines = readlines(datFile)

points = []
elements = []

i=4

while lines[i][1:4] == "NODE"
    values = split(lines[i])
    push!(points,[parse(Int,values[2]),parse(Float64,values[3]),parse(Float64,values[5])])
    global i+=1
end
i+=1
while lines[i][1:4] == "PLAT"
    values = split(lines[i])
    push!(elements,[parse(Int,values[3]),parse(Int,values[4]),parse(Int,values[5]),parse(Int,values[6])])
    global i+=4
end

nodFile = glob("*.nod")[1]
lines = readlines(nodFile)

Nsteps = length(findall(l->l[2:5]=="STEP",lines))

steps = [ [] for i in 1:Nsteps ]

i=3

for stepn in 1:Nsteps
    while i <= length(lines) && lines[i][2:5] != "STEP"
        values = split(lines[i])
        push!(steps[stepn],defFactor .* [parse(Float64,values[2]),parse(Float64,values[4])])
        global i+=1
    end
    global i+=2
end

fldFile = glob("*.fld")[1]
lines = readlines(fldFile)

i=3
stress = [[] for i in 1:Nsteps]

for stepn in 1:Nsteps
    while i <= length(lines) && lines[i][2:5] != "STEP"
        sigxx = 0; sigyy = 0; sigxy =0
        for j in 0:3
            values = split(lines[i+j])
            sigxx += parse(Float64,values[6])/4
            sigyy += parse(Float64,values[7])/4
            sigxy += parse(Float64,values[8])/4
        end
        push!(stress[stepn], eigvals([sigxx sigxy; sigxy sigyy]))
        global i+=4
    end
    global i+=2
end

sigmax = maximum([ stress[i][j][2] for j=1:length(elements),i=1:Nsteps ])
sigmin = minimum([ stress[i][j][1] for j=1:length(elements),i=1:Nsteps ])
colors = [[] for i in 1:Nsteps ]
for stepn in 1:Nsteps
    for elem in stress[stepn]
        push!(colors[stepn],[colorElem(elem[1],sigmin,sigmax),colorElem(elem[2],sigmin,sigmax)])
    end
end


for stepn in 1:Nsteps
    for i in 1:length(points)
        steps[stepn][i][1] += points[i][2]
        steps[stepn][i][2] += points[i][3]
    end
end


pltColor1 = plot(xlabel="kgf/cm2",yaxis=([], false),margin = 5Plots.mm,border=false)
pltColor2 = plot(xlabel="kgf/cm2",yaxis=([], false),margin = 5Plots.mm,border=false)
x1 = reverse([ sigmin/101*i for i in 0:101 ])
x2 = reverse([ sigmax/101*i for i in 0:101 ])
for i=1:101
    plot!(pltColor1, Shape([x1[i],x1[i+1],x1[i+1],x1[i]],[0,0,1,1]),color=cmap1[i],label=false)
    plot!(pltColor2, Shape([x2[i],x2[i+1],x2[i+1],x2[i]],[0,0,1,1]),color=cmap2[i],label=false)
end

l = [ @layout [ a b; c{0.05h} d{0.05h}] for i in 1:Nsteps]

# plt0 = plot(size=(800,800),axis=nothing,xlims=(-10,58),ylims=(-10,58))
# polygons0 = []
# for element in elements
#     P = [ points[findall(p->p[1]==element[i],points)[1]] for i in 1:4 ]
#     X = (p->p[2]).(P); Y = (p->p[3]).(P)
#     push!(polygons0, Shape(X, Y))
#     push!(polygons0, Shape(24*(H-1).-X, Y))
#     push!(polygons0, Shape(X, 24*(H-1).-Y))
#     push!(polygons0, Shape(24*(H-1).-X, 24*(H-1).-Y))
# end
# for polygon in polygons0
#     plot!(plt0,polygon,label=false,color=:grey)
# end
# png(plt0,"step000.png")

for stepn in 1:2:Nsteps
    plt1 = plot(size=(800,800),axis=([], false),xlims=(-10,58),ylims=(-10,58),title=L"\sigma_{I}",border=false)
    plt2 = plot(size=(800,800),axis=([], false),xlims=(-10,58),ylims=(-10,58),title=L"\sigma_{II}",border=false)
    Xmax = maximum(getindex.(steps[stepn],1)); Ymax = maximum(getindex.(steps[stepn],2))
    polygons = []
    for element in elements
        P = [ steps[stepn][findall(p->p[1]==element[i],points)[1]] for i in 1:4 ]
        X = (p->p[1]).(P); Y = (p->p[2]).(P)
        push!(polygons, Shape(X, Y))
        push!(polygons, Shape(2*Xmax.-X, Y))
        push!(polygons, Shape(X, 2*Ymax.-Y))
        push!(polygons, Shape(2*Xmax.-X, 2*Ymax.-Y))
    end
    for i in 1:length(elements)
        plot!(plt1,polygons[4*i-3],label=false,color=colors[stepn][i][1])
        plot!(plt1,polygons[4*i-2],label=false,color=colors[stepn][i][1])
        plot!(plt1,polygons[4*i-1],label=false,color=colors[stepn][i][1])
        plot!(plt1,polygons[4*i],label=false,color=colors[stepn][i][1])
        plot!(plt2,polygons[4*i-3],label=false,color=colors[stepn][i][2])
        plot!(plt2,polygons[4*i-2],label=false,color=colors[stepn][i][2])
        plot!(plt2,polygons[4*i-1],label=false,color=colors[stepn][i][2])
        plot!(plt2,polygons[4*i],label=false,color=colors[stepn][i][2])
    end
    plt = plot(plt1,plt2,pltColor1,pltColor2,layout=l[stepn],size = (1600,900))
    # plt = plot(plt1,plt2,layout=(1,2),size = (1600,840))
    png(plt,"step$(lpad(stepn,3,"0")).png")
end
