using Glob
using Plots

H=3

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

steps = [ [] for i in 1:5 ]

i=3

for stepn in 1:5
    while i <= length(lines) && lines[i][2:5] != "STEP"
        values = split(lines[i])
        push!(steps[stepn],[parse(Float64,values[2]),parse(Float64,values[4])])
        global i+=1
    end
    global i+=2
end

for stepn in 1:5
    for i in 1:length(points)
        steps[stepn][i][1] += points[i][2]
        steps[stepn][i][2] += points[i][3]
    end
end

plt = plot(size=(800,800),axis=nothing,xlims=(-10,58),ylims=(-10,58))
polygons = []
for element in elements
    P = [ points[findall(p->p[1]==element[i],points)[1]] for i in 1:4 ]
    X = (p->p[2]).(P); Y = (p->p[3]).(P)
    push!(polygons, Shape(X, Y))
    push!(polygons, Shape(24*(H-1).-X, Y))
    push!(polygons, Shape(X, 24*(H-1).-Y))
    push!(polygons, Shape(24*(H-1).-X, 24*(H-1).-Y))
end
for polygon in polygons
    plot!(plt,polygon,label=false,color=:grey)
end
png(plt,"step0.png")

for stepn in 1:5
    plt = plot(size=(800,800),axis=nothing,xlims=(-10,58),ylims=(-10,58))
    polygons = []
    for element in elements
        P = [ steps[stepn][findall(p->p[1]==element[i],points)[1]] for i in 1:4 ]
        X = (p->p[1]).(P); Y = (p->p[2]).(P)
        push!(polygons, Shape(X, Y))
        push!(polygons, Shape(24*(H-1).-X, Y))
        push!(polygons, Shape(X, 24*(H-1).-Y))
        push!(polygons, Shape(24*(H-1).-X, 24*(H-1).-Y))
    end
    for polygon in polygons
        plot!(plt,polygon,label=false,color=:grey)
    end
    png(plt,"step$stepn.png")
end
