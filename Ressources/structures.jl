# Data structures and methods

struct Point
    n::Int
    x::Float64
    y::Float64
    z::Float64
    Point(n,coordinates::Vector{Float64}) = length(coordinates) == 3 ? new(n,coordinates[1],coordinates[2],coordinates[3]) : error("Point should have 3 coordinates.")
end

struct Element
    n::Int
    points::Vector{Int}
    Element(n,points) = length(points) == 4 ? new(n,points) : error("Wrong number of points provided : should be 4.")
end

struct Material
    isECC::Bool
    compressive::Float64
    tensile::Float64
end

struct Skeleton
    size::Int
    nodes::Array{Bool,2}
    links::Array{Bool,3}
end

Base.copy(s::Skeleton) = Skeleton(s.size,s.nodes,s.links)

function randomSkeleton(H::Int; p=0.3)
    potLinks = trues(H,H,4)
    nodes = falses(H, H)
    links = falses(H, H, 4)
    potLinks[:,1,1] = falses(H)
    potLinks[:,H,3:4] = falses(H,2)
    potLinks[H,:,:] = falses(H,4)
    randomGeometry(H,nodes,links,potLinks, p=p)
    return Skeleton(H,nodes,links)
end

function alteredSkeleton(skeleton::Skeleton, linkCoor, linkValue)
    nodes = falses(skeleton.size, skeleton.size)
    links = skeleton.links
    links[linkCoor] = linkValue
    for link in findall(links)
        if link[1]!=0 && link[2]!=0
            nodes[link[1],link[2]] = true
        end
        if link[3]==1
            nodes[link[1]+1,link[2]-1] = true
        elseif link[3]==2
            nodes[link[1]+1,link[2]] = true
        elseif link[3]==3
            nodes[link[1]+1,link[2]+1] = true
        else
            nodes[link[1],link[2]+1] = true
        end
    end
    Skeleton(skeleton.size,nodes,links)
end

struct Model
    points::Vector{Point}
    elements::Vector{Element}
    loadPoints::Vector{Vector{Int}}
    boundaries::Vector{String}
    unitSize::Float64
    size::Int
    weight::Float64
    material::Material
end


function modelFromSkeleton(skeleton::Skeleton,material::Material,unitSize, basePoints, baseElements, nodeWeights, linkWeights)
    points = []
    for index in findall(skeleton.nodes)
        i = index[1]; j = index[2]
        nodeNumber = j + skeleton.size*(i-1)
        startingPoint = 145*(nodeNumber-1)
        for pointNumber in 1:37
            push!(points,Point(startingPoint + pointNumber, [unitSize * (12*(j-1) + basePoints[pointNumber].x), 0.0, unitSize * (12*(i-1) + basePoints[pointNumber].z)]))
        end
        for link in findall(skeleton.links[i, j,:])
            if link == 1
                for pointNumber in 38:70
                    push!(points,Point( startingPoint + pointNumber, [unitSize * (12*(j-1) + basePoints[pointNumber].x), 0.0, unitSize * (12*(i-1) + basePoints[pointNumber].z)]))
                end
            elseif link == 2
                for pointNumber in 71:91
                    push!(points,Point( startingPoint + pointNumber, [unitSize * (12*(j-1) + basePoints[pointNumber].x), 0.0, unitSize * (12*(i-1) + basePoints[pointNumber].z)]))
                end
            elseif link == 3
                for pointNumber in 92:124
                    push!(points,Point( startingPoint + pointNumber, [unitSize * (12*(j-1) + basePoints[pointNumber].x), 0.0, unitSize * (12*(i-1) + basePoints[pointNumber].z)]))
                end
            else
                for pointNumber in 125:145
                    push!(points,Point( startingPoint + pointNumber, [unitSize * (12*(j-1) + basePoints[pointNumber].x), 0.0, unitSize * (12*(i-1) + basePoints[pointNumber].z)]))
                end
            end
        end
    end
    boundaries = ["" for i in 1:(145*skeleton.size^2)]
    loadPoints = []
    for point in points
        restraint = collect("010")
        if point.z == 0
            push!(loadPoints, [point.n,-1])
            restraint[3] = '1'
        elseif point.z == 12*(skeleton.size-1)*unitSize
            push!(loadPoints, [point.n,1])
            restraint[3] = '1'
        elseif point.x == 12*(skeleton.size-1)*unitSize
            restraint[1] = '1'
        end
        boundaries[point.n] = join(restraint)
    end
    elements = []
    for index in findall(skeleton.nodes)
        nodeNumber = index[2] + skeleton.size*(index[1]-1)
        startingPoint = 145*(nodeNumber-1)
        startingElement = 112*(nodeNumber-1)
        for elementNumber in 1:32
            push!(elements,Element( startingElement + elementNumber, startingPoint .+ baseElements[elementNumber].points))
        end
        for link in findall(skeleton.links[index[1], index[2],:])
            if link == 1
                for elementNumber in 33:56
                    push!(elements,Element( startingElement + elementNumber, startingPoint .+ baseElements[elementNumber].points))
                end
            elseif link == 2
                for elementNumber in 57:72
                    push!(elements,Element( startingElement + elementNumber, startingPoint .+ baseElements[elementNumber].points))
                end
            elseif link == 3
                for elementNumber in 73:96
                    push!(elements,Element( startingElement + elementNumber, startingPoint .+ baseElements[elementNumber].points))
                end
            else
                for elementNumber in 97:112
                    push!(elements,Element( startingElement + elementNumber, startingPoint .+ baseElements[elementNumber].points))
                end
            end
        end
    end
    weight = sum(skeleton.nodes .* nodeWeights) + sum(skeleton.links .* linkWeights)
    return Model(points, elements, loadPoints, boundaries, unitSize, skeleton.size, weight, material)
end

function fullScaleModelFromModel(model::Model)
    basePoints = model.points
    baseElements = model.elements
    elements = deepcopy(baseElements)
    points = deepcopy(basePoints)
    currentPoint = basePoints[end].n+1
    currentElement = baseElements[end].n+1
    size = model.size
    unitSize = model.unitSize

    correspondance1 = [ basePoint.n for basePoint in basePoints ]
    correspondance2 = [ basePoint.n for basePoint in basePoints ]

    for i in 1:length(basePoints)
        if basePoints[i].x != 12*(size-1)*unitSize
            correspondance2[i] = currentPoint
            push!(points, Point(currentPoint, [24*(size-1)*unitSize-basePoints[i].x,basePoints[i].y,basePoints[i].z]))
            currentPoint += 1
        end
    end
    for baseElement in baseElements
        push!(elements, Element(currentElement, [ correspondance2[findall(x->x==baseElement.points[i],correspondance1)[1]] for i in 1:4 ]))
        currentElement += 1
    end

    basePoints = deepcopy(points)
    baseElements = deepcopy(elements)

    correspondance1 = [ basePoint.n for basePoint in basePoints ]
    correspondance2 = [ basePoint.n for basePoint in basePoints ]
    currentPoint = basePoints[end].n+1
    currentElement = baseElements[end].n+1

    for i in 1:length(basePoints)
        if basePoints[i].z != 12*(size-1)*unitSize
            correspondance2[i] = currentPoint
            push!(points, Point(currentPoint, [basePoints[i].x,basePoints[i].y,24*(size-1)*unitSize-basePoints[i].z]))
            currentPoint += 1
        end
    end
    for baseElement in baseElements
        push!(elements, Element(currentElement, [ correspondance2[findall(x->x==baseElement.points[i],correspondance1)[1]] for i in 1:4 ]))
        currentElement += 1
    end
    Model(points,elements,model.loadPoints, model.boundaries,unitSize,2*size-1,model.weight,model.material)
end

mutable struct Simulation
    filename::String
    model::Model
    dEpsilon::Float64
    step::Int
    exit::Int
    strain::Vector{Float64}
    stress::Vector{Float64}
end

function emptySimulation(iter,dEpsilon::Float64)
    return Simulation("metamat$iter",Model([],[],[],[],0,0,0,Material(true,0,0)),dEpsilon,0,1,[],[])
end

function runSimulation(simulation::Simulation)
    filename = simulation.filename
    run(`mkdir $filename`)
    datFile = open("$filename.dat","w")
    write(datFile, "Metamaterial Project\n1202000001000200000000002     0.700     0.000     0.000         0\nNODE\n")
    for point in simulation.model.points
        write(datFile, nodeLine(point,simulation.model.boundaries[point.n]))
    end
    write(datFile, "ELEM\n")
    for element in simulation.model.elements
        write(datFile,elementLine(element,simulation.model.material))
    end
    write(datFile, "LOAD\n\n")
    close(datFile)

    auxFile = open("$filename-restart.aux","w")
    lines = readlines("$filename.dat",keep=true)
    for i=1:length(lines)
        if i == 2
            write(auxFile,"1212000001000200000000002     0.700     0.000     0.000         0\n")
        else
            write(auxFile,lines[i])
        end
    end
    close(auxFile)
    runSteps(simulation)
    while simulation.stress[end]/maximum(simulation.stress) > 0.5 && simulation.exit == 0 && simulation.step <= 150
        runSteps(simulation)
    end

    mechFiles = glob("$(filename)/MECHIFI*")
    run(`rm $mechFiles`)

    if simulation.step > 150
        simulation.exit = 1
    end

    if simulation.exit == 1
        #println("Simulation $filename failed. Starting over.")
        run(`rm -r $filename`)
    end
end
