# Data structures and methods

################################################################################
###############                   GEOMETRY                   ###################
################################################################################

mutable struct Point
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

################################################################################
###############                   MATERIAL                   ###################
################################################################################

struct Material
    isECC::Bool
    compressive::Float64
    tensile::Float64
    tensilePeak::Float64
    peakStrain::Float64
    crackStrainRatio::Float64
    Material(compressive,tensile) = new(true,compressive,tensile,0.040,0.045,0.8)
    Material(compressive,tensile,tensilePeak,peakStrain,crackStrainRatio) = new(true,compressive,tensile,tensilePeak,peakStrain,crackStrainRatio)
    Material(isECC,compressive,tensile,tensilePeak,peakStrain,crackStrainRatio) = new(isECC,compressive,tensile,tensilePeak,peakStrain,crackStrainRatio)
end

Base.copy(m::Material) = Material(m.isECC,m.compressive,m.tensile,m.tensilePeak,m.peakStrain,m.crackStrainRatio)

function randomMaterial()
    compressive = round(30+20*rand(),digits=2)
    tensile = round(3+2*rand(),digits=2)
    tensilePeak = round(0.02+0.02*rand(),digits=3)
    peakStrain = max(tensilePeak+0.005,round(0.02+0.03*rand(),digits=3))
    crackStrainRatio = round(0.7+0.25*rand(),digits=2)
    Material(compressive,tensile,tensilePeak,peakStrain,crackStrainRatio)
end

################################################################################
###############        METAMATERIAL FRAME : "SKELETON"       ###################
################################################################################

struct Skeleton
    size::Int # size of the pattern
    nodes::Array{Bool,2}
    links::Array{Bool,3}
end

Base.copy(s::Skeleton) = Skeleton(s.size,deepcopy(s.nodes),deepcopy(s.links))
Base.copy(t::Tuple{Skeleton,Material}) = (copy(t[1]),Material(t[2].isECC,t[2].compressive,t[2].tensile,t[2].tensilePeak,t[2].peakStrain,t[2].crackStrainRatio))

function randomSkeleton(H::Int; p=0.5)
    potLinks = trues(H,H,4)
    potLinks[:,1,1] = falses(H)
    potLinks[:,H,3:4] = falses(H,2)
    potLinks[H,:,:] = falses(H,4)
    nodes,links = randomGeometry(H,potLinks, p=p)
    return Skeleton(H,nodes,links)
end

function alteredSkeleton(skeleton::Skeleton, linkCoor, linkValue) # Alters an existing skeleton by changing the value of link at linkCoor
    nodes = falses(skeleton.size, skeleton.size)
    links = skeleton.links
    links[linkCoor] += linkValue
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

################################################################################
###############            MODEL FOR FEM SIMULATION          ###################
################################################################################

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
    forbiddenPoints = []
    for point in points
        restraint = collect("010")
        if !( -1e-9 < point.x < 12*(skeleton.size-1)+1e-9 ) || !( -1e-9 < point.z < 12*(skeleton.size-1)+1e-9 )
            push!(forbiddenPoints, point.n)
        elseif point.z == 0
            push!(loadPoints, [point.n,1])
            restraint[3] = '1'
        elseif point.z == 12*(skeleton.size-1)*unitSize
            push!(loadPoints, [point.n,-1])
            restraint[3] = '1'
        end
        if point.x == 12*(skeleton.size-1)*unitSize
            restraint[1] = '1'
        end
        if point.x == 0
            restraint[1] = '1'
        end
        boundaries[point.n] = join(restraint)
    end
    filter!(p -> !(p.n in forbiddenPoints), points)
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
    filter!(e -> length(findall(in(forbiddenPoints),e.points))==0, elements)
    weight = sum(skeleton.nodes .* nodeWeights) + sum(skeleton.links .* linkWeights)
    return Model(points, elements, loadPoints, boundaries, unitSize, skeleton.size, weight, material)
end

function fullScaleModelFromModel(model::Model) # Performs X and Y symmetry to result in 4 times the initial geometry
    basePoints = deepcopy(model.points)
    baseElements = deepcopy(model.elements)
    elements = deepcopy(baseElements)
    points = deepcopy(basePoints)
    currentPoint = maximum( (p -> p.n).(points) )+1
    currentElement = maximum( (e-> e.n).(elements) )+1
    size = model.size
    unitSize = model.unitSize
    loadPoints = []

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
    currentPoint = maximum( (p -> p.n).(points) )+1
    currentElement = maximum( (e-> e.n).(elements) )+1

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

    boundaries = [ "" for i in 1:maximum((p->p.n).(points)) ]

    for point in points
        restraint = collect("010")
        if point.z == 0
            push!(loadPoints, [point.n,1])
            restraint[3] = '1'
        elseif point.z == 24*(size-1)*unitSize
            push!(loadPoints, [point.n,-1])
            restraint[3] = '1'
        end
        boundaries[point.n] = join(restraint)
    end

    Model(points, elements, loadPoints, boundaries, unitSize, 2*size-1, model.weight, model.material)
end

function flatFullScaleModelFromModel(model::Model)

    basePoints = deepcopy(model.points)
    baseElements = deepcopy(model.elements)
    elements = deepcopy(baseElements)
    points = deepcopy(basePoints)
    size = model.size
    unitSize = model.unitSize
    loadPoints = []

    for step in 1:2
        baseElements = deepcopy(elements)
        basePoints = deepcopy(points)
        currentPoint = maximum( (p -> p.n).(points) )+1
        currentElement = maximum( (e-> e.n).(elements) )+1


        correspondance1 = [ basePoint.n for basePoint in basePoints ]
        correspondance2 = [ basePoint.n for basePoint in basePoints ]

        for i in 1:length(basePoints)
            if basePoints[i].x != 12*2^step*(size-1)*unitSize
                correspondance2[i] = currentPoint
                push!(points, Point(currentPoint, [12*2^step*(size-1)*unitSize-basePoints[i].x,basePoints[i].y,basePoints[i].z]))
                currentPoint += 1
            end
        end
        for baseElement in baseElements
            push!(elements, Element(currentElement, [ correspondance2[findall(x->x==baseElement.points[i],correspondance1)[1]] for i in 1:4 ]))
            currentElement += 1
        end
    end

    boundaries = [ "" for i in 1:maximum((p->p.n).(points)) ]

    for point in points
        restraint = collect("010")
        if point.z == 0
            push!(loadPoints, [point.n,1])
            restraint[3] = '1'
        elseif point.z == 24*unitSize
            push!(loadPoints, [point.n,-1])
            restraint[3] = '1'
        end
        if point.x == 0
            restraint[1] = '1'
        end
        boundaries[point.n] = join(restraint)
    end

    Model(points, elements, loadPoints, boundaries, unitSize, size, model.weight, model.material)
end

function zeroSymmetry(model::Model)
    basePoints = deepcopy(model.points)
    baseElements = deepcopy(model.elements)
    elements = deepcopy(baseElements)
    points = deepcopy(basePoints)
    size = model.size
    unitSize = model.unitSize
    loadPoints = []
    currentPoint = maximum( (p -> p.n).(points) )+1
    currentElement = maximum( (e-> e.n).(elements) )+1


    correspondance1 = [ basePoint.n for basePoint in basePoints ]
    correspondance2 = [ basePoint.n for basePoint in basePoints ]

    for i in 1:length(basePoints)
        if basePoints[i].x != 0
            correspondance2[i] = currentPoint
            push!(points, Point(currentPoint,[-basePoints[i].x,basePoints[i].y,basePoints[i].z]))
            currentPoint += 1
        end
    end
    for baseElement in baseElements
        push!(elements, Element(currentElement, [ correspondance2[findall(x->x==baseElement.points[i],correspondance1)[1]] for i in 1:4 ]))
        currentElement += 1
    end

    boundaries = [ "" for i in 1:maximum((p->p.n).(points)) ]

    for point in points
        restraint = collect("010")
        if point.z == 0
            push!(loadPoints, [point.n,1])
            restraint[3] = '1'
        elseif point.z == 24*unitSize
            push!(loadPoints, [point.n,-1])
            restraint[3] = '1'
        end
        boundaries[point.n] = join(restraint)
    end

    Model(points, elements, loadPoints, boundaries, unitSize, size, model.weight, model.material)
end


################################################################################
#############   SIMULATION : HOLDS DATA THROUGHOUT SIMULATION  #################
################################################################################


mutable struct Simulation
    filename::String
    model::Model
    dEpsilon::Float64 # Strain corresponding to each displacement step in the simulation
    step::Int # Current simulation step
    exit::Int # 0 : success, 1 : fail
    strain::Vector{Float64}
    stress::Vector{Float64}
end

function emptySimulation(iter,dEpsilon::Float64)
    return Simulation("metamat$iter",Model([],[],[],[],0,0,0,Material(0,0)),dEpsilon,0,1,[],[])
end

function runSimulation(simulation::Simulation;limit=0.2,output="none",angle=0,direction=1)

    # Files preparation

    filename = simulation.filename
    run(`mkdir $filename`)

    datFile = open("$filename.dat","w") # Initial start top of file (without steps)
    if output =="none"
        write(datFile, "Metamaterial Project\n1202000011000200000000002     0.700     0.000     0.000         0\nNODE\n")
    else
        write(datFile, "Metamaterial Project\n1202000011000000000000002     0.700     0.000     0.000         0\nNODE\n")
    end
    for point in simulation.model.points
        write(datFile, nodeLine(point,simulation.model.boundaries[point.n]))
    end
    write(datFile, "ELEM\n")
    for element in simulation.model.elements
        write(datFile,elementLine(element,simulation.model.material))
    end
    write(datFile, "LOAD\n\n")
    close(datFile)

    auxFile = open("$filename-restart.aux","w") # Restart function top of file (without steps)
    lines = readlines("$filename.dat",keep=true)
    for i=1:length(lines)
        if i == 2
            if output =="none"
                write(auxFile,"1212000011000200000000002     0.700     0.000     0.000         0\n")
            else
                write(auxFile,"1212000011000000000000002     0.700     0.000     0.000         0\n")
            end
        else
            write(auxFile,lines[i])
        end
    end
    close(auxFile)

    # COM3 Simulation

    runSteps(simulation,angle,direction)
    while simulation.strain[end] < limit && simulation.exit == 0 #&& simulation.stress[end]/maximum(simulation.stress) > 0.25 && ( simulation.stress[end]/maximum(simulation.stress) > 0.4 || std(simulation.stress[end-5:end]) > 0.35 )
        runSteps(simulation,angle,direction)
    end

    # Cleanup

    mechFiles = glob("$(filename)/MECHIFI*")
    run(`rm $mechFiles`)
    addSteps(simulation.model.loadPoints, 0, simulation.dEpsilon,
    simulation.model.unitSize, simulation.model.size, simulation.filename, restart = true, nSteps = simulation.step,angle=angle,direction=direction) # Creates the full .dat file with all steps for later re-running purposes
    run(`mv $filename.dat $filename`)
    run(`mv $filename/aux.nod $filename/$filename-MECH.nod`)
    run(`mv $filename/aux.crk $filename/$filename-MECH.crk`)
    run(`mv $filename/aux.fld $filename/$filename-MECH.fld`)

    if simulation.exit == 1
        run(`rm -r $filename`)
    end

end
