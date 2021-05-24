# Tools for generating / Plotting geometries


function makeNodeWeights(H,unitSize)
    nodeWeights = [4*ones(1,H);8*ones(H-2,H);4*ones(1,H)]
    nodeWeights[:,1] = 4*ones(H); nodeWeights[:,end] = 4*ones(H)
    nodeWeights[1,1] = 2; nodeWeights[1,end] = 2; nodeWeights[end,1] = 2; nodeWeights[end,end] = 2
    nodeWeights = unitSize^2 .* nodeWeights
end

function makeLinkWeights(H,unitSize)
    linkWeights = zeros(H,H,4)
    linkWeights[:,:,1] = (20*sqrt(2)-2)*unitSize^2 .* ones(H,H); linkWeights[:,:,3] = (20*sqrt(2)-2)*unitSize^2 .* ones(H,H)
    linkWeights[:,:,2] = (20-4*sqrt(2))*unitSize^2 .* ones(H,H); linkWeights[:,:,4] = (20-4*sqrt(2))*unitSize^2 .* ones(H,H)
    linkWeights
end

function makeBasePoints()
    basePointsData = readdlm("Ressources/points.csv",',')
    basePoints = [ Point(Int(basePointsData[i,1]),[basePointsData[i,2],0.0,basePointsData[i,3]]) for i in 1:length(basePointsData[:,1])]
end

function makeBaseElements(H)
    baseElementsData = readdlm("Ressources/elements.csv",',',Int)
    for index in findall(x -> x > 1000, baseElementsData)
        if baseElementsData[index] < 2000
            baseElementsData[index] = baseElementsData[index]%1000 + (H-1)*145
        elseif baseElementsData[index] < 3000
            baseElementsData[index] = baseElementsData[index]%2000 + H*145
        elseif baseElementsData[index] < 4000
            baseElementsData[index] = baseElementsData[index]%3000 + (H+1)*145
        else
            baseElementsData[index] = baseElementsData[index]%4000 + 145
        end
    end
    baseElements = [ Element(i,baseElementsData[i,2:end]) for i in baseElementsData[:,1] ]
end


function randomGeometry(H,nodes,links,potLinks; p=0.7)
    for i in 1:H-1
        i == 1 ? nodes[1,:] = rand(H) .< p : nothing
        nodes[1,1] = true
        for j in 1:H
            if nodes[i,j]
                potLinksCount = length(findall(potLinks[i,j,:]))
                for l in 1:rand(1:potLinksCount)
                    index = rand(findall(potLinks[i,j,:]))
                    links[i,j,index] = true; potLinks[i,j,index] = false
                    if index == 1
                        nodes[i+1,j-1] = true
                    elseif index == 2
                        nodes[i+1,j] = true
                    elseif index == 3
                        nodes[i+1,j+1] = true
                        potLinks[i,j+1,1] = false
                    else
                        nodes[i,j+1] = true
                    end
                end
            end
        end
    end

    for j in 1:H-1
        if nodes[H,j] && rand(Bool)
            links[H,j,4] = true
            nodes[H,j+1] = true
        end
    end
end

function plotGeometry(skeleton::Skeleton)
    barPlot = []
    H = skeleton.size
    for i in 1:H, j in 1:H
        for index in findall(skeleton.links[i,j,:])
            if index == 1
                push!(barPlot, [[j, j-1], [i,i+1]])
                push!(barPlot, [[j, j-1], [2*(H)-i,2*(H)-(i+1)]])
                push!(barPlot, [[2*(H)-j, 2*(H)-(j-1)], [i,(i+1)]])
                push!(barPlot, [[2*(H)-j, 2*(H)-(j-1)], [2*(H)-i,2*(H)-(i+1)]])
            elseif index == 2
                push!(barPlot, [[j, j], [i,i+1]])
                push!(barPlot, [[j, j], [2*(H)-i,2*(H)-(i+1)]])
                push!(barPlot, [[2*(H)-j, 2*(H)-j], [i,i+1]])
                push!(barPlot, [[2*(H)-j, 2*(H)-j], [2*(H)-i,2*(H)-(i+1)]])
            elseif index == 3
                push!(barPlot, [[j, j+1], [i,i+1]])
                push!(barPlot, [[j, j+1], [2*(H)-i,2*(H)-(i+1)]])
                push!(barPlot, [[2*(H)-j, 2*(H)-(j+1)], [i,i+1]])
                push!(barPlot, [[2*(H)-j, 2*(H)-(j+1)], [2*(H)-i,2*(H)-(i+1)]])
            else
                push!(barPlot, [[j, j+1], [i,i]])
                push!(barPlot, [[j, j+1], [2*(H)-i,2*(H)-i]])
                push!(barPlot, [[2*(H)-j, 2*(H)-(j+1)], [i,i]])
                push!(barPlot, [[2*(H)-j, 2*(H)-(j+1)], [2*(H)-i,2*(H)-i]])
            end
        end
    end
    barPlot
end

function plotGeometry(model::Model)
    H = model.size
    polygons = []
    points = model.points
    for element in model.elements
        P = [ points[findall(p->p.n==element.points[i],points)[1]] for i in 1:4 ]
        X = (p->p.x).(P); Y = (p->p.z).(P)
        push!(polygons, Shape(X, Y))
        push!(polygons, Shape(24*(H-1).-X, Y))
        push!(polygons, Shape(X, 24*(H-1).-Y))
        push!(polygons, Shape(24*(H-1).-X, 24*(H-1).-Y))
    end
    polygons
end

function plotModel(model)
    polygons = []
    points = model.points
    for element in model.elements
        P = [ points[findall(p->p.n==element.points[i],points)[1]] for i in 1:4 ]
        push!(polygons, Shape((p->p.x).(P), (p->p.z).(P)))
    end
    plt = plot(size=(800,800)) #axis = nothing
    for polygon in polygons
        plot!(plt,polygon,label=false,color=:grey)
    end
    png(plt,"modelPlot.png")
end
