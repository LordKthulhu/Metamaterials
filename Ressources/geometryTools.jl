# Tools for generating / Plotting geometries

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

function plotGeometry(barPlots,iter,links,H)
    for i in 1:H, j in 1:H
        for index in findall(links[i,j,:])
            if index == 1
                push!(barPlots[iter], [[j, j-1], [i,i+1]])
                push!(barPlots[iter], [[j, j-1], [2*(H)-i,2*(H)-(i+1)]])
                push!(barPlots[iter], [[2*(H)-j, 2*(H)-(j-1)], [i,(i+1)]])
                push!(barPlots[iter], [[2*(H)-j, 2*(H)-(j-1)], [2*(H)-i,2*(H)-(i+1)]])
            elseif index == 2
                push!(barPlots[iter], [[j, j], [i,i+1]])
                push!(barPlots[iter], [[j, j], [2*(H)-i,2*(H)-(i+1)]])
                push!(barPlots[iter], [[2*(H)-j, 2*(H)-j], [i,i+1]])
                push!(barPlots[iter], [[2*(H)-j, 2*(H)-j], [2*(H)-i,2*(H)-(i+1)]])
            elseif index == 3
                push!(barPlots[iter], [[j, j+1], [i,i+1]])
                push!(barPlots[iter], [[j, j+1], [2*(H)-i,2*(H)-(i+1)]])
                push!(barPlots[iter], [[2*(H)-j, 2*(H)-(j+1)], [i,i+1]])
                push!(barPlots[iter], [[2*(H)-j, 2*(H)-(j+1)], [2*(H)-i,2*(H)-(i+1)]])
            else
                push!(barPlots[iter], [[j, j+1], [i,i]])
                push!(barPlots[iter], [[j, j+1], [2*(H)-i,2*(H)-i]])
                push!(barPlots[iter], [[2*(H)-j, 2*(H)-(j+1)], [i,i]])
                push!(barPlots[iter], [[2*(H)-j, 2*(H)-(j+1)], [2*(H)-i,2*(H)-i]])
            end
        end
    end
end

function writeNodes(nodes, links, datFile, basePoints, loadPoints, H, unitSize)
    for index in findall(nodes)
        i = index[1]; j = index[2]
        nodeNumber = j + H*(i-1)
        startingPoint = 145*(nodeNumber-1)
        for pointNumber in 1:37
            write(datFile,nodeLine( startingPoint + pointNumber, unitSize * (12*(j-1) + basePoints[pointNumber,2]), 0, unitSize * (12*(i-1) + basePoints[pointNumber,3]), loadPoints ))
        end
        for link in findall(links[i, j,:])
            if link == 1
                for pointNumber in 38:70
                    write(datFile,nodeLine( startingPoint + pointNumber, unitSize * (12*(j-1) + basePoints[pointNumber,2]), 0, unitSize * (12*(i-1) + basePoints[pointNumber,3]), loadPoints ))
                end
            elseif link == 2
                for pointNumber in 71:91
                    write(datFile,nodeLine( startingPoint + pointNumber, unitSize * (12*(j-1) + basePoints[pointNumber,2]), 0, unitSize * (12*(i-1) + basePoints[pointNumber,3]), loadPoints ))
                end
            elseif link == 3
                for pointNumber in 92:124
                    write(datFile,nodeLine( startingPoint + pointNumber, unitSize * (12*(j-1) + basePoints[pointNumber,2]), 0, unitSize * (12*(i-1) + basePoints[pointNumber,3]), loadPoints ))
                end
            else
                for pointNumber in 125:145
                    write(datFile,nodeLine( startingPoint + pointNumber, unitSize * (12*(j-1) + basePoints[pointNumber,2]), 0, unitSize * (12*(i-1) + basePoints[pointNumber,3]), loadPoints ))
                end
            end
        end
    end
end

function writeElements(nodes, links, datFile, baseElements, H)
    for index in findall(nodes)
        nodeNumber = index[2] + H*(index[1]-1)
        startingPoint = 145*(nodeNumber-1)
        startingElement = 112*(nodeNumber-1)
        for elementNumber in 1:32
            write(datFile,elementLine( startingElement + elementNumber, startingPoint .+ baseElements[elementNumber,2:end] ))
        end
        for link in findall(links[index[1], index[2],:])
            if link == 1
                for elementNumber in 33:56
                    write(datFile,elementLine( startingElement + elementNumber, startingPoint .+ baseElements[elementNumber,2:end] ))

                end
            elseif link == 2
                for elementNumber in 57:72
                    write(datFile,elementLine( startingElement + elementNumber, startingPoint .+ baseElements[elementNumber,2:end] ))
                end
            elseif link == 3
                for elementNumber in 73:96
                    write(datFile,elementLine( startingElement + elementNumber, startingPoint .+ baseElements[elementNumber,2:end] ))
                end
            else
                for elementNumber in 97:112
                    write(datFile,elementLine( startingElement + elementNumber, startingPoint .+ baseElements[elementNumber,2:end] ))
                end
            end
        end
    end
end

function writeElementsPlain(nodes, links, datFile, baseElements, H)
    for index in findall(nodes)
        nodeNumber = index[2] + H*(index[1]-1)
        startingPoint = 145*(nodeNumber-1)
        startingElement = 112*(nodeNumber-1)
        for elementNumber in 1:32
            write(datFile,elementLinePlain( startingElement + elementNumber, startingPoint .+ baseElements[elementNumber,2:end] ))
        end
        for link in findall(links[index[1], index[2],:])
            if link == 1
                for elementNumber in 33:56
                    write(datFile,elementLinePlain( startingElement + elementNumber, startingPoint .+ baseElements[elementNumber,2:end] ))

                end
            elseif link == 2
                for elementNumber in 57:72
                    write(datFile,elementLinePlain( startingElement + elementNumber, startingPoint .+ baseElements[elementNumber,2:end] ))
                end
            elseif link == 3
                for elementNumber in 73:96
                    write(datFile,elementLinePlain( startingElement + elementNumber, startingPoint .+ baseElements[elementNumber,2:end] ))
                end
            else
                for elementNumber in 97:112
                    write(datFile,elementLinePlain( startingElement + elementNumber, startingPoint .+ baseElements[elementNumber,2:end] ))
                end
            end
        end
    end
end
