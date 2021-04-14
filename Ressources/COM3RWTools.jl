### Tools for Reading from / Writing to COM3 files


################################################################################
#####                           WRITING TOOLS                              #####
################################################################################

function nodeLine(point::Point,restraint::String)
    sn = string(point.n)
    sx = @sprintf("%.3f",point.x); sy = @sprintf("%.3f",point.y); sz = @sprintf("%.3f",point.z)
    "NODE " * " "^(5-length(sn)) * sn * " "^(10-length(sx)) * sx * " "^(10-length(sy)) * sy * " "^(10-length(sz)) * sz * "                    " * restraint * "\n"
end

function elementLine(element::Element, material::Material)

    sn = string(element.n); P = string.(element.points); sTensile = @sprintf("%.2f",10*material.tensile); sCompressive = @sprintf("%.2f",10*material.compressive)
    material.isECC ? mat="3" : mat="0"

    "PLAT " *" "^(5-length(sn)) * sn * " "^(5-length(P[1])) * P[1] * " "^(5-length(P[2])) * P[2] * " "^(5-length(P[3])) * P[3] *
    " "^(5-length(P[4])) * P[4] * "    0    0    0    0    0    0    0    0    0\n" *
    "    0    0    0    0    0    0    0    "*mat*"       0.0       1.0    3000.0 0.0025000       0.0                                                                                                         0.005     0.046       0.8\n" *
    " "^(10-length(sCompressive)) * sCompressive * " "^(10-length(sTensile)) * sTensile * "      0.17    0\n" *
    " 2100000.0    4000.0   0.00000 2100000.0    4000.0   0.00000       0.0       0.0       0.0\n"
end

function loadLine(n,F)
    sn = string(n)
    sf = @sprintf("%.5f",F)
    "LOAD " * " "^(5-length(sn)) * sn * "   0.00000   0.00000" * " "^(10-length(sf)) * sf * "\n"
end

function addSteps(loadPoints, step, dEpsilonMax, unitSize, filename; restart = true, nSteps = 5)
    restart ? run(`cp $(filename)-restart.aux $(filename).dat`) : nothing
    datFile = open("$filename.dat","a")
    for i=1:nSteps
        timeStr = @sprintf("%.4f",0.001*(step+i))
        write(datFile, "STEP " * " "^(5-length(string(step+i)))*string(step+i) * " "^(10-length(timeStr))*timeStr * "     0.000     0.000     0.000                                             0.000\n")
        for load in loadPoints
            write(datFile, loadLine(load[1], 12*(H-1)*unitSize*dEpsilonMax * load[2]))
        end
    end
    write(datFile, "END\n")
    close(datFile)
end

################################################################################
#####                           READING TOOLS                              #####
################################################################################

function forceSteps(filename)
    force = 0
    forceList = []
    for line in readlines(filename * "/" * filename * "-MECH.nod",keep=true)[3:end]
        if line[2:5] == "STEP"
            push!(forceList,force)
        elseif line[2:5] == "NODE"
            force = 0
        else
            force += abs( parse( Float64, split(line)[7] ))
        end
    end
    push!(forceList,force)
    forceList
end
