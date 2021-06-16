
mutable struct MetamatEnv2 <: AbstractEnvironment
    state::Tuple{Skeleton,Material}
    reward::Float64
    maxsteps::Int
    index::Int
end


MetamatEnv2(;H=3, maxsteps = 20, index = 1) = MetamatEnv2((randomSkeleton(H),randomMaterial()), 1e-9, maxsteps, index)

reset!(env::MetamatEnv2) = (env.state = (randomSkeleton(env.state[1].size), randomMaterial()); env.reward = 1e-9; env)

function actions(env::MetamatEnv2, s)
    H = s[1].size
    potLinks = trues(H,H,4)
    potLinks[:,1,1] = falses(H)
    potLinks[:,H,3:4] = falses(H,2)
    potLinks[H,:,1:3] = falses(H,3)
    for index in findall(s[1].links[1:end-1,1:end-1,3])
        potLinks[index[1],index[2]+1,1] = false
    end
    for index in findall(s[1].links[1:end-1,2:end,1])
        potLinks[index[1],index[2],3] = false
    end

    mat = s[2]
    potMat = []
    mat.compressive > 20 ? push!(potMat, [-1,0,0,0,0]) : nothing
    mat.compressive < 60 ? push!(potMat, [1,0,0,0,0]) : nothing
    mat.tensile > 2 ? push!(potMat, [0,-1,0,0,0]) : nothing
    mat.tensile < 6 ? push!(potMat, [0,1,0,0,0]) : nothing
    mat.tensilePeak > 0.002 ? push!(potMat, [0,0,-1,0,0]) : nothing
    mat.tensilePeak < 0.006 ? push!(potMat, [0,0,1,0,0]) : nothing
    mat.peakStrain > 0.02 ? push!(potMat, [0,0,0,-1,0]) : nothing
    mat.peakStrain < 0.06 ? push!(potMat, [0,0,0,1,0]) : nothing
    mat.crackStrainRatio > 0.7 ? push!(potMat, [0,0,0,0,-1]) : nothing
    mat.crackStrainRatio < 0.95 ? push!(potMat, [0,0,0,0,1]) : nothing

    append!( [ ( [ ( CartesianIndex(i,j,k)== l ? -1 : 0 ) for i=1:H,j=1:H,k=1:4 ], [ 0, 0, 0, 0, 0 ])  for l in findall(potLinks .* s[1].links) ],
            [ ( [ ( CartesianIndex(i,j,k)== l ? 1 : 0 ) for i=1:H,j=1:H,k=1:4 ], [ 0, 0, 0, 0, 0 ]) for l in findall(potLinks .& (s[1].links .== 0)) ],
            [ ( zeros(H,H,4), pot) for pot in potMat ] )
end

maxsteps(env::MetamatEnv2) = env.maxsteps

function step!(env::MetamatEnv2, s, a)
    s = state(env)
    skeleton = copy(s[1])
    material = copy(s[2])
    act = findall(x->x!=0, a[1])
    if length(act) > 0
        skeleton = alteredSkeleton(s[1],act,a[1][act])
    else
        material = Material( ( [material.compressive, material.tensile, material.tensilePeak, material.peakStrain, material.crackStrainRatio] .+ (a[2] .* [10,1,0.01,0.1,0.05]) )... )
    end
    s′ = (skeleton,material)
    if ( s′[1].links[:,1,:] == falses(s′[1].size,4) || s′[1].links[:,s′[1].size-1,3:4] == falses(s′[1].size,2) )
        env.reward = 0.0
        env.reward, s′
    else
        model = modelFromSkeleton(s′[1], s′[2], 1, makeBasePoints(), makeBaseElements(s′[1].size), makeNodeWeights(s′[1].size,1), makeLinkWeights(s′[1].size,1))
        simulation = emptySimulation("current$(env.index)",2e-4)
        i=0
        while simulation.exit == 1 && i<3
            simulation = emptySimulation("current$(env.index)",2e-4)
            simulation.model = model
            runSimulation(simulation)
            run(`rm -r metamatcurrent$(env.index) metamatcurrent$(env.index)-restart.aux`)
            i+=1
        end
        env.reward = ( simulation.exit == 1 ? 0.0 : 0.1*energy(simulation.strain,simulation.stress)/(simulation.model.weight/310 + material.compressive/60 + material.tensile/6 + material.peakStrain/0.06) )
        env.reward, s′
    end
end

function finished(env::MetamatEnv2, s′)
    env.reward < 1e-10
end
