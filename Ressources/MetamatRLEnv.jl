
mutable struct MetamatEnv <: AbstractEnvironment
    state::Skeleton
    reward::Float64
    maxsteps::Int
    index::Int
end


MetamatEnv(;H=3, maxsteps = 20, index = 1) = MetamatEnv(randomSkeleton(H), 1e-9, maxsteps, index)

reset!(env::MetamatEnv) = (env.state = randomSkeleton(env.state.size); env.reward = 1e-9; env)

function actions(env::MetamatEnv, s)
    H = s.size
    potLinks = trues(H,H,4)
    potLinks[:,1,1] = falses(H)
    potLinks[:,H,3:4] = falses(H,2)
    potLinks[H,:,1:3] = falses(H,3)
    for index in findall(s.links[1:end-1,1:end-1,3])
        potLinks[index[1],index[2]+1,1] = false
    end
    for index in findall(s.links[1:end-1,2:end,1])
        potLinks[index[1],index[2],3] = false
    end
    append!( [ [ ( CartesianIndex(i,j,k)== l ? -1 : 0 ) for i=1:3,j=1:3,k=1:4 ] for l in findall(potLinks .* s.links) ],
            [ [ ( CartesianIndex(i,j,k)== l ? 1 : 0 ) for i=1:3,j=1:3,k=1:4 ] for l in findall(potLinks .& (s.links .== 0)) ])
    #append!( [ [i,0] for i in findall(potLinks .& s.links) ], [ [i,1] for i in findall(potLinks .& (s.links .== 0) )] )
end

maxsteps(env::MetamatEnv) = env.maxsteps

function step!(env::MetamatEnv, s, a)
    s = state(env)
    act = findall(x->x!=0, a)
    s′ = alteredSkeleton(s,act,a[act])
    if ( s′.links[:,1,:] == falses(s′.size,4) || s′.links[:,s′.size-1,3:4] == falses(s′.size,2) )
        env.reward = 0.0
        env.reward, s′
    else
        material = Material(45.0,4.8)
        model = modelFromSkeleton(s′, material, 1, makeBasePoints(), makeBaseElements(s′.size), makeNodeWeights(s′.size,1), makeLinkWeights(s′.size,1))
        simulation = emptySimulation("current$(env.index)",2e-4)
        i=0
        while simulation.exit == 1 && i<3
            simulation = emptySimulation("current$(env.index)",2e-4)
            simulation.model = model
            runSimulation(simulation)
            run(`rm -r metamatcurrent$(env.index) metamatcurrent$(env.index)-restart.aux`)
            i+=1
        end
        env.reward = ( simulation.exit == 1 ? 0.0 : 10*energy(simulation.strain,simulation.stress)/simulation.model.weight )
        env.reward, s′
    end
end

function finished(env::MetamatEnv, s′)
    env.reward < 1e-10
end
