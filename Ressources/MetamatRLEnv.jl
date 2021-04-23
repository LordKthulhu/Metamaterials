
mutable struct MetamatEnv <: AbstractEnvironment
    state::Skeleton
    reward::Float64
    maxsteps::Int
end


MetamatEnv(;H=3, maxsteps = 10) = MetamatEnv(randomSkeleton(H), 1e-9, maxsteps)

reset!(env::MetamatEnv) = (env.state = randomSkeleton(env.state.size); env.reward = 1e-9; env)

function actions(env::MetamatEnv, s)
    H = s.size
    potLinks = trues(H,H,4)
    potLinks[:,1,1] = falses(H)
    potLinks[:,H,3:4] = falses(H,2)
    potLinks[H,:,:] = falses(H,4)
    for index in findall(s.links[:,1:end-1,3])
        potLinks[index[1],index[2]+1,1] = false
    end
    for index in findall(s.links[:,2:end,3])
        potLinks[index[1],index[2],3] = false
    end
    append!( [ [i,0] for i in findall(potLinks .& s.links) ], [ [i,1] for i in findall(potLinks .& (s.links .== 0) )] )
end

maxsteps(env::MetamatEnv) = env.maxsteps

function step!(env::MetamatEnv, s, a)
    s = state(env)
    s′ = alteredSkeleton(s,a[1],a[2])
    material = Material(true,45.0,4.8)
    model = modelFromSkeleton(s, material, 1, makeBasePoints(), makeBaseElements(s.size), makeNodeWeights(s.size,1), makeLinkWeights(s.size,1))
    simulation = emptySimulation("current$(Threads.threadid())",2e-4)
    simulation.model = model
    runSimulation(simulation)
    run(`rm -r metamatcurrent$(Threads.threadid()) metamatcurrent$(Threads.threadid())-restart.aux`)
    env.reward = simulation.exit == 1 ? 0.0 : 1000*energy(simulation.strain,simulation.stress)/simulation.model.weight
    env.reward, s′
end

function finished(env::MetamatEnv, s′)
    env.reward < 1e-10
end
