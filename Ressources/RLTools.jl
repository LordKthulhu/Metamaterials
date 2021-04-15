module MetamatEnv
using Reinforce

import Reinforce: reset!, actions, finished, step!, state

export
  MetamatEnv,
  reset!,
  step!,
  actions,
  finished

mutable struct MetamatEnv <: AbstractEnvironment
    state::Skeleton
    reward::Float64
    maxsteps::Int
end


MetamatEnv(;H=3, maxsteps = 10) = MetamatEnv(randomSkeleton(H), 0.0, maxsteps)

reset!(env::MetamatEnv) = (env.state = randomSkeleton(env.state.size); env.reward = 0.0; env)

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
    s.links[a[1]] = a[2]
    checkCompatibility!(s, a[1], a[2])

    material = Material(45.0,4.8)
    model = modelFromSkeleton(s, material, unitSize, nodeWeights, linkWeights)
    simulation = emptySimulation("current",dEpsilon)
    runSimulation(simulation)
    run(`rm -r metamatcurrent`)
    env.reward = energy(simulation.strain,simulation.stress)
    env.reward, s
end

function finished(env::MetamatEnv, s′)
    false
end

end
