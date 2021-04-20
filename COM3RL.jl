push!(LOAD_PATH,"/mnt/c/Users/Pierre/Metamaterials/Ressources")
using MetamatModule
using DelimitedFiles
using Plots
using Printf
using LaTeXStrings
using Reinforce

const H,iterations,randomMat,parameters = parseArguments()
const unitSize = 1
const dEpsilon = 2e-4
const nodeWeights = makeNodeWeights(H,unitSize)
const linkWeights = makeLinkWeights(H,unitSize)
const basePoints = makeBasePoints()
const baseElements = makeBaseElements(H)

ENV["GKSwstype"] = "100"

mutable struct MetamatPolicy <: AbstractPolicy end

function epsGreedy(s,A,theta,N,t)
    A[findmax([ ( Q(s.links,a,theta) + sqrt(20*log(t)/ evalN(N,s.links,a)) ) for a in A ])[2]]
end

function Q(state,a,theta)
    sum( theta[:,:,:,a[1][1],a[1][2],a[1][3],a[2]+1] .* state )
end

function evalN(N,state,a)
    sum( N[:,:,:,a[1][1],a[1][2],a[1][3],a[2]+1] .* state )
end

learning(oldTheta,theta) = sum( (oldTheta .- theta).^2 )#/(32*H^4)

env = [ MetamatEnv(H=H) for i in 1:Threads.nthreads() ]

N = zeros(H,H,4,H,H,4,2)
#theta = randn(H,H,4,H,H,4,2)
thetas = [deepcopy(theta)]

fullHistory = [ [] for i in 1:Threads.nthreads() ]
for iter in 10:20
    Reinforce.action(π::MetamatPolicy, r, s, A) = epsGreedy(s,actions(env[1],s),theta,N,iter)

    Threads.@threads for run in 1:Threads.nthreads()
        history = []
        R = run_episode(env[run], MetamatPolicy()) do (s,a,r,s′)
           push!(history, [deepcopy(s.links),a,r,deepcopy(s′.links)])
           #@info "On thread $(Threads.threadid()), reward is $r"
        end

        g = reverse(cumsum(reverse(getindex.(history[2:end],3).-getindex.(history[1:end-1],3))))
        i = 0
        for (s,a,r,s′) in history[2:end]
            i += 1
            for j in findall(s)
                N[j[1],j[2],j[3],a[1][1],a[1][2],a[1][3],a[2]+1] += 1
                theta[j[1],j[2],j[3],a[1][1],a[1][2],a[1][3],a[2]+1] += 1/N[j[1],j[2],j[3],a[1][1],a[1][2],a[1][3],a[2]+1] * (g[i]-Q(s,a,theta))
            end
        end
        push!(fullHistory[Threads.threadid()],history)
        @info "On thread $(Threads.threadid()), average reward is $(mean(getindex.(history,3)))"
    end

    @info "Theta was improved (squared error) by $(learning(thetas[end],theta))"
    push!(thetas,deepcopy(theta))
end
