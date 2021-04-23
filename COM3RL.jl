push!(LOAD_PATH,pwd()*"/Ressources")
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
    if rand() < 0.7
        A[findmax([ ( Q(s.links,a,theta)  ) for a in A ])[2]] #+ sqrt(2*log(t)/ (1+log(evalN(N,s.links,a))) )) for a in A ])[2]]
    else
        rand(A)
    end
end

function Q(state,a,theta)
    sum( theta[:,:,:,a[1][1],a[1][2],a[1][3],a[2]+1] .* state )
end

function evalN(N,state,a)
    sum( N[:,:,:,a[1][1],a[1][2],a[1][3],a[2]+1] .* state )
end

learning(oldTheta,theta) = sum( (oldTheta .- theta).^2 )#/(32*H^4)

env = [ MetamatEnv(H=H,maxsteps=5) for i in 1:Threads.nthreads() ]

#N = zeros(H,H,4,H,H,4,2)
#theta = randn(H,H,4,H,H,4,2)
#thetas = [deepcopy(theta)]

fullHistory = [ [] for i in 1:Threads.nthreads() ]
for iter in 1:100
    Reinforce.action(π::MetamatPolicy, r, s, A) = epsGreedy(s,actions(env[1],s),theta,N,iter)

    Threads.@threads for run in 1:Threads.nthreads()
        history = []
        delta = []
        R = run_episode(env[run], MetamatPolicy()) do (s,a,r,s′)
           push!(history, [deepcopy(s.links),a,r,deepcopy(s′.links)])
           #@info "On thread $(Threads.threadid()), reward is $r"
        end

        g = reverse(cumsum(reverse(getindex.(history[2:end],3).-getindex.(history[1:end-1],3))))
        i = 0
        for (s,a,r,s′) in history[2:end]
            i += 1
            features = findall(s)
            error = (g[i]-Q(s,a,theta))
            for j in features
                #N[j[1],j[2],j[3],a[1][1],a[1][2],a[1][3],a[2]+1] += 1
                theta[j[1],j[2],j[3],a[1][1],a[1][2],a[1][3],a[2]+1] += 0.5/length(features) * error #1/(1+log(N[j[1],j[2],j[3],a[1][1],a[1][2],a[1][3],a[2]+1])) * (g[i]-Q(s,a,theta))
            end
            push!(delta,error)
        end
        push!(fullHistory[Threads.threadid()],history)
        if length(delta)>0
            @info "On thread $(Threads.threadid()), reward improvement is $(history[end][3]-history[1][3]), errors on value function : $(sum(abs.(delta)))"
        end
    end


    @info "Theta was improved (squared error) by $(learning(thetas[end],theta))"
    push!(thetas,deepcopy(theta))
end
