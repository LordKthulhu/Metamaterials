push!(LOAD_PATH,pwd()*"/Ressources")
using MetamatModule
using DelimitedFiles
using Plots
using Printf
using LaTeXStrings
using Reinforce
using Knet: Knet, dir, accuracy, progress, sgd, load, save, gc, Param, Data, minibatch, nll, relu, sigm, training, dropout, KnetArray, gpu
import Knet: param, param0, xavier
using StatsBase
using Base.Iterators: flatten
using IterTools: ncycle, takenth
using FileIO

param(d...; init=xavier, atype=atype())=Param(atype(init(d...)))
param0(d...; atype=atype())=param(d...; init=zeros, atype=atype)
xavier(o,i) = ( s = sqrt(2/(i+o)); 2s .* rand(o,i) .- s)
atype() = (gpu() >= 0 ? KnetArray{Float32} : Array{Float32})

const H=3
const unitSize = 1
const dEpsilon = 2e-4
const nodeWeights = makeNodeWeights(H,unitSize)
const linkWeights = makeLinkWeights(H,unitSize)
const basePoints = makeBasePoints()
const baseElements = makeBaseElements(H)

const indexes = [4,5,7,8,10,11,13,14,16,17,19,20,22,23,28,29,30,31,32,33]
const γ = 0.7
const threads = 40
const iterations = 5000

progressb = pBar(iterations,"Computing progess... ")


channel1 = Channel{Int64}(100)
channel2 = Channel(100)
channel3 = Channel(100)


ENV["GKSwstype"] = "100"

struct Chain
    layers; λ1; λ2
    Chain(layers...; λ1=0, λ2=0) = new(layers, λ1, λ2)
end

(c::Chain)(x) = (for l in c.layers; x = l(x); end; x)
(c::Chain)(d::Data) = mean(c(x,y) for (x,y) in d)
function (c::Chain)(x,y)
    loss = sum((c(x) .- hcat(y...)).^2)/length(y)
    if training() # Only apply regularization during training, only to weights, not biases.
        c.λ1 != 0 && (loss += c.λ1 * sum(sum(abs, l.w) for l in c.layers))
        c.λ2 != 0 && (loss += c.λ2 * sum(sum(abs2,l.w) for l in c.layers))
    end
    return loss
end

struct Layer; w; b; f; pdrop; end
Layer(i::Int,o::Int,f=relu; pdrop=0) = Layer(param(o,i),param0(o),f,pdrop)
(l::Layer)(x) = l.f.(l.w * dropout(x,l.pdrop) .+ l.b)

mutable struct MetamatPolicy <: AbstractPolicy end

function epsGreedy(s,A,Q;eps=0.95)
    if rand() < eps
        A[findmax([ Q( 0.5 .* reshape( s.links .+ a, (1,:))[indexes] ) for a in A ])[2]]
    else
        rand(A)
    end
end

function putIter(iterations::Int)
    for i in 1:iterations
        put!(channel1,i)
    end
end

function ep(env)
    while true
        iter = take!(channel1)
        shistory = []
        rhistory = []
        R = run_episode(env, MetamatPolicy()) do (s,a,r,s′)
           put!(channel2, (copy(s),a,r,copy(s′)))
           push!(shistory,copy(s′))
           push!(rhistory,r)
        end
        put!(channel3,(iter,shistory,rhistory))
        next!(progressb)
    end
end

function nnSGD(NN,previousNN,D)
    counter = 0
    while true
        s,a,r,s′ = take!(channel2)
        counter += 1
        push!(D,(s,a,r,s′))
        samples = sample(D,min(length(D),1000),replace=false)
        x = [ (sample[4].links) for sample in samples ]
        y = [ sample[3] + γ*findmax([ previousNN( 0.5 .* reshape( s′.links .+ a, (1,:))[indexes] ) for a in actions(env[1],s′) ])[1][1] for sample in samples ]
        dtrn = minibatch(Int.( hcat((e->(vcat(e...))).(x)...))[indexes,:], y, length(samples))
        learning = ( NN(dtrn) for x in takenth(sgd(NN,ncycle(dtrn,50)), length(dtrn)))
        learning = reshape(collect(Float32,flatten(learning)),(1,:))
        if counter % 10 == 0
            previousNN = deepcopy(NN)
        end
    end
end

function plotGeom(iterations)
    i = 1
    while i<=iterations
        iter, s, r = take!(channel3)
        run(`mkdir metamat$iter`)
        for i in 1:length(s)
            barPlots = plotGeometry(s[i])
            plt = plot(showaxis=false,size=(400,400),title="$(r[i])")
            for line in barPlots
                plot!(plt, line[1], line[2], lw =3, lc = :black, label = false)
            end
            png(plt,"metamat$iter/metamat"*lpad(iter,3,"0")*"-"*lpad(i,2,"0")*".png")
        end
        i += 1
        println("")
        @info "Simulation $iter, max improvement : $(maximum(r)-r[1])."
    end
end

env = [ MetamatEnv(H=H,maxsteps=10,index=i) for i in 1:threads ]

NN = Chain(Layer(20,20), Layer(20,10), Layer(10,1,identity); λ1=4f-5)

if isfile("NN.jld2")
    NN = FileIO.load("NN.jld2")["NN"]
else
    skeletons = load("Batch_09-06-2021_11:39:03/skeletons.jld")["skeletons"]
    results = readdlm("Batch_09-06-2021_11:39:03/results.csv",',')
    results = 10 .* results[4,:] ./ results[1,:]
    dtrn = minibatch( 0.5 .* Int.( hcat((e->(vcat(e...))).(skeletons)...))[indexes,:], results, 10; xtype = atype())
    r = ( NN(dtrn) for x in takenth(progress(sgd(NN,ncycle(dtrn,1000))), length(dtrn)))
    r = reshape(collect(Float32,flatten(r)),(1,:))
end

NN = Chain(NN.layers[1],Layer(20,20),Layer(20,20),NN.layers[2],NN.layers[3]; λ1=4f-5)

global previousNN = deepcopy(NN)

D = []

Reinforce.action(π::MetamatPolicy, r, s, A) = epsGreedy(s,actions(env[1],s),NN)

@async putIter(iterations)

for i in 1:threads
    @async ep(env[i])
end

@async try
    nnSGD(NN,previousNN,D)
catch err
    bt = catch_backtrace()
    println()
    showerror(stderr, err, bt)
end

t = @async plotGeom(iterations)

while !istaskdone(t)
    sleep(10)
end

isfile("NN.jld2") ? run(`rm NN.jld2`) : nothing
save("NN.jld2","NN",NN)
