using Distributions,Combinatorics,PyPlot
include("utils/cp_helpers.jl")


function sigmoid(x; a = 1, s=10, t=0.5)
    return a ./ ( a .+ exp.(-s .* (x .- t) ) )
end

function edge_prob_1(N, edge, ϑ ; q = 10)
    μ = sum( ((N .- edge) ./ N).^q ) .^ (1/q)
    ξ = ϑ(length(edge))
    prob = (1/ϑ(2)) * ξ * sigmoid(μ)
    # prob = sigmoid((1/ϑ(2)) * ξ * μ)
    return prob, μ, ξ
end
function edge_prob_2(N, edge, ϑ ; q = 10)
    μ = sum( ((N .- edge) ./ N).^q ) .^ (1/q)
    ξ = ϑ(length(edge))
    # prob = (1/ϑ(2)) * ξ * sigmoid(μ)
    prob = sigmoid((1/ϑ(2)) * ξ * μ)
    return prob, μ, ξ
end
function edge_prob_3(N, edge, ϑ ; q = 10)
    μ = sum( ((N .- edge) ./ N).^q ) .^ (1/q)
    ξ = ϑ(length(edge))
    prob =  sigmoid(μ, a=ξ^2)
    # prob = sigmoid((1/ϑ(2)) * ξ * μ)
    return prob, μ, ξ
end

function incidence_matrix(edges)
    Is = []
    Js = []
    for (e,edge) in enumerate(edges)
        append!(Is,collect(edge))
        append!(Js,[e for _ in edge])
    end
    return sparse(Is,Js,1)
end


function sample_logistic_CP_hypergraph(N,edge_prob,ϑ)
    edges = []
    all_edges = []
    pp = []
    for edge in powerset(1:N, 2)
        pr, μ, ξ = edge_prob(N,edge,ϑ)
        # @show (edge , μ, ξ)
        push!(pp,  pr)
        push!(all_edges, edge)
        # println(edge)
        d = Binomial(1,pr)
        if rand(d) == 1
            push!(edges,edge)
        end
    end
    B = incidence_matrix(edges)

    return myHypergraph(Set.(edges),vec(ones(length(edges),1)),B,N,length(edges)), all_edges, pp
end


function compute_profiles_from_score(x,H)
    cp_profile = Dict()
    # t1 = @elapsed cp_profile["1"] = compute_cp_profile_1(x,H)
    # println("Time computing cp-profile 1: $t1")
    # t2 = @elapsed cp_profile["2"] = compute_cp_profile_2(x,H)
    # println("Time computing cp-profile 2: $t2")
    t3 = @elapsed cp_profile["3"] = compute_cp_profile_3(x,H)
    println("Time computing cp-profile 3: $t3")
    # t4 = @elapsed cp_profile["4"] = compute_cp_profile_4(x,H)
    # println("Time computing cp-profile 4: $t4")
    return cp_profile
end


function plot_edge_distribution(edges, Pr; type="plot")
    a = 0
    for l in unique(length.(edges))
        b = a+ sum(length.(edges).==l)
        if type == "semilogx"
            semilogx(a:b-1,Pr[length.(edges).==l])
        elseif type == "loglog"
            loglog(a:b-1,Pr[length.(edges).==l])
        elseif type == "semilogy"
            semilogy(a:b-1,Pr[length.(edges).==l])
        else
            plot(a:b-1,Pr[length.(edges).==l], "-o")
        end
        a = b
    end
end



N = 10
plot_type = "semilogy"

ϑ(x) = (1 ./ x) .^ 1
figure();
H,edges,Pr = sample_logistic_CP_hypergraph(N,edge_prob_2,ϑ)
plot_edge_distribution(edges, Pr; type=plot_type)
ylabel(L"P(e\in E)",fontsize=18)
yticks(fontsize=14)
xticks(fontsize=14)
# xticks(0:length(Pr)-1,edges,rotation=90)
# xticks(0:20,edges[1:21],rotation=90)
# xlabel(L"edges",fontsize=18)
# text(-2.5,-.1,L"e=",fontsize=18)


ϑ(x) = (1 ./ x) .^ 3
N = 15
H,edges,Pr = sample_logistic_CP_hypergraph(N,edge_prob_2,ϑ)
num_trials = 5
x_profile = []
y_profile = []
for trial in 1:num_trials 
    # N = 20 # number of nodes    

    B = H.incidence
    m = H.num_edges
    println("Number of edges = $m")
    x, x_array, x_er_array = Hypergraph_NSM(H,ϑ=ϑ)
    A =  B*(H.weights .* B')
    A = A - Diagonal(A)
    y, y_array, y_er_array = Graph_NSM(A)

    x_cp_profile = compute_profiles_from_score(x,H)#,ϑ=x->x./x)
    push!(x_profile, x_cp_profile["3"])
    y_cp_profile = compute_profiles_from_score(y,H)#,ϑ=x->x./x)
    push!(y_profile, y_cp_profile["3"])
end
xx = mean(x_profile)
yy = mean(y_profile)

figure()
l = length(xx)
xrange = range(0,1,length=l)

plot(xrange,xx,"r",linewidth=2)
plot(xrange,yy,"k")
ylabel(L"\frac{\# \{e : S\subseteq e\}}{\# \{e : S\cap e \neq \emptyset\}}", fontsize=18)
