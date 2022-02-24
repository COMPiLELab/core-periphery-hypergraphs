using Distributions,Combinatorics
include("cp_helpers.jl")


# N = 10 # number of nodes


function sigmoid(x; s=10, t=0.5)
    return 1 ./ ( 1 .+ exp.(-s .* (x .- t) ) )
end

function edge_prob(N, edge; ϑ = x -> 1 ./ x, q = 10)
    μ = sum( ((N .- edge) ./ N).^q ) .^ (1/q)
    ξ = ϑ(length(edge))
    # prob = 2 * ξ * sigmoid(μ)
    prob =  sigmoid(ξ * μ)

    return prob
end

# pp = []
# for edge in powerset(1:N, 2)
#     push!(pp,edge_prob(N,edge))
# end
# figure()
# plot(pp)

function incidence_matrix(edges)
    Is = []
    Js = []
    for (e,edge) in enumerate(edges)
        append!(Is,collect(edge))
        append!(Js,[e for _ in edge])
    end
    return sparse(Is,Js,1)
end


function sample_logistic_CP_hypergraph(N,edge_prob)
    edges = []
    for edge in powerset(1:N, 2)
        d = Binomial(1,edge_prob(N,edge))
        if rand(d) == 1
            push!(edges,edge)
        end
    end
    B = incidence_matrix(edges)

    return myHypergraph(Set.(edges),ones(length(edges),1),B,N,length(edges))
end
