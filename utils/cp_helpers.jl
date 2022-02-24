using SparseArrays,LinearAlgebra

struct myHypergraph
    edges::Array{Set{Int64}}
    weights::Vector{Float64}
    incidence::SparseArrays.SparseMatrixCSC{Float64,Int64} #Incidence matrix
    num_nodes::Int64
    num_edges::Int64
end

    
    
function μ(x,e,α)
    val = 0
    for i in e
        val += abs(x[i])^α
    end
    return val^(1/α)
end
    

function cp_loss(x, H::myHypergraph; ϑ = x -> 1 ./ x, α = 10)
    edges = H.edges
    weights = H.weights
    val = 0
    for (e,w) in zip(edges, weights)
        val += w * β(length(e)) * μ(x,e,α)
    end
    return val
end


function pm_iterator(x, H::myHypergraph, ϑ = x -> 1 ./ x, α = 10)
    B = H.incidence
    weights = H.weights
    y = abs.(x).^(α-2) .* x
    z = B * (weights .* ϑ(length.(H.edges)) .* (B'*abs.(x).^α).^(1/α-1)  )
    return y .* z
end


function  Hypergraph_NSM(H::myHypergraph; ϑ = x -> 1 ./ x, α = 10, p = 11, tol = 1e-8, x0 = "ones", verbose = false, maxiter = 200)
# Compute x with nonnegative entries that maximizes
# f(x) = sum_{e ∈ E} w(e)ϑ(e) μ_e(x)

    n  = H.num_nodes;

    if x0 == "ones"
        x0 = ones(n,1);
    end

    if verbose
        println("Nonlinear Power Method for Hypergraph CP:");
        println("-------------------------------");
        println("alpha:\t\t$α\np:\t\t\t$p\ntol:\t\t$tol");
    end

    pp = p/(p-1);
    x_array = x0./norm(x0,pp); er_array = [];
    
    for k = 1 : maxiter
        y = pm_iterator(vec(x0), H, ϑ, α)
        y = y./norm(y,pp);
        x = y.^(pp/p);
        x_array = [x_array x];
        er_array = [er_array; norm(x-x0)];
        if er_array[end] < tol || isnan(er_array[end])
            if verbose println("Num iter:\t$k"); end
            break;
        else
            x0 = copy(x);
        end
    end

    return x_array[:,end], x_array, er_array
end





function γ(S::Set{Int64},H::myHypergraph; ϑ = x -> 1 ./ x)

    num = 0
    den = 0
    for e in H.edges        
        if length(intersect(e,S)) > 0 #(S ∩ e) is not empty
            num += ϑ(length(e))
        end
        den += ϑ(length(e))
    end
    return num / den
    # numerator = 0
    # denominator = 0
    # for e in H.edges
    #     i = intersect(e,S)
    #     println("e=$e, S=$S, intersect = $i")
    #     if length(intersect(e,S)) == min(length(e),length(S)) # (S ∩ e) ∈ {S,e}
    #         numerator += 1/ϑ(length(e))
    #     end
    # 
    #     if length(intersect(e,S)) > 0 #(S ∩ e) is not empty
    #         denominator += 1/ϑ(length(e))
    #     end
    # end
    # return numerator/denominator    
end

function γ(S,v)
    g = 0
    for i in S
        g += v[i]
    end
    return g / sum(v)
end

function compute_cp_profile_1(x,H; ϑ = x -> 1 ./ x)
    # x_cp_profile = []; S = Set{Int64}();
    # B = H.incidence
    # 
    # v = B*vec(ϑ(length.(H.edges)))
    # 
    # for ind in sortperm(x)
    #     push!(S,ind)
    #     if length(S)>1
    #         push!(x_cp_profile,γ(S,H))
    #     end
    # end
    # return x_cp_profile
    x_cp_profile = []; S = Set{Int64}();
    edges = H.edges
    for ind in sortperm(x)
        push!(S,ind)
        if length(S)>1
            a = sum( (length.(intersect.(Ref(S),edges)).>0) .* ϑ(length.(edges))) #num of edges that have at least one node in S
            gamma = a
            push!(x_cp_profile,gamma)
        end
    end
    return x_cp_profile ./ sum(ϑ(length.(edges)))
end

function compute_cp_profile_2(x,H; ϑ = x -> 1 ./ x)
    x_cp_profile = []; S = Set{Int64}();
    edges = H.edges
    for ind in sortperm(x)
        push!(S,ind)
        if length(S)>1
            gamma = sum(length.(intersect.(Ref(S),edges)) ./ length.(edges))
            push!(x_cp_profile,gamma)
        end
    end
    return x_cp_profile
end

function compute_cp_profile_3(x,H; ϑ = x -> 1 ./ x)
    x_cp_profile = []; S = Set{Int64}();
    edges = H.edges
    for ind in sortperm(x)
        push!(S,ind)
        if length(S)>1
            a = sum(length.(intersect.(Ref(S),edges)).>0 ) #num of edges that have at least one node in S
            b = sum(issubset.(edges,Ref(S)) ) #num of edges that are completely contained in S
            gamma = b / a
            push!(x_cp_profile,gamma)
        end
    end
    return x_cp_profile
end

function compute_cp_profile_4(x,H; ϑ = x -> 1 ./ x)
    x_cp_profile = []; S = Set{Int64}();
    edges = H.edges
    for ind in sortperm(x)
        push!(S,ind)
        if length(S)>1
            a = sum( length.(intersect.(Ref(S),edges)) ./ length.(edges) ) #sum_e size(e ∩ S)/size(e)
            b = sum(length.(intersect.(Ref(S),edges)).>0 ) #num of edges that have at least one node in S
            gamma = a / b
            push!(x_cp_profile,gamma)
        end
    end
    return x_cp_profile
end



function F(A,x,a)
    i,j = findnz(A);
    m,n = size(A);
    aa = 1/a - 1;
    z = abs.(x).^a;
    Z = sparse(i,j, (z[i]+z[j]).^aa, m,n);
    nabla = 2 .* (abs.(x).^(a-1)) .* sign.(x) .* sum(A.*Z,dims=2);
    return nabla
end


function  Graph_NSM(A; a= 10, p=11, tol = 1e-8, x0 = "ones", verbose = false, maxiter = 200)
## Compute the x with nonnegative entries that maximizes
## f(x) = sum_{ij} A(i,j) mu_a(x(i),x(j))

    n  = size(A,1);

    if x0 == "ones"
        x0 = ones(n,1);
    end

    if verbose
        println("Nonlinear Power Method for Graph CP:");
        println("-------------------------------");
        println("alpha:\t\t$a\np:\t\t\t$p\ntol:\t\t$tol");
    end

    pp = p/(p-1);
    x_array = x0./norm(x0,pp); er_array = [];
    for k = 1 : maxiter
        y = F(A,x0,a);
        y = y./norm(y,pp);
        x = y.^(pp/p);
        x_array = [x_array x];
        er_array = [er_array; norm(x-x0)];
        if er_array[end] < tol || isnan(er_array[end])
            if verbose println("Num iter:\t$k"); end
            break;
        else
            x0 = copy(x);
        end
    end

    return x_array[:,end], x_array, er_array
end

function ranks(A)
# % RANKS : Replace numbers by their ranks.
# %
# % B = ranks(A)
# % returns a matrix the same shape as A, whose elements are
# % in the same relative order, but are integers from 1:k,
# % where k is the number of different values in A.
# %
    B = copy(A);
    p = sortperm(B[:]);
    D = ones(diff([-Inf; B[p]]));
    C = cumsum(D);
    B[p] = C;
    return B
end


function BorgattiEverett(A::SparseMatrixCSC, max_iter::Int64=10000,
                               tol::Float64=1e-8)
    n = size(A, 1)
    d = vec(sum(A, dims=2))
    c = rand(n)
    c[d .== 0] .= 0
    c /= norm(c, 2)
    for iter = 1:max_iter
        num = A * c
        denom = sum(c .^ 2) .- c .^ 2
        next_c = num ./ denom
        next_c /= norm(next_c, 2)
        diff = norm(next_c - c, 2)
        c = next_c
        if diff < tol; break; end
    end
    return c, sortperm(c, rev=true) #sort(collect(1:n), by= v -> c[v], rev=true)
end
