using SparseArrays, LinearAlgebra

using JLD2

include("utils/read_data_helpers.jl")
include("utils/cp_helpers.jl")
include("utils/UMHS_FT.jl")




function cp_analysis(dataset)
    
    t0 = @elapsed (B,weights,edges_dict,addresses,core_nodes) = read_email_data(dataset)
    println("$dataset: time reading = $t0")

    
    edges = collect(keys(edges_dict))
    num_edges = length(edges)
    num_nodes = size(B,1)
    max_edge_size = maximum(length.(edges))
    println("num nodes = $num_nodes, num hyperedges = $num_edges, largest hyperedge = $max_edge_size")

    # Hypergraph NSM ===========================================
    w = copy(weights)
    H = myHypergraph(Set.(edges), w, B, num_nodes, num_edges)
    t1 = @elapsed (x, x_array, x_er_array) = Hypergraph_NSM(H)
    println("Time cp hypergraph = $t1")
    t3 = @elapsed x_cp_profile = compute_cp_profile_3(x,H)
    println("Time cp profile x = $t3")
    # ==========================================================
    
    
    A =  B*(weights .* B')
    A = A - Diagonal(A)
    
    # Graph NSM ================================================
    t2 = @elapsed (y, y_array, y_er_array) = Graph_NSM(A)
    println("Time cp cliquegraph = $t2")
    t4 = @elapsed y_cp_profile = compute_cp_profile_3(y,H)
    println("Time cp profile y = $t4")
    # ==========================================================

    # BorgattiEverett ==========================================
    t6 = @elapsed (z,ranking) = BorgattiEverett(Int64.(A))
    println("Time Borgatti Everett = $t6")
    t7 = @elapsed z_cp_profile = compute_cp_profile_3(z,H)
    println("Time cp profile z = $t7")
    # ==========================================================
    
    # UMHS =====================================================
    nhits = 5
    t5 = @elapsed (s,ranking) = UMHS_FT(edges, num_nodes, nhits)
    println("Time UMHS $nhits nhits = $t5")
    t8 = @elapsed s_cp_profile = compute_cp_profile_3(s,H)
    println("Time cp profile s = $t8")
    # ==========================================================
            
    return Dict("hypergraph_info" => Dict("num_nodes" => num_nodes, "num_edges" => num_edges, "largest_edge" => max_edge_size), 
                "cp_hyper" => Dict("cp_score" => x, "cp_profile" => x_cp_profile, "time_compute_score" => t1), 
                "cp_clique" => Dict("cp_score" => y, "cp_profile" => y_cp_profile, "time_compute_score" => t2),
                "BE" => Dict("cp_score" => z, "cp_profile" => z_cp_profile, "time_compute_scpre" => t6),
                "umhs" => Dict("nhits"=>nhits, "cp_score" => s, "cp_profile" => s_cp_profile, "time_compute_scpre" => t5),
                "hypergraph" => Dict("incidence"=>B, "clique_adj" => A, "weights"=>weights))
    
end


function main(folder_path)
    datasets = readdir(folder_path)
    # datasets=["Enron", "ask-ubuntu", "math"] #"W3C", "Avocado",
    results = Dict()
    for (idx,dataset) in enumerate(datasets)

        results[dataset] = cp_analysis(dataset)        
        @save "./results/results_email_data_umhs_be_1-$idx.jld2" results
        println("  ")
    end
    return results
end


results = main("./data_emails/")
println("ALL DONE")
# @save "./results/results-$type-data-all.jld2" results

# figure(); 
# plot(results["Enron"]["hyper"][2],"r"); plot(results["Enron"]["clique"][2],"k")
# plot(resutls["W3C"]["hyper"][1],"r"); plot(resutls["W3C"]["clique"][1],"k")
