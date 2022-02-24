using DelimitedFiles

function read_simplex_data(dataset::String)
    
    read(filename::String) = convert(Vector{Int64}, readdlm(filename, Int64)[:, 1])
    simplices = read("data/$(dataset)/$(dataset)-simplices.txt")
    nverts = read("data/$(dataset)/$(dataset)-nverts.txt")
    
    edges = Dict()
    Is = []
    Js = []
    w = []
    curr_ind = 1
    
    
    for nvert in nverts
        edge = Set(simplices[curr_ind:(curr_ind + nvert - 1)])
        curr_ind += nvert
        
        if (edge in keys(edges))
            edges[edge] = edges[edge]+1
        else
            edges[edge] = 1
        end
        
    end
    
    for (e,edge) in enumerate(keys(edges))
        append!(Is,edge)
        append!(Js,[e for _ in edge])
        append!(w,edges[edge])
    end
    
    m = length(w)
    
    return sparse(Is,Js,1), Float64.(w), edges
end

function read_email_data(dataset::String)
    
    edges = Dict()
    open("data_emails/$(dataset)/$(dataset).txt") do f
        for (e,line) in enumerate(eachline(f))
            nodes = [parse(Int64, v) for v in split(line, ' ')]
            edge = Set(nodes)
            if (edge in keys(edges))
                edges[edge] = edges[edge]+1
            else
                edges[edge] = 1
            end
        end
    end
    
    Is = []
    Js = []
    w = []
    for (e,edge) in enumerate(keys(edges))
        append!(Is,edge)
        append!(Js,[e for _ in edge])
        append!(w,edges[edge])
    end
    m = length(w)

    addresses = Dict()
    open("data_emails/$(dataset)/addresses-$(dataset).txt") do f
        for (e,line) in enumerate(eachline(f))
            a = split(line, ' ')
            node = parse(Int64, a[1])
            address = a[2]
            addresses[node] = address
        end
    end
    
    core_nodes = []
    read(filename::String) = convert(Vector{Int64}, readdlm(filename, Int64)[:, 1])
    core_nodes = read("data_emails/$(dataset)/core-$(dataset).txt")
    
    return sparse(Is,Js,1),  Float64.(w), edges, addresses, core_nodes
end

function read_citation_data(dataset::String)
    edges = []
    open("./data_citation/$dataset/hyperedges-$dataset.txt") do f
    for line in eachline(f)
        nodes = [parse(Int64, v) for v in split(line, [',', '\t'])]
        push!(edges, nodes)
    end
    end
    
    Is = []
    Js = []
    w = []
    for (e,edge) in enumerate(edges)
        append!(Is,edge)
        append!(Js,[e for _ in edge])
        append!(w,1)
    end
    m = length(w)
    
    return sparse(Is,Js,1),  Float64.(w), edges
end
