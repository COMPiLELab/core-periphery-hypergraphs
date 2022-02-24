using PyCall, PyPlot, Colors, SparseArrays

htx = pyimport("hypernetx")
nx = pyimport("networkx")
mpl = pyimport("matplotlib")

include("utils/cp_helpers.jl")
include("utils/UMHS_FT.jl")

function incidence_matrix(edges)
    Is = []
    Js = []
    for (e,edge) in enumerate(edges)
        append!(Is,collect(edge))
        append!(Js,[e for _ in edge])
    end
    return sparse(Is,Js,1)
end

function py_hypernetx_hypergraph(edges)
    edges_dict = Dict()
    for (id,e) in enumerate(edges) 
        edges_dict["$id"] = collect(e)
    end
    return htx.Hypergraph(edges_dict)
end

function cp_intersection_profile(x,cp)
    profile = []
    xx = sortperm(x,rev=true)
    for i = 1 : 2*length(cp)
        push!(profile, length(intersect(xx[1:i],cp))/i)
    end
    return profile
end

function hypercircle(edge_sizes)
    edges = []
    count = 1
    for (i,s) in enumerate(edge_sizes)
        push!(edges, collect(count:count+s-1))
        count = count+s-1
    end
    edges[end][end] = edges[1][1]
    
    core_nodes = Set([([[e[1],e[end]] for e in edges]...)...])

    u = 2*π/length(edge_sizes)
    pos = Dict()
    count = 1
    for (i,s) in enumerate(edge_sizes)
        pos[count] = 10 .*[cos(u*(i-1)), sin(u*(i-1))]
        count = count+s-1
    end
    for e in edges
        a = pos[e[1]]
        b = pos[e[end]]
        alpha = 1/(length(e[2:end-1])+1)
        for (i,node) in enumerate(e[2:end-1])
            pos[node] = i*alpha*(b-a) + a
        end
    end
    return Set.(edges), pos, collect(core_nodes)
    
end

function hyperplane()
    e1 = [1,2,3,4,5,6,7,8,9,10,11]
    e2 = [1,12]
    e3 = [1,13]
    e4 = [13,14]
    e5 = [13,15]
    edges = Set.([e1,e2,e3,e4,e5])
    pos = Dict()
    pos[1] = [0,0]
    pos[2] = [1,-.1].+[1,0]
    pos[3] = [2,-.2].+[1,0]
    pos[4] = [3,-.3].+[1,0]
    pos[5] = [4,-.4].+[1,0]
    pos[6] = [5,-.5].+[1,0]
    pos[7] = [-1,-.1].-[1,0]
    pos[8] = [-2,-.2].-[1,0]
    pos[9] = [-3,-.3].-[1,0]
    pos[10] = [-4,-.4].-[1,0]
    pos[11] = [-5,-.5].-[1,0]
    pos[12] = [0,2]
    pos[13] = [0,-5]
    pos[14] = [-1.6,-6.5]
    pos[15] = [1.6,-6.5]
    return edges, pos
end


## Hyperplane ## Uncomment this section to work with hyperplane
edges, my_pos = hyperplane()
n_nodes = maximum(maximum.(collect.(edges)))
n_edges = length(edges)
weights = ones(length(edges))
B = incidence_matrix(edges)
H = myHypergraph(edges,weights,B,n_nodes,n_edges)
H_py = py_hypernetx_hypergraph(edges)

## Hypercircle ## Uncomment this section to work with the hypercicle
# edge_sizes = [3,4,5,6,15]
# edges, my_pos, core_nodes = hypercircle(edge_sizes)
# n_nodes = maximum(maximum.(collect.(edges)))
# n_edges = length(edges)
# weights = ones(length(edges))
# B = incidence_matrix(edges)
# H = myHypergraph(edges,weights,B,n_nodes,n_edges)
# H_py = py_hypernetx_hypergraph(edges)


# Plot plain hypergraph
figure(); 
node_radius = 2
node_labels_kwargs=Dict("fontsize" => 12)
edges_kwargs=Dict("edgecolors" => "black")
myblue = [0, 0.604, 0.976]
nodes_kwargs=Dict("facecolors" => [  "#"*hex(RGB(myblue...)) for _ in H_py] )
htx.draw(H_py,pos=my_pos,label_alpha=0,node_radius=node_radius,with_node_labels=true,with_edge_labels=false, nodes_kwargs=nodes_kwargs, edges_kwargs=edges_kwargs, node_labels_kwargs=node_labels_kwargs)
axis("on");xticks([]);yticks([])

## Compute core scores and cp_profiles
x, x_array, x_er_array = Hypergraph_NSM(H)
A =  B*(weights .* B')
A = A - Diagonal(A)
y, y_array, y_er_array = Graph_NSM(A)
z,sort_z = BorgattiEverett(A); z = z ./ norm(z,Inf)
x_cp_profile = compute_cp_profile_3(x,H)#,ϑ=x->x./x)
y_cp_profile = compute_cp_profile_3(y,H)#,ϑ=x->x./x)
z_cp_profile = compute_cp_profile_3(z,H)





ftlegend = 13
ftticks = 13
fttitle = 17
ftlabels = 15


## Three panels figure: hypergraph cp, clique cp, cp profile
figure(); 
node_radius = 2
node_labels_kwargs=Dict("fontsize" => 12)
edges_kwargs=Dict("edgecolors" => "black")
myblue = [0, 0.604, 0.976]
myorange = [1.000, 0.498, 0.055]
# node_color_string(score) = "#"*hex(RGB((myblue.*score + [.8,.8,.8].*(1-score))...))
cmap = mpl.cm.get_cmap("Blues")
im = imshow([x x] .+ 0.1,cmap=cmap)
clim(0, 1)
plt.clf()
node_color_string(score) = cmap(score+.05) #"#"*hex(RGBA(collect(cmap(score))...))
with_node_labels = false

subplot(131)
nodes_kwargs=Dict("facecolors" => [  node_color_string(x[v]) for v in H_py] )
htx.draw(H_py,pos=my_pos,node_radius=node_radius,with_node_labels=with_node_labels,with_edge_labels=false, nodes_kwargs=nodes_kwargs, edges_kwargs=edges_kwargs, node_labels_kwargs=node_labels_kwargs)
axis("on");xticks([]);yticks([])
colorbar(im)
title("HyperNSM core score", fontsize=fttitle)

subplot(132)
nodes_kwargs=Dict("facecolors" => [  node_color_string(y[v]) for v in H_py] )
htx.draw(H_py,pos=my_pos,node_radius=node_radius,with_node_labels=with_node_labels,with_edge_labels=false, nodes_kwargs=nodes_kwargs, edges_kwargs=edges_kwargs, node_labels_kwargs=node_labels_kwargs)
axis("on");xticks([]);yticks([])
colorbar(im)
title("GraphNSM core score", fontsize=fttitle)

subplot(133)
plot(x_cp_profile,"-o",color="#"*hex(RGB(myblue...)), linewidth=2, label="HyperNSM")
plot(y_cp_profile, "k-x", linewidth=1, label="GraphNSM")
xticks(fontsize=ftticks)
yticks(fontsize=ftticks)
legend(fontsize=ftlegend)
xticks(3:5:length(x_cp_profile)-1 , 5:5:length(x_cp_profile)+1)
title("CP Profile", fontsize=fttitle)


savefig("hypercycle.pdf", format="pdf",  dpi=600)

# subplot(133)
# nodes_kwargs=Dict("facecolors" => [  node_color_string(z[v]) for v in H_py] )
# htx.draw(H_py,pos=my_pos,node_radius=node_radius,with_node_labels=with_node_labels,with_edge_labels=false, nodes_kwargs=nodes_kwargs, edges_kwargs=edges_kwargs, node_labels_kwargs=node_labels_kwargs)
# axis("on");xticks([]);yticks([])
# colorbar(im)
# title("Borgatti-Everett core score", fontsize=fttitle)

1; 
