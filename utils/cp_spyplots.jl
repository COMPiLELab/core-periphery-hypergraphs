using MATLAB
using JLD2

include("read_data_helpers.jl")
include("cp_helpers.jl")

function compute_thresholds(a)
    as = sort(a)
    step = Int64(round(length(a)/6))
    k = zeros(5,1)
    for i = 1 : 5
        k[i] = as[(i-1)*step + 1]
    end
    return vec(k) 
end

# @load "results/results_real_data_1-12.jld2"
@load "results/results-email-data-1-2.jld2"

datasets = ["congress-bills-full",      #time reading,n,m = 3.215227542, 1718, 105733
            "email-Eu-full",            #time reading,n,m = 1.009586095, 1005, 25148  
            "NDC-classes-full",         #time reading,n,m = 0.192053563, 1161, 1090
            "contact-high-school",      #time reading,n,m = 1.012831498, 327, 7818
            "contact-primary-school",   #time reading,n,m = 0.439013101, 242, 12704
            "NDC-substances-full",      #time reading,n,m = 0.564391331, 5556, 10273
            "email-Enron-full",         #time reading,n,m = 0.042804309, 148, 1514
            "tags-ask-ubuntu"]          #time reading,n,m = 4.08292788, 3029, 147222

datasets = ["phs-email-W3C", 
            "phs-email-Enron"]

dataset = datasets[1]
t0 = @elapsed (B,weights,edges_dict,addresses,core) = read_email_data(dataset) #read_simplex_data(dataset)
println("$dataset: time reading = $t0")

A =  B*(weights .* B')
A = A - Diagonal(A)

thresholds = compute_thresholds(unique(findnz(A)[3]))

k_mat = mxarray(thresholds)
A_mat = mxarray(A)
y_mat = mxarray(results[dataset]["cp_clique"]["cp_score"])
x_mat = mxarray(results[dataset]["cp_hyper"]["cp_score"])

x = results[dataset]["cp_hyper"]["cp_score"]
id = sortperm(x, rev=true)
[addresses[id[i]]  for i = 1:40]

mat"close all"
mat"subplot(121)"
mat"threshold_spy($A_mat,$y_mat,$k_mat)"
mat"subplot(122)"
mat"threshold_spy($A_mat,$x_mat,$k_mat)"
t = mxarray("title"); mat"sgtitle($t)"


function cp_intersection_profile(x,cp)
    profile = []
    xx = sortperm(x,rev=true)
    for i = 1 : 2*length(cp)
        push!(profile, length(intersect(xx[1:i],cp))/i)
    end
    return profile
end

figure()
x = results[dataset]["cp_hyper"]["cp_score"]
semilogx(cp_intersection_profile(x,core), "r", linewidth=2)
y = results[dataset]["cp_clique"]["cp_score"]
semilogx(cp_intersection_profile(y,core),"k")

# using ScikitLearn
# @sk_import metrics: average_precision_score
# 
# aorder = 1:length(x)
# resultx = sort(aorder,by=v->x[findall(aorder.==v)],rev=true)
# resultx01=zeros(length(x),1)
# for node in resultx
#      resultx01[findfirst(aorder.==node)]=1
# end
# precision_nsm=average_precision_score(core,resultx01)
# recovered_nsm=length(intersect(resultx,core))/length(x)

# 
