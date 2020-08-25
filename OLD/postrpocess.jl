using  Statistics
using PyPlot
using DataFrames
using CSV

cd("C:/Users/u0125813/Desktop/Leuven/KU_Leuven/JULIA/Projects/CASE_STUDY_2/DATA/")

results = CSV.read("results.csv")
price_distr = CSV.read("price_distr.csv")

function t_interval(lenght, resolution)
    #DECLARATION OF TIME INTERVAL
    #e.g.
    #lenght = 24 h
    #resolution = 1 h or 0.25h(15mins)
    x = Int(lenght/ resolution)
    t_steps = [i for i in range(1, stop = x)]
end
t_steps = t_interval(24,1)
SW = results.SW
nms = results.model
#PLOTTING SWs
using PyCall, PyPlot
pygui(true)
sns = pyimport("seaborn")
sns.set()
sns.set_style("ticks")
sns.set_context("paper", font_scale=2.0)
fig = figure()
sns.barplot(nms, SW)
plot(nms, SW, marker="D",linestyle=":",linewidth =3.5, markersize = 8,color ="k",)
plt.ylabel("Total Social Welfare [EUR]")
grid()
sns.despine()


#PLOTTING THE PRFITS AS BARPLOTS
fig = figure()
bar_width = 0.35
indexes = 1:2:length(nms)*2
plt.bar(indexes, results.PROFIT_STR_GEN, label = "GEN", bar_width )
plt.bar(indexes .+ bar_width, results.PROFIT_ESS, label = "ESS", bar_width)
plt.bar(indexes .- bar_width, results.PROFIT_WIND_STR, label = "WIND", bar_width)
plt.xticks(indexes, nms)
grid()
legend()
# title("Profits of strategic agents")
plt.ylabel("PROFIT [EUR]")
# plt.xlabel("Models")


# long_results = stack(results,[:PROFIT_STR_GEN, :PROFIT_WIND_STR, :PROFIT_ESS])
# delete!(long_results, :SW)
# delete!(long_results, :Avg_MP)
# delete!(long_results, :Max_MP)
# delete!(long_results, :Min_MP)
# delete!(long_results, :Var_MP)
# sns.catplot("model", "value", "variable", data = long_results)
# plot(long_results)

# sns = pyimport("seaborn")
# sns.set()
# sns.set_style("whitegrid")
# sns.set_style("ticks")
# sns.set_context("paper")
# sns.despine()
fig = figure()
results.Avg_MP[7] - results.Avg_MP[11]

#Price distribution
sns.set_style("ticks")
sns.set_context("paper", font_scale=2.5)
fig = figure()
plot(t_steps, price_distr.lambda_MCP, label="MCP", linestyle="--",linewidth =4.5, marker="D",markersize = 8)
plot(t_steps, price_distr.lambda_MPEC, label="MPEC",linestyle="-.",linewidth =4.5,  marker="x",markersize = 12)
plot(t_steps, price_distr.lambda_EPEC_G4, label="EPEC",linestyle=":",linewidth =4.5,  marker="o",markersize = 8)
# plot(t_steps, price_distr.lambda_EPEC_peakstr, label="EPEC_peakload_GEN",linestyle=":",linewidth =4.5,  marker="o",markersize = 8,)
x_stk = 0:2:24
y_stk = 0:10:95
plt.xticks(x_stk)
plt.xlim((1,24))
plt.yticks(y_stk)
plt.ylabel("Market Price (EUR)")
plt.xlabel("Time (h)")
legend()
grid()



#Plotting the profits of the strategicegic agents
fig = figure()
plot(nms, results.PROFIT_STR_GEN, label="STR_GEN", linestyle=":",linewidth =3.5, marker="D",markersize = 8,)
plot(nms, results.PROFIT_ESS, label="STR_ESS", linewidth =3.5,  marker="*",markersize = 12)
plot(nms, results.PROFIT_WIND_STR, label="STR_WIND",linestyle=":",linewidth =3.5,  marker="x",markersize = 8,)
legend()
grid()
title("Profits of strategic agents")
plt.ylabel("Profit[EUR]")
plt.xlabel("Models")


res_dat = CSV.read("wind_PV.csv", header=true, delim= ";")
load_dat = CSV.read("load.csv",header=true, delim= ";")

#PRINTING THE INPUT PROTILES
fig = figure()
plot(t_steps, res_dat.Wind, label="Wind", linestyle=":",linewidth =3.5)
plot(t_steps, res_dat.PV, label="PV", linestyle="--",linewidth =3.5,  marker="o",markersize = 8)
legend(loc = "upper left")
grid()
x_stk = 0:2:24
plt.xticks(x_stk)
plt.xlim((1,24))
plt.ylabel("Renewable Generation (MW)")
plt.xlabel("Time (h)")

ax2 = plt.twinx()
ax2.plot(t_steps, load_dat.load, label="Load", color ="k",  linewidth =3.5)
ax2.legend(loc = "upper right")
ax2.set_ylabel("Load (MW)")


#MERIT ORDER CURVE

generators_dat = CSV.read("generators_data_original.csv", header=true, delim= ";")
MC = values(generators_dat[1,2:9])
cap = values(generators_dat[2,2:9])./1000
capacity = []
costs = []
for i in 1:length(cap)
    push!(capacity, cap[i])
    push!(costs, MC[i])
end
cap_cum  = cumsum(capacity)

sns = pyimport("seaborn")
sns.set()
sns.set_style("ticks")
sns.set_context("paper", font_scale=2.0)
plt.bar(0, costs[1], capacity[1], align ="edge",alpha=0.9, label = "GEN1")
plt.bar(cap_cum[1], costs[2], capacity[2], align ="edge",alpha=0.9, label = "GEN2")
plt.bar(cap_cum[2], costs[3], capacity[3], align ="edge",alpha=0.9, label = "GEN3")
plt.bar(cap_cum[3], costs[4], capacity[4], align ="edge",alpha=0.9, label = "GEN4")
plt.bar(cap_cum[4], costs[5], capacity[5], align ="edge",alpha=0.9, label = "GEN5")
plt.bar(cap_cum[5], costs[6], capacity[6], align ="edge",alpha=0.9, label = "GEN6")
plt.bar(cap_cum[6], costs[7], capacity[7], align ="edge",alpha=0.9, label = "GEN7")
plt.bar(cap_cum[7], costs[8], capacity[8], align ="edge",alpha=0.9, label = "GEN8")
plt.xticks(cap_cum)
plt.yticks(costs)
xlim((0,cap_cum[end]))
ylim((0,95))
grid()
legend()
# title("Profits of strategic agents")
plt.ylabel("Variable Cost (EUR/MWh)")
plt.xlabel("Capacity (GW)")
# sns.despine()
