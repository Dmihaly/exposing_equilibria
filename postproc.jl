using Statistics
using Query
using CSV
using DataFrames
using LinearAlgebra
using PyCall, PyPlot
# using Plots
pygui(true)


simu_EPEC_quantity_empty = CSV.read("C:/Users/u0125813/Desktop/Leuven/KU_Leuven/EPECs/scr_nlp/RESULTS/CASE2/simu_EPEC_quantity_empty.csv"; copycols = true)
simu_EPEC_quantity_competitive = CSV.read("C:/Users/u0125813/Desktop/Leuven/KU_Leuven/EPECs/scr_nlp/RESULTS/CASE2/simu_EPEC_quantity_competitive.csv"; copycols = true)
simu_EPEC_quantity_collusive = CSV.read("C:/Users/u0125813/Desktop/Leuven/KU_Leuven/EPECs/scr_nlp/RESULTS/CASE2/simu_EPEC_quantity_collusive.csv"; copycols = true)
simu_EPEC_quantity_ESS_fav = CSV.read("C:/Users/u0125813/Desktop/Leuven/KU_Leuven/EPECs/scr_nlp/RESULTS/CASE2/simu_EPEC_quantity_ESS_fav.csv"; copycols = true)
simu_EPEC_quantity_PC = CSV.read("C:/Users/u0125813/Desktop/Leuven/KU_Leuven/EPECs/scr_nlp/RESULTS/CASE2/simu_EPEC_quantity_PC.csv"; copycols = true)


simu_EPEC_quantity_empty_confirmed = simu_EPEC_quantity_empty |> @filter(_.confirmedNE == true) |> DataFrame
simu_EPEC_quantity_competitive_confirmed = simu_EPEC_quantity_competitive |> @filter(_.confirmedNE == true) |> DataFrame
simu_EPEC_quantity_collusive_confirmed = simu_EPEC_quantity_collusive |> @filter(_.confirmedNE == true) |> DataFrame
simu_EPEC_quantity_ESS_fav_confirmed = simu_EPEC_quantity_ESS_fav |> @filter(_.confirmedNE == true) |> DataFrame

######RESULTS OF THE PRICE CONSISTENT MODEL#########
avg_SW_PC = Statistics.mean(simu_EPEC_quantity_PC.SW[:])
avg_PS_PC = Statistics.mean(simu_EPEC_quantity_PC.PS[:])
avg_profit_ESS_PC = Statistics.mean(simu_EPEC_quantity_PC.ESS[:])
avg_profit_GEN_PC = Statistics.mean(simu_EPEC_quantity_PC.GEN[:])
avg_profit_WIND_PC = Statistics.mean(simu_EPEC_quantity_PC.WIND_STR[:])



######RESULTS OF THE MODEL WITH EMPTY OBJECTIVE#########
avg_SW_empty = Statistics.mean(simu_EPEC_quantity_empty.SW[:])
avg_PS_empty = Statistics.mean(simu_EPEC_quantity_empty.PS[:])
avg_profit_ESS_empty = Statistics.mean(simu_EPEC_quantity_empty.ESS[:])
avg_profit_GEN_empty = Statistics.mean(simu_EPEC_quantity_empty.GEN[:])
avg_profit_WIND_empty = Statistics.mean(simu_EPEC_quantity_empty.WIND_STR[:])

norm_profit_ESS_empty = simu_EPEC_quantity_empty.ESS[:] / avg_profit_ESS_empty
norm_profit_GEN_empty = simu_EPEC_quantity_empty.GEN[:] / avg_profit_GEN_empty
norm_profit_WIND_empty = simu_EPEC_quantity_empty.WIND_STR[:] / avg_profit_WIND_empty

norms_empty = vcat(norm_profit_ESS_empty, norm_profit_GEN_empty, norm_profit_WIND_empty)


avg_SW_empty_confirmed = Statistics.mean(simu_EPEC_quantity_empty_confirmed.SW[:])
avg_PS_empty_confirmed = Statistics.mean(simu_EPEC_quantity_empty_confirmed.PS[:])
avg_profit_ESS_empty_confirmed = Statistics.mean(simu_EPEC_quantity_empty_confirmed.ESS[:])
avg_profit_GEN_empty_confirmed = Statistics.mean(simu_EPEC_quantity_empty_confirmed.GEN[:])
avg_profit_WIND_empty_confirmed = Statistics.mean(simu_EPEC_quantity_empty_confirmed.WIND_STR[:])

norm_profit_ESS_empty_confirmed = simu_EPEC_quantity_empty_confirmed.ESS[:] / avg_profit_ESS_empty_confirmed
norm_profit_GEN_empty_confirmed = simu_EPEC_quantity_empty_confirmed.GEN[:] / avg_profit_GEN_empty_confirmed
norm_profit_WIND_empty_confirmed = simu_EPEC_quantity_empty_confirmed.WIND_STR[:] / avg_profit_WIND_empty_confirmed

norms_empty_confirmed = vcat(norm_profit_ESS_empty_confirmed, norm_profit_GEN_empty_confirmed, norm_profit_WIND_empty_confirmed)

######RESULTS OF THE MODEL WITH COMPETITIVE OBJECTIVE#########
avg_SW_competitive = Statistics.mean(simu_EPEC_quantity_competitive.SW[:])
avg_PS_competitive = Statistics.mean(simu_EPEC_quantity_competitive.PS[:])
avg_profit_ESS_competitive = Statistics.mean(simu_EPEC_quantity_competitive.ESS[:])
avg_profit_GEN_competitive = Statistics.mean(simu_EPEC_quantity_competitive.GEN[:])
avg_profit_WIND_competitive = Statistics.mean(simu_EPEC_quantity_competitive.WIND_STR[:])

norm_profit_ESS_competitive = simu_EPEC_quantity_competitive.ESS[:] / avg_profit_ESS_empty
norm_profit_GEN_competitive = simu_EPEC_quantity_competitive.GEN[:] / avg_profit_GEN_empty
norm_profit_WIND_competitive = simu_EPEC_quantity_competitive.WIND_STR[:] / avg_profit_WIND_empty

norms_competitive = vcat(norm_profit_ESS_competitive, norm_profit_GEN_competitive, norm_profit_WIND_competitive)

avg_SW_competitive_confirmed = Statistics.mean(simu_EPEC_quantity_competitive_confirmed.SW[:])
avg_PS_competitive_confirmed = Statistics.mean(simu_EPEC_quantity_competitive_confirmed.PS[:])
avg_profit_ESS_competitive_confirmed = Statistics.mean(simu_EPEC_quantity_competitive_confirmed.ESS[:])
avg_profit_GEN_competitive_confirmed = Statistics.mean(simu_EPEC_quantity_competitive_confirmed.GEN[:])
avg_profit_WIND_competitive_confirmed = Statistics.mean(simu_EPEC_quantity_competitive_confirmed.WIND_STR[:])

norm_profit_ESS_competitive_confirmed = simu_EPEC_quantity_competitive_confirmed.ESS[:] / avg_profit_ESS_empty_confirmed
norm_profit_GEN_competitive_confirmed = simu_EPEC_quantity_competitive_confirmed.GEN[:] / avg_profit_GEN_empty_confirmed
norm_profit_WIND_competitive_confirmed = simu_EPEC_quantity_competitive_confirmed.WIND_STR[:] / avg_profit_WIND_empty_confirmed

norms_competitive_confirmed = vcat(norm_profit_ESS_competitive_confirmed, norm_profit_GEN_competitive_confirmed, norm_profit_WIND_competitive_confirmed)

######RESULTS OF THE MODEL WITH COLLUSIVE OBJECTIVE#########
avg_SW_collusive = Statistics.mean(simu_EPEC_quantity_collusive.SW[:])
avg_PS_collusive = Statistics.mean(simu_EPEC_quantity_collusive.PS[:])
avg_profit_ESS_collusive = Statistics.mean(simu_EPEC_quantity_collusive.ESS[:])
avg_profit_GEN_collusive = Statistics.mean(simu_EPEC_quantity_collusive.GEN[:])
avg_profit_WIND_collusive = Statistics.mean(simu_EPEC_quantity_collusive.WIND_STR[:])

norm_profit_ESS_collusive = simu_EPEC_quantity_collusive.ESS[:] / avg_profit_ESS_empty
norm_profit_GEN_collusive = simu_EPEC_quantity_collusive.GEN[:] / avg_profit_GEN_empty
norm_profit_WIND_collusive = simu_EPEC_quantity_collusive.WIND_STR[:] / avg_profit_WIND_empty

norms_collusive = vcat(norm_profit_ESS_collusive, norm_profit_GEN_collusive, norm_profit_WIND_collusive)


avg_SW_collusive_confirmed = Statistics.mean(simu_EPEC_quantity_collusive_confirmed.SW[:])
avg_PS_collusive_confirmed = Statistics.mean(simu_EPEC_quantity_collusive_confirmed.PS[:])
avg_profit_ESS_collusive_confirmed = Statistics.mean(simu_EPEC_quantity_collusive_confirmed.ESS[:])
avg_profit_GEN_collusive_confirmed = Statistics.mean(simu_EPEC_quantity_collusive_confirmed.GEN[:])
avg_profit_WIND_collusive_confirmed = Statistics.mean(simu_EPEC_quantity_collusive_confirmed.WIND_STR[:])

norm_profit_ESS_collusive_confirmed = simu_EPEC_quantity_collusive_confirmed.ESS[:] / avg_profit_ESS_empty_confirmed
norm_profit_GEN_collusive_confirmed = simu_EPEC_quantity_collusive_confirmed.GEN[:] / avg_profit_GEN_empty_confirmed
norm_profit_WIND_collusive_confirmed = simu_EPEC_quantity_collusive_confirmed.WIND_STR[:] / avg_profit_WIND_empty_confirmed

norms_collusive_confirmed = vcat(norm_profit_ESS_collusive_confirmed, norm_profit_GEN_collusive_confirmed, norm_profit_WIND_collusive_confirmed)


######RESULTS OF THE MODEL WITH ESS FAVORED OBJECTIVE#########
avg_SW_ESS_fav = Statistics.mean(simu_EPEC_quantity_ESS_fav.SW[:])
avg_PS_ESS_fav = Statistics.mean(simu_EPEC_quantity_ESS_fav.PS[:])
avg_profit_ESS_ESS_fav = Statistics.mean(simu_EPEC_quantity_ESS_fav.ESS[:])
avg_profit_GEN_ESS_fav = Statistics.mean(simu_EPEC_quantity_ESS_fav.GEN[:])
avg_profit_WIND_ESS_fav = Statistics.mean(simu_EPEC_quantity_ESS_fav.WIND_STR[:])

norm_profit_ESS_ESS_fav = simu_EPEC_quantity_ESS_fav.ESS[:] / avg_profit_ESS_empty
norm_profit_GEN_ESS_fav = simu_EPEC_quantity_ESS_fav.GEN[:] / avg_profit_GEN_empty
norm_profit_WIND_ESS_fav = simu_EPEC_quantity_ESS_fav.WIND_STR[:] / avg_profit_WIND_empty

norms_ESS_fav = vcat(norm_profit_ESS_ESS_fav, norm_profit_GEN_ESS_fav, norm_profit_WIND_ESS_fav)


avg_SW_ESS_fav_confirmed = Statistics.mean(simu_EPEC_quantity_ESS_fav_confirmed.SW[:])
avg_PS_ESS_fav_confirmed = Statistics.mean(simu_EPEC_quantity_ESS_fav_confirmed.PS[:])
avg_profit_ESS_ESS_fav_confirmed = Statistics.mean(simu_EPEC_quantity_ESS_fav_confirmed.ESS[:])
avg_profit_GEN_ESS_fav_confirmed = Statistics.mean(simu_EPEC_quantity_ESS_fav_confirmed.GEN[:])
avg_profit_WIND_ESS_fav_confirmed = Statistics.mean(simu_EPEC_quantity_ESS_fav_confirmed.WIND_STR[:])

norm_profit_ESS_ESS_fav_confirmed = simu_EPEC_quantity_ESS_fav_confirmed.ESS[:] / avg_profit_ESS_empty_confirmed
norm_profit_GEN_ESS_fav_confirmed = simu_EPEC_quantity_ESS_fav_confirmed.GEN[:] / avg_profit_GEN_empty_confirmed
norm_profit_WIND_ESS_fav_confirmed = simu_EPEC_quantity_ESS_fav_confirmed.WIND_STR[:] / avg_profit_WIND_empty_confirmed

norms_ESS_fav_confirmed = vcat(norm_profit_ESS_ESS_fav_confirmed, norm_profit_GEN_ESS_fav_confirmed, norm_profit_WIND_ESS_fav_confirmed)



#PLOTTING SWs
SW = [SW_MCP, avg_SW_empty, avg_SW_competitive, avg_SW_collusive, avg_SW_ESS_fav]
SW_confirmed = [SW_MCP, avg_SW_empty_confirmed, avg_SW_competitive_confirmed, avg_SW_collusive_confirmed, avg_SW_ESS_fav_confirmed]
PS = [PS_MCP, avg_PS_empty, avg_PS_competitive, avg_PS_collusive, avg_PS_ESS_fav]
PS_confirmed = [PS_MCP, avg_PS_empty_confirmed, avg_PS_competitive_confirmed, avg_PS_collusive_confirmed, avg_PS_ESS_fav_confirmed]
nms = ["empty", "competitive", "collusive", "ESS_fav"]

function calc_percentage(arr::Array)
    SW_perc = []
    for i in 2:length(arr)
        # if i  == 1
        #     push!(SW_perc, 0)
        # else
            push!(SW_perc, ((arr[i] -  arr[1])/arr[1])*100)
        # end
    end
    return SW_perc
end

SW_perc = calc_percentage(SW)
SW_perc_confirmed = calc_percentage(SW_confirmed)

PS_perc = calc_percentage(PS)
PS_perc_confirmed = calc_percentage(PS_confirmed)


cd("C:/Users/u0125813/Desktop/Leuven/KU_Leuven/EPECs/scr_nlp/RESULTS/")
###FOR SW######
sns = pyimport("seaborn")
sns.set()
sns.set_style("ticks")
sns.set_context("paper", font_scale=1.5)
bar_width = 0.5
fig = figure()
indexes = 1:2:length(nms)*2
plt.bar(indexes, SW_perc, label = "all", bar_width, align="center")
plt.bar(indexes .+ bar_width, SW_perc_confirmed, label = "confirmed", bar_width, align="center")
plt.xticks(indexes, nms)
plt.grid()
plt.legend()
plt.ylabel("Total social welfare difference compared to perfect competition [%]")


###FOR PROD SURPLUS######
sns = pyimport("seaborn")
sns.set()
sns.set_style("ticks")
sns.set_context("paper", font_scale=1.5)
bar_width = 0.5
fig = figure()
indexes = 1:2:length(nms)*2
plt.bar(indexes, PS_perc, label = "all", bar_width)
plt.bar(indexes .+ bar_width, PS_perc_confirmed, label = "confirmed", bar_width)
plt.xticks(indexes, nms)
plt.grid()
legend()
plt.ylabel("Total producer surplus difference compared to perfect competition [%]")


sns = pyimport("seaborn")
sns.set()
sns.set_style("ticks")
sns.set_context("paper", font_scale=1.3)
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.scatter(simu_EPEC_quantity_empty.ESS, simu_EPEC_quantity_empty.GEN, simu_EPEC_quantity_empty.WIND_STR,s=40, label = "empty", marker ="+")
ax.scatter(simu_EPEC_quantity_competitive.ESS, simu_EPEC_quantity_competitive.GEN, simu_EPEC_quantity_competitive.WIND_STR,s=40, label = "competitive")
ax.scatter(simu_EPEC_quantity_collusive.ESS, simu_EPEC_quantity_collusive.GEN, simu_EPEC_quantity_collusive.WIND_STR,s=40, label = "collusive", marker = "x")
ax.scatter(simu_EPEC_quantity_ESS_fav.ESS, simu_EPEC_quantity_ESS_fav.GEN, simu_EPEC_quantity_ESS_fav.WIND_STR,s=40, label = "ESS_fav", marker = "p")
ax.set_xlabel("Profit of ESS")
ax.set_ylabel("Profit of CG")
ax.set_zlabel("Profit of RG")
ax.view_init(30, 185)
ax.legend()
plt.show()
plt.savefig("3Dscatter_profits",bbox_inches="tight",dpi=600)


sns = pyimport("seaborn")
sns.set()
sns.set_style("ticks")
sns.set_context("paper", font_scale=1.3)
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.scatter(simu_EPEC_quantity_empty_confirmed.ESS, simu_EPEC_quantity_empty_confirmed.GEN, simu_EPEC_quantity_empty_confirmed.WIND_STR,s=40, label = "empty", marker ="+")
ax.scatter(simu_EPEC_quantity_competitive_confirmed.ESS, simu_EPEC_quantity_competitive_confirmed.GEN, simu_EPEC_quantity_competitive_confirmed.WIND_STR,s=40, label = "competitive")
ax.scatter(simu_EPEC_quantity_collusive_confirmed.ESS, simu_EPEC_quantity_collusive_confirmed.GEN, simu_EPEC_quantity_collusive_confirmed.WIND_STR,s=40, label = "collusive", marker = "x")
ax.scatter(simu_EPEC_quantity_ESS_fav_confirmed.ESS, simu_EPEC_quantity_ESS_fav_confirmed.GEN, simu_EPEC_quantity_ESS_fav_confirmed.WIND_STR,s=40, label = "ESS_fav", marker = "p")
ax.set_xlabel("Profit of ESS")
ax.set_ylabel("Profit of CG")
ax.set_zlabel("Profit of RG")
ax.view_init(30, 185)
ax.legend()
plt.show()
plt.savefig("3Dscatter_profits_conf",bbox_inches="tight",dpi=600)


#Social Welfare vs Producer Surplus
sns = pyimport("seaborn")
sns.set()
sns.set_style("ticks")
sns.set_context("paper", font_scale=1.0)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(simu_EPEC_quantity_empty.SW, simu_EPEC_quantity_empty.PS, s=60, label = "empty")
ax.scatter(simu_EPEC_quantity_competitive.SW, simu_EPEC_quantity_competitive.PS, s=60, label = "competitive")
ax.scatter(simu_EPEC_quantity_collusive.SW, simu_EPEC_quantity_collusive.PS, s=60, label = "collusive")
ax.scatter(simu_EPEC_quantity_ESS_fav.SW, simu_EPEC_quantity_ESS_fav.PS,s=60, label = "ESS_fav")
ax.set_xlabel("Social Welfare")
ax.set_ylabel("Producers Surplus")
ax.legend()
plt.show()
plt.grid()




sns = pyimport("seaborn")
sns.set()
sns.set_style("ticks")
sns.set_context("paper", font_scale=1.0)
# data = [simu_EPEC_quantity_empty.SW, simu_EPEC_quantity_competitive.SW, simu_EPEC_quantity_collusive.SW, simu_EPEC_quantity_ESS_fav.SW]
fig, ax = plt.subplots()
ax.boxplot([simu_EPEC_quantity_empty.SW, simu_EPEC_quantity_empty_confirmed.SW], positions = [1,2], widths = 0.9)
ax.boxplot([simu_EPEC_quantity_competitive.SW, simu_EPEC_quantity_competitive_confirmed.SW], positions = [4,5], widths = 0.9)
ax.boxplot([simu_EPEC_quantity_collusive.SW, simu_EPEC_quantity_collusive_confirmed.SW], positions = [7,8], widths = 0.9)
ax.boxplot([simu_EPEC_quantity_ESS_fav.SW, simu_EPEC_quantity_ESS_fav_confirmed.SW], positions = [10,11], widths = 0.9)
plt.grid()
ax.set_xticklabels(nms)
ax.set_xticks([1.5, 4.5, 7.5, 10.5])
# ax.set_xlabel("Social Planner objectives")
ax.set_ylabel("Social Welfare")
plt.savefig("boxplot_SW",bbox_inches="tight",dpi=600)


# data = [simu_EPEC_quantity_empty_confirmed.SW, simu_EPEC_quantity_competitive_confirmed.SW, simu_EPEC_quantity_collusive_confirmed.SW, simu_EPEC_quantity_ESS_fav_confirmed.SW]
# fig, ax = plt.subplots()
# ax.boxplot(data, widths=0.3, whis=[5,95], showfliers=true)
# plt.grid()
# ax.set_xticklabels(nms)
# # ax.set_xlabel("Social Planner objectives")
# ax.set_ylabel("Social Welfare")
# plt.savefig("boxplot_SW_conf",bbox_inches="tight",dpi=600)


fig, ax = plt.subplots()
ax.boxplot([simu_EPEC_quantity_empty.PS, simu_EPEC_quantity_empty_confirmed.PS], positions = [1,2], widths = 0.9)
ax.boxplot([simu_EPEC_quantity_competitive.PS, simu_EPEC_quantity_competitive_confirmed.PS], positions = [4,5], widths = 0.9)
ax.boxplot([simu_EPEC_quantity_collusive.PS, simu_EPEC_quantity_collusive_confirmed.PS], positions = [7,8], widths = 0.9)
ax.boxplot([simu_EPEC_quantity_ESS_fav.PS, simu_EPEC_quantity_ESS_fav_confirmed.PS], positions = [10,11], widths = 0.9)
plt.grid()
ax.set_xticklabels(nms)
ax.set_xticks([1.5, 4.5, 7.5, 10.5])
ax.set_ylabel("Producer Surplus")
plt.savefig("boxplot_PS",bbox_inches="tight",dpi=600)


###########PS BLOXPLOT####################
fig, ax = plt.subplots()
font_size = 15
plt.rc("font", size=font_size)
plt.rc("xtick", labelsize=font_size)
plt.rc("ytick", labelsize=font_size)
color_a = "black"
box1 = ax.boxplot([simu_EPEC_quantity_empty.PS/ avg_PS_empty], positions = [1], widths = 0.9)
plt.setp(box1["boxes"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box1["caps"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box1["whiskers"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box1["medians"], color="orange", linestyle="-", linewidth =2.0)

box2 = ax.boxplot([simu_EPEC_quantity_empty_confirmed.PS/ avg_PS_empty_confirmed], positions = [2], widths = 0.9)
plt.setp(box2["boxes"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box2["caps"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box2["whiskers"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box2["medians"], color="orange", linestyle="-", linewidth =2.0)

box3 = ax.boxplot([simu_EPEC_quantity_competitive.PS / avg_PS_empty], positions = [4], widths = 0.9)
plt.setp(box3["boxes"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box3["caps"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box3["whiskers"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box3["medians"], color="orange", linestyle="-", linewidth =2.0)

box4 = ax.boxplot([simu_EPEC_quantity_competitive_confirmed.PS / avg_PS_empty_confirmed], positions = [5], widths = 0.9)
plt.setp(box4["boxes"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box4["caps"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box4["whiskers"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box4["medians"], color="orange", linestyle="-", linewidth =2.0)

box5 = ax.boxplot([simu_EPEC_quantity_collusive.PS / avg_PS_empty ], positions = [7], widths = 0.9)
plt.setp(box5["boxes"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box5["caps"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box5["whiskers"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box5["medians"], color="orange", linestyle="-", linewidth =2.0)

box6 = ax.boxplot([simu_EPEC_quantity_collusive_confirmed.PS / avg_PS_empty_confirmed], positions = [8], widths = 0.9)
plt.setp(box6["boxes"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box6["caps"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box6["whiskers"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box6["medians"], color="orange", linestyle="-", linewidth =2.0)

box7 = ax.boxplot([simu_EPEC_quantity_ESS_fav.PS/ avg_PS_empty], positions = [10], widths = 0.9)
plt.setp(box7["boxes"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box7["caps"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box7["whiskers"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box7["medians"], color="orange", linestyle="-", linewidth =2.0)

box8 = ax.boxplot([simu_EPEC_quantity_ESS_fav_confirmed.PS/ avg_PS_empty_confirmed], positions = [11], widths = 0.9)
plt.setp(box8["boxes"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box8["caps"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box8["whiskers"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box8["medians"], color="orange", linestyle="-", linewidth =2.0)

plt.grid()
ax.set_xticklabels(nms)
ax.set_xticks([1.5, 4.5, 7.5, 10.5])
# ax.set_ylabel("Producer Surplus")
plt.savefig("boxplot_PS",bbox_inches="tight",dpi=600)


###########SW BLOXPLOT####################
fig, ax = plt.subplots()
font_size = 15
plt.rc("font", size=font_size)
plt.rc("xtick", labelsize=font_size)
plt.rc("ytick", labelsize=font_size)
color_a = "black"
box1 = ax.boxplot([simu_EPEC_quantity_empty.SW/ avg_SW_empty], positions = [1], widths = 0.9)
plt.setp(box1["boxes"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box1["caps"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box1["whiskers"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box1["medians"], color="orange", linestyle="-", linewidth =2.0)

box2 = ax.boxplot([simu_EPEC_quantity_empty_confirmed.SW/ avg_SW_empty_confirmed], positions = [2], widths = 0.9)
plt.setp(box2["boxes"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box2["caps"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box2["whiskers"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box2["medians"], color="orange", linestyle="-", linewidth =2.0)

box3 = ax.boxplot([simu_EPEC_quantity_competitive.SW / avg_SW_empty], positions = [4], widths = 0.9)
plt.setp(box3["boxes"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box3["caps"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box3["whiskers"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box3["medians"], color="orange", linestyle="-", linewidth =2.0)

box4 = ax.boxplot([simu_EPEC_quantity_competitive_confirmed.SW / avg_SW_empty_confirmed], positions = [5], widths = 0.9)
plt.setp(box4["boxes"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box4["caps"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box4["whiskers"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box4["medians"], color="orange", linestyle="-", linewidth =2.0)

box5 = ax.boxplot([simu_EPEC_quantity_collusive.SW / avg_SW_empty ], positions = [7], widths = 0.9)
plt.setp(box5["boxes"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box5["caps"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box5["whiskers"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box5["medians"], color="orange", linestyle="-", linewidth =2.0)

box6 = ax.boxplot([simu_EPEC_quantity_collusive_confirmed.SW / avg_SW_empty_confirmed], positions = [8], widths = 0.9)
plt.setp(box6["boxes"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box6["caps"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box6["whiskers"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box6["medians"], color="orange", linestyle="-", linewidth =2.0)

box7 = ax.boxplot([simu_EPEC_quantity_ESS_fav.SW/ avg_SW_empty], positions = [10], widths = 0.9)
plt.setp(box7["boxes"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box7["caps"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box7["whiskers"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box7["medians"], color="orange", linestyle="-", linewidth =2.0)

box8 = ax.boxplot([simu_EPEC_quantity_ESS_fav_confirmed.SW/ avg_SW_empty_confirmed], positions = [11], widths = 0.9)
plt.setp(box8["boxes"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box8["caps"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box8["whiskers"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box8["medians"], color="orange", linestyle="-", linewidth =2.0)

plt.grid()
ax.set_xticklabels(nms)
ax.set_xticks([1.5, 4.5, 7.5, 10.5])
# plt.plot([], c=color_a, linestyle="--", linewidth =1.5,  label="not-confirmed")
# plt.plot([], c=color_a, linestyle="-", linewidth =1.5,  label="confirmed")
# plt.legend()
plt.savefig("boxplot_SW",bbox_inches="tight",dpi=600)



###########NORMALIZED PROFIT DEVIATIONS#############
fig, ax = plt.subplots()
font_size = 15
plt.rc("font", size=font_size)
plt.rc("xtick", labelsize=font_size)
plt.rc("ytick", labelsize=font_size)
color_a = "black"
color_b = "blue"
box1 = ax.boxplot(norms_empty, positions = [1], widths = 0.9, showfliers=true)
plt.setp(box1["boxes"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box1["caps"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box1["whiskers"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box1["medians"], color="orange", linestyle="-", linewidth =2.0)

box2 = ax.boxplot(norms_empty_confirmed, positions = [2], widths = 0.9)
plt.setp(box2["boxes"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box2["caps"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box2["whiskers"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box2["medians"], color="orange", linestyle="-", linewidth =2.0)

box3 = ax.boxplot(norms_competitive, positions = [4], widths = 0.9, showfliers=true)
plt.setp(box3["boxes"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box3["caps"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box3["whiskers"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box3["medians"], color="orange", linestyle="-", linewidth =2.0)

box4 = ax.boxplot(norms_competitive_confirmed, positions = [5], widths = 0.9)
plt.setp(box4["boxes"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box4["caps"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box4["whiskers"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box4["medians"], color="orange", linestyle="-", linewidth =2.0)

box5 = ax.boxplot(norms_collusive, positions = [7], widths = 0.9, showfliers=true)
plt.setp(box5["boxes"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box5["caps"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box5["whiskers"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box5["medians"], color="orange", linestyle="-", linewidth =2.0)

box6 = ax.boxplot(norms_collusive_confirmed, positions = [8], widths = 0.9)
plt.setp(box6["boxes"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box6["caps"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box6["whiskers"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box6["medians"], color="orange", linestyle="-", linewidth =2.0)

box7 = ax.boxplot(norms_ESS_fav, positions = [10], widths = 0.9, showfliers=true)
plt.setp(box7["boxes"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box7["caps"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box7["whiskers"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box7["medians"], color="orange", linestyle="-", linewidth =2.0)

box8 = ax.boxplot(norms_ESS_fav_confirmed, positions = [11], widths = 0.9)
plt.setp(box8["boxes"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box8["caps"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box8["whiskers"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box8["medians"], color="orange", linestyle="-", linewidth =2.0)

plt.grid()
ax.set_xticklabels(nms)
ax.set_xticks([1.5, 4.5, 7.5, 10.5])
plt.plot([], c=color_a, linestyle="--", linewidth =1.5,  label="not-confirmed")
plt.plot([], c=color_a, linestyle="-", linewidth =1.5,  label="confirmed")
plt.legend()

# ax.set_ylabel("Normalized deviation in strategic agent's profit")
plt.savefig("boxplot_SA_profits",bbox_inches="tight",dpi=600)


###########ESS NORMALIZED PROFIT DEVIATION#############
fig, ax = plt.subplots()
font_size = 15
plt.rc("font", size=font_size)
plt.rc("xtick", labelsize=font_size)
plt.rc("ytick", labelsize=font_size)
color_a = "black"
color_b = "blue"
box1 = ax.boxplot(norm_profit_ESS_empty, positions = [1], widths = 0.9, showfliers=true)
plt.setp(box1["boxes"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box1["caps"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box1["whiskers"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box1["medians"], color="orange", linestyle="-", linewidth =2.0)

box2 = ax.boxplot(norm_profit_ESS_empty_confirmed, positions = [2], widths = 0.9)
plt.setp(box2["boxes"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box2["caps"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box2["whiskers"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box2["medians"], color="orange", linestyle="-", linewidth =2.0)

box3 = ax.boxplot(norm_profit_ESS_competitive, positions = [4], widths = 0.9, showfliers=true)
plt.setp(box3["boxes"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box3["caps"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box3["whiskers"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box3["medians"], color="orange", linestyle="-", linewidth =2.0)

box4 = ax.boxplot(norm_profit_ESS_competitive_confirmed, positions = [5], widths = 0.9)
plt.setp(box4["boxes"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box4["caps"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box4["whiskers"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box4["medians"], color="orange", linestyle="-", linewidth =2.0)

box5 = ax.boxplot(norm_profit_ESS_collusive, positions = [7], widths = 0.9, showfliers=true)
plt.setp(box5["boxes"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box5["caps"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box5["whiskers"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box5["medians"], color="orange", linestyle="-", linewidth =2.0)

box6 = ax.boxplot(norm_profit_ESS_collusive_confirmed, positions = [8], widths = 0.9)
plt.setp(box6["boxes"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box6["caps"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box6["whiskers"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box6["medians"], color="orange", linestyle="-", linewidth =2.0)

box7 = ax.boxplot(norm_profit_ESS_ESS_fav, positions = [10], widths = 0.9, showfliers=true)
plt.setp(box7["boxes"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box7["caps"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box7["whiskers"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box7["medians"], color="orange", linestyle="-", linewidth =2.0)

box8 = ax.boxplot(norm_profit_ESS_ESS_fav_confirmed, positions = [11], widths = 0.9)
plt.setp(box8["boxes"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box8["caps"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box8["whiskers"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box8["medians"], color="orange", linestyle="-", linewidth =2.0)

plt.grid()
ax.set_xticklabels(nms)
ax.set_xticks([1.5, 4.5, 7.5, 10.5])
plt.plot([], c=color_a, linestyle="--", linewidth =1.5,  label="not-confirmed")
plt.plot([], c=color_a, linestyle="-", linewidth =1.5,  label="confirmed")
plt.legend()

# ax.set_ylabel("Normalized deviation in strategic agent's profit")
plt.savefig("boxplot_ESS_profits",bbox_inches="tight",dpi=600)


###########GEN NORMALIZED PROFIT DEVIATION#############
fig, ax = plt.subplots()
font_size = 15
plt.rc("font", size=font_size)
plt.rc("xtick", labelsize=font_size)
plt.rc("ytick", labelsize=font_size)
color_a = "black"
color_b = "blue"
box1 = ax.boxplot(norm_profit_GEN_empty, positions = [1], widths = 0.9, showfliers=true)
plt.setp(box1["boxes"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box1["caps"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box1["whiskers"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box1["medians"], color="orange", linestyle="-", linewidth =2.0)

box2 = ax.boxplot(norm_profit_GEN_empty_confirmed, positions = [2], widths = 0.9)
plt.setp(box2["boxes"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box2["caps"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box2["whiskers"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box2["medians"], color="orange", linestyle="-", linewidth =2.0)

box3 = ax.boxplot(norm_profit_GEN_competitive, positions = [4], widths = 0.9, showfliers=true)
plt.setp(box3["boxes"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box3["caps"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box3["whiskers"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box3["medians"], color="orange", linestyle="-", linewidth =2.0)

box4 = ax.boxplot(norm_profit_GEN_competitive_confirmed, positions = [5], widths = 0.9)
plt.setp(box4["boxes"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box4["caps"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box4["whiskers"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box4["medians"], color="orange", linestyle="-", linewidth =2.0)

box5 = ax.boxplot(norm_profit_GEN_collusive, positions = [7], widths = 0.9, showfliers=true)
plt.setp(box5["boxes"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box5["caps"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box5["whiskers"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box5["medians"], color="orange", linestyle="-", linewidth =2.0)

box6 = ax.boxplot(norm_profit_GEN_collusive_confirmed, positions = [8], widths = 0.9)
plt.setp(box6["boxes"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box6["caps"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box6["whiskers"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box6["medians"], color="orange", linestyle="-", linewidth =2.0)

box7 = ax.boxplot(norm_profit_GEN_ESS_fav, positions = [10], widths = 0.9, showfliers=true)
plt.setp(box7["boxes"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box7["caps"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box7["whiskers"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box7["medians"], color="orange", linestyle="-", linewidth =2.0)

box8 = ax.boxplot(norm_profit_GEN_ESS_fav_confirmed, positions = [11], widths = 0.9)
plt.setp(box8["boxes"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box8["caps"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box8["whiskers"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box8["medians"], color="orange", linestyle="-", linewidth =2.0)

plt.grid()
ax.set_xticklabels(nms)
ax.set_xticks([1.5, 4.5, 7.5, 10.5])
plt.plot([], c=color_a, linestyle="--", linewidth =1.5,  label="not-confirmed")
plt.plot([], c=color_a, linestyle="-", linewidth =1.5,  label="confirmed")
plt.legend()

# ax.set_ylabel("Normalized deviation in strategic agent's profit")
plt.savefig("boxplot_GEN_profits",bbox_inches="tight",dpi=600)

###########WIND NORMALIZED PROFIT DEVIATION#############
fig, ax = plt.subplots()
font_size = 15
plt.rc("font", size=font_size)
plt.rc("xtick", labelsize=font_size)
plt.rc("ytick", labelsize=font_size)
color_a = "black"
color_b = "blue"
box1 = ax.boxplot(norm_profit_WIND_empty, positions = [1], widths = 0.9, showfliers=true)
plt.setp(box1["boxes"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box1["caps"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box1["whiskers"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box1["medians"], color="orange", linestyle="-", linewidth =2.0)

box2 = ax.boxplot(norm_profit_WIND_empty_confirmed, positions = [2], widths = 0.9)
plt.setp(box2["boxes"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box2["caps"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box2["whiskers"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box2["medians"], color="orange", linestyle="-", linewidth =2.0)

box3 = ax.boxplot(norm_profit_WIND_competitive, positions = [4], widths = 0.9, showfliers=true)
plt.setp(box3["boxes"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box3["caps"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box3["whiskers"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box3["medians"], color="orange", linestyle="-", linewidth =2.0)

box4 = ax.boxplot(norm_profit_WIND_competitive_confirmed, positions = [5], widths = 0.9)
plt.setp(box4["boxes"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box4["caps"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box4["whiskers"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box4["medians"], color="orange", linestyle="-", linewidth =2.0)

box5 = ax.boxplot(norm_profit_WIND_collusive, positions = [7], widths = 0.9, showfliers=true)
plt.setp(box5["boxes"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box5["caps"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box5["whiskers"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box5["medians"], color="orange", linestyle="-", linewidth =2.0)

box6 = ax.boxplot(norm_profit_WIND_collusive_confirmed, positions = [8], widths = 0.9)
plt.setp(box6["boxes"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box6["caps"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box6["whiskers"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box6["medians"], color="orange", linestyle="-", linewidth =2.0)

box7 = ax.boxplot(norm_profit_WIND_ESS_fav, positions = [10], widths = 0.9, showfliers=true)
plt.setp(box7["boxes"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box7["caps"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box7["whiskers"], color=color_a,  linestyle="--", linewidth =2.0)
plt.setp(box7["medians"], color="orange", linestyle="-", linewidth =2.0)

box8 = ax.boxplot(norm_profit_WIND_ESS_fav_confirmed, positions = [11], widths = 0.9)
plt.setp(box8["boxes"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box8["caps"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box8["whiskers"], color=color_a,  linestyle="-", linewidth =2.0)
plt.setp(box8["medians"], color="orange", linestyle="-", linewidth =2.0)

plt.grid()
ax.set_xticklabels(nms)
ax.set_xticks([1.5, 4.5, 7.5, 10.5])
plt.plot([], c=color_a, linestyle="--", linewidth =1.5,  label="not-confirmed")
plt.plot([], c=color_a, linestyle="-", linewidth =1.5,  label="confirmed")
plt.legend()

# ax.set_ylabel("Normalized deviation in strategic agent's profit")
plt.savefig("boxplot_WIND_profits",bbox_inches="tight",dpi=600)





data = [simu_EPEC_quantity_empty_confirmed.PS, simu_EPEC_quantity_competitive_confirmed.PS, simu_EPEC_quantity_collusive_confirmed.PS, simu_EPEC_quantity_ESS_fav_confirmed.PS]
fig, ax = plt.subplots()
ax.boxplot(data, widths=0.3, whis=[5,95], showfliers=true)
plt.grid()
ax.set_xticklabels(nms)
# ax.set_xlabel("Social Planner objectives")
ax.set_ylabel("Producer Surplus")
plt.savefig("boxplot_PS_conf",bbox_inches="tight",dpi=600)


using StatsPlots
using Plots

# using LaTeXStrings
theme(:solarized_light)

@df df Plots.boxplot(["empty"], simu_EPEC_quantity_empty.SW[:], fmt = :png,  legend = false, dpi=6000)
@df df Plots.boxplot!(["competitive"], simu_EPEC_quantity_competitive.SW[:], fmt = :png,  legend = false, dpi=600)
@df df Plots.boxplot!(["collusive"], simu_EPEC_quantity_collusive.SW[:], fmt = :png,  legend = false, dpi=600)
@df df Plots.boxplot!(["ESS_fav"], simu_EPEC_quantity_ESS_fav.SW[:], fmt = :png,  legend = false, dpi=600)

Plots.savefig("SW_box")



using Pandas
df_2 = Pandas.DataFrame()
df_2["variables"] = df[:,:variables]

sns = pyimport("seaborn")
sns.set()
sns.set_style("ticks")
sns.set_context("paper", font_scale=1.0)

sns.boxplot(data = df)

df = DataFrames.stack(simu_EPEC_quantity_empty[:,3:5])
df[:, :variable] = string(df[:, :variable])
