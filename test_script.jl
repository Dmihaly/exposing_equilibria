using Revise
using JuMP
using Ipopt
using PyPlot
using DataFrames
using StructArrays
using CSV
using LinearAlgebra
using KNITRO
using Gurobi


include("MPPDCs/MPPDC_quantity.jl") #The MPPDC model with quantity bidding
include("MPECs/MPEC_quantity.jl") #The MPEC model with quantity bidding
include("MILPs/MILP_quantity.jl") #The MILP model with quantity bidding
include("MPPDCs/NLPs/NLP_primal_MPPDC_quantity.jl") #The NLP model with quantity bidding
include("MPPDCs/NLPs/NLP_dual_MPPDC_quantity_mod.jl") #The NLP model with quantity bidding
include("MPPDCs/NLPs/NLP_dual_MPPDC_quantity.jl") #The NLP model with quantity bidding

###############DECLARE PARTICIPATING AGENTS###########

cd("C:/Users/u0125813/Desktop/Leuven/KU_Leuven/EPECs/scr_nlp/DATA/")

generators_dat = CSV.read("generators_data.csv", header=true, delim= ";")
MC_g = generators_dat[1,2:8]
MC_g_str = generators_dat[1,9]
#CAPACITIES OF GENERATORS
capacity_g = generators_dat[2,2:9]
capacity_g = values(capacity_g) ./ 100
capacity_g_nonstr = values(capacity_g[1:7])
capacity_g_str = values(capacity_g[8])

#####SIMPLIFIED VERSION OF GENERATOR DATA######
# generators_dat = CSV.read("generators_data_simpl.csv", header=true, delim= ";")
# MC_g = generators_dat[1,2:4]
# MC_g_str = generators_dat[1,5]
# #CAPACITIES OF GENERATORS
# capacity_g = generators_dat[2,2:5]./100
# capacity_g_nonstr = values(capacity_g[1:3])
# capacity_g_str = values(capacity_g[4])


#RAMPING LIMIT OF THE STRATEGIC GENERATOR:
Rᴳ = 1.0

res_dat = CSV.read("wind_PV.csv", header=true, delim= ";")
eight_res  = res_dat[1:3:end,:] #selecting every 3rd of the res_dat
renewable_nonstr = ["PV"]
renewable_str = ["WIND"]
capacity_re_nonstr = eight_res.PV/100
capacity_re_str =  eight_res.Wind/100

#LOADING THE DEMAND
load_dat = CSV.read("load.csv",header=true, delim= ";")
load = load_dat.load[1:3:end,:]/100


eff_ch = eff_dch = generators_dat.ESS[end]
energy_max = generators_dat.ESS[end-2]/100
power_max = generators_dat.ESS[end-4]/100

MC_re_str = MC_re_nonstr = 0 #Zero marginal cost

str_stor = ["ESS_STR"]
str_gens = ["GEN_STR", "WIND_STR"]
str_agents = vcat(str_stor, str_gens)

generators = names(MC_g)
renewables = [:PV]

non_str_gens = vcat(generators)
non_str_res = vcat(renewables)
non_str_agents = vcat(generators, renewables)

#NOTE:ATM There is no non-str ESS involved in the pool
#MARKET PRICE LIMITS
p_cap = 100
p_floor = 0


function t_interval(lenght, resolution)
    t_steps = [i for i in range(1, stop = lenght, step = resolution)]
    return t_steps
end

t_steps = t_interval(8,1)

function compute_profits()
    lambda = value.(NLP[:λ][:])
    #NON-STRATEGIC GENERATORS
    gen_1 = value.(NLP[:g]["GEN1",:])
    gen_2 = value.(NLP[:g]["GEN2",:])
    gen_3 = value.(NLP[:g]["GEN3",:])
    gen_5 = value.(NLP[:g]["GEN5",:])
    gen_6 = value.(NLP[:g]["GEN6",:])
    gen_7 = value.(NLP[:g]["GEN7",:])
    gen_8 = value.(NLP[:g]["GEN8",:])

    #DEMAND
    d = value.(NLP[:d])
    #NON-STRATEGIC RES
    res_PV = value.(NLP[:w][:PV,:])
    #DISPATCHED QUANTITIES FROMS STR AGENTS
    dch = value.(NLP[:dch_st])
    ch = value.(NLP[:ch_st])
    gen = value.(NLP[:g_st])
    res = value.(NLP[:w_st])
    #CALCULATING PROFITS
    earned_ESS = sum(dch[t] * lambda[t] for t in t_steps)
    paid_ESS = sum(ch[t] * lambda[t] for t in t_steps)
    profit_ESS =  earned_ESS - paid_ESS
    # println("PROFIT_ESS: ", profit_ESS)

    earned_GEN = sum(gen[t] * lambda[t] for t in t_steps)
    paid_GEN = sum(dat_str.GEN_STR.MC_g[t] * gen[t] for t in t_steps)
    profit_GEN =  earned_GEN - paid_GEN
    # println("PROFIT_GEN: ", profit_GEN)

    earned_RES = sum(res[t] * lambda[t] for t in t_steps)
    paid_RES = sum(dat_str.WIND_STR.MC_g[t] * res[t] for t in t_steps)
    profit_RES =  earned_RES - paid_RES
    # println("PROFIT_RES: ", profit_RES)

    return profit_ESS, profit_GEN, profit_RES
end
function compute_profits_MPPDC()
    lambda = value.(MPPDC[:λ][:])
    #NON-STRATEGIC GENERATORS
    gen_1 = value.(MPPDC[:g]["GEN1",:])
    gen_2 = value.(MPPDC[:g]["GEN2",:])
    gen_3 = value.(MPPDC[:g]["GEN3",:])
    gen_5 = value.(MPPDC[:g]["GEN5",:])
    gen_6 = value.(MPPDC[:g]["GEN6",:])
    gen_7 = value.(MPPDC[:g]["GEN7",:])
    gen_8 = value.(MPPDC[:g]["GEN8",:])

    #DEMAND
    d = value.(MPPDC[:d])
    #NON-STRATEGIC RES
    res_PV = value.(MPPDC[:w][:PV,:])
    #DISPATCHED QUANTITIES FROMS STR AGENTS
    dch = value.(MPPDC[:dch_st])
    ch = value.(MPPDC[:ch_st])
    gen = value.(MPPDC[:g_st])
    res = value.(MPPDC[:w_st])
    #CALCULATING PROFITS
    earned_ESS = sum(dch[t] * lambda[t] for t in t_steps)
    paid_ESS = sum(ch[t] * lambda[t] for t in t_steps)
    profit_ESS =  earned_ESS - paid_ESS
    # println("PROFIT_ESS: ", profit_ESS)

    earned_GEN = sum(gen[t] * lambda[t] for t in t_steps)
    paid_GEN = sum(dat_str.GEN_STR.MC_g[t] * gen[t] for t in t_steps)
    profit_GEN =  earned_GEN - paid_GEN
    # println("PROFIT_GEN: ", profit_GEN)

    earned_RES = sum(res[t] * lambda[t] for t in t_steps)
    paid_RES = sum(dat_str.WIND_STR.MC_g[t] * res[t] for t in t_steps)
    profit_RES =  earned_RES - paid_RES
    # println("PROFIT_RES: ", profit_RES)

    return profit_ESS, profit_GEN, profit_RES
end


###########INIT BID, TECH LIMIT DFs#########
function return_data()
    #DATA OF THE STR AGENTS, DIFFERENTIATE BETWEEN BIDS AND TRUE VALUES
    dat_str = DataFrame()
    for i in range(1, stop = length(str_gens)), l in range(1, stop = length(str_stor))
        #This needs to be revisited....
        if i ==  1 #if str agent is GEN
            dat_str[i] =  StructArray(
            P_max = ones(length(t_steps)) * capacity_g_str,
            G = zeros(length(t_steps)),
            c_g = zeros(length(t_steps)),
            MC_g = ones(length(t_steps)) * MC_g_str)
        elseif i == 2 #if str agent is RE
            dat_str[i] =  StructArray(
            P_max = capacity_re_str,
            G = zeros(length(t_steps)),
            c_g = zeros(length(t_steps)),
            MC_g = zeros(length(t_steps)))
        end

        dat_str[i + l] =  StructArray(
        E_max = ones(length(t_steps)) * energy_max,
        P_max = ones(length(t_steps)) * power_max,
        DCH = zeros(length(t_steps)),
        MC_dch = ones(length(t_steps)) * p_cap,
        c_dch = zeros(length(t_steps)),
        CH = zeros(length(t_steps)),
        c_ch = zeros(length(t_steps)),
        MC_ch = ones(length(t_steps)) * p_floor)
    end

    #NOTE: For non-str agents set the quantity bid equal to the maximum capacity
    #NON STR AGENTS THE BIDS ARE EQUAL TO THE TRUE VALUES
    dat_nonstr = DataFrame()
    for i in range(1, stop = length(non_str_gens)), k in range(1, stop = length(non_str_res))
        dat_nonstr[i] =  StructArray(
        G = ones(length(t_steps)) * capacity_g_nonstr[i],
        c_g = ones(length(t_steps)) * MC_g[1+i-1])
        dat_nonstr[i + k] = StructArray(
        G = capacity_re_nonstr,
        c_g = ones(length(t_steps)) * MC_re_nonstr #as zero marginal cost
        )
    end

    #NOW RENAME THEM
    names!(dat_str, vcat(Symbol.(str_gens), Symbol.(str_stor)))
    names!(dat_nonstr, Symbol.(non_str_agents))

    return dat_str, dat_nonstr
end

dat_str, dat_nonstr = return_data()
#
optimizer = (with_optimizer(Ipopt.Optimizer,
                                max_iter = 3000,
                                print_level = 0
                                # tol = 1.0e-16,
                                # acceptable_tol = 1.0e-10
                                ))

optimizer = with_optimizer(KNITRO.Optimizer,
                # algorithm = 1,
                hessian_no_f=1,
                datacheck=0,
                outlev=0,
                # ms_enable= 1,
                # ms_maxsolves = 50
                )

# optimizer = with_optimizer(Gurobi.Optimizer)


# optimizer = with_optimizer(Gurobi.Optimizer, NonConvex = 2)

#SETTING THE CURRENT AGENT
#NOTE: ACCEPTED STRINGS: "ESS_STR", "GEN_STR", "WIND_STR"
current = "ESS_STR"
bidding = "quantity"
model_type = "MPEC"
model_type = "MPPDC"
model_type = "MILP"
objective_type = "dual"
flag_regularize = true

MPPDC = build_MPPDC_quantity(t_steps, non_str_gens,
non_str_res, current, optimizer, objective_type)
# optimize!(MPPDC)
# objective_value(MPPDC)

println("SOLUTION METHOD:
                ....Agent: $current
                ....Objective_type: $objective_type")
#
MPPDC = build_MILP_quantity(t_steps, non_str_gens,
non_str_res, current, optimizer, objective_type)
# JuMP.fix(MPPDC[:M1], 100; force = true)
# JuMP.fix(MPPDC[:M1], 100; force = true)
NLP = build_NLP_dual_quantity_mod(t_steps, non_str_gens,
non_str_res, str_agents, current, optimizer)

# optimize!(NLP)
# objective_value(NLP)

# JuMP.fix(MPPDC[:ϵ], 0.00001; force = true)
optimize!(MPPDC)
objective_value(MPPDC)
profit_ESS, profit_GEN, profit_RES = compute_profits_MPPDC()


flag_sequanctial_red = true
#NOTE: MIND THE STEP SIZE
results_MPPDC_ESS = DataFrame(ϵ = [], iter = [], profit = [])
results_MPPDC_GEN = DataFrame(ϵ = [], iter = [], profit = [])
results_MPPDC_WIND = DataFrame(ϵ = [], iter = [], profit = [])

#FOR MPPDC
if flag_sequanctial_red == true
    #Rebuild the problem
    MPPDC = build_MPPDC_quantity(t_steps, non_str_gens,
    non_str_res, current, optimizer, objective_type)
    global iter_limit = 100
    global iterations = 0
    global tol = 0.001
    global epsilon_new = 1000
    JuMP.fix(MPPDC[:ϵ], epsilon_new; force = true)

    optimize!(MPPDC)
    objective_value(MPPDC)
    profit_ESS, profit_GEN, profit_RES = compute_profits_MPPDC()
    println("PROFITS:
                    PROFIT_ESS: $profit_ESS,
                    PROFIT_GEN: $profit_GEN,
                    PROFIT_RES: $profit_RES")
    global iterations += 1
    push!(results_MPPDC_ESS, [epsilon_new, iterations, profit_ESS])
    ############THE SMAL SCRIPT SNIPPET FOR UPDATING THE REGULARIZATION TERM
    while epsilon_new >= tol && iterations < iter_limit
        #########SETTING THE WARM START POINT TO THE PREVIOUS SOLUTION#########
        #if terminition status is local optimal, otherwise should not use the solution as starting point
        term_stat = termination_status(MPPDC)
        # if term_stat == MOI.OPTIMAL || term_stat == MOI.LOCALLY_SOLVED
        set_start_value.(all_variables(MPPDC), value.(all_variables(MPPDC)))
        ###################FIXING EPRSILONS TO THE REDUCED PARAMS##############
        global epsilon_new = value.(MPPDC[:ϵ])*0.1
        # else
        #     set_start_value.(all_variables(MPPDC), prev_vals)
        #     global epsilon_new = value.(MPPDC[:ϵ])*0.99
        # end
        fix(MPPDC[:ϵ], epsilon_new; force = true)
        println("\\
        ######### \\
        THE NEW RELAXATION TERMS: \\
                    ϵ: $epsilon_new and \\
                    #ITERATIONS: $iterations \\
                    #TERMINATION: $term_stat
        #########")
        ################OPTIMIZING############
        optimize!(MPPDC)
        profit_ESS, profit_GEN, profit_RES = compute_profits_MPPDC()
        println("PROFITS:
                        PROFIT_ESS: $profit_ESS,
                        PROFIT_GEN: $profit_GEN,
                        PROFIT_RES: $profit_RES")

        ##############CACLCULATING PROFIT############
        # profit = compute_profit()
        # println("PROFIT: ", profit)
        ############PROGRESSING IN TIME##############
        global iterations += 1
        push!(results_MPPDC_ESS, [epsilon_new, iterations, profit_ESS])
        global prev_vals = value.(all_variables(MPPDC))
    end
end



using PyCall, PyPlot
pygui(true)
figure()
#plot here the line with the optimal
plt.plot(results_MPPDC_ESS.iter[1:5], results_MPPDC_ESS.profit[1:5], label="ϵ_init = 1", color ="k", linestyle="-", linewidth =1.5, marker="^", markersize = 12)
plt.plot(results_MPPDC_ESS.iter[6:11], results_MPPDC_ESS.profit[6:11], label="ϵ_init = 10", color ="k", linestyle="--", linewidth =1.5, marker="D", markersize = 10)
plt.plot(results_MPPDC_ESS.iter[12:18], results_MPPDC_ESS.profit[12:18], label="ϵ_init = 100", color ="k", linestyle=":", linewidth =1.5, marker="s", markersize = 10)
plt.plot(results_MPPDC_ESS.iter[19:26], results_MPPDC_ESS.profit[19:26], label="ϵ_init = 1000", color ="k", linestyle="-.", linewidth =1.5, marker="p", markersize = 10)
plt.plot(range(1,8,step = 1), ones(8).*235.579, label = "MILP", color = "r", linewidth = 2.5)
# plt.xticks(range(0.02, 0.24, step = 0.02))
plt.ylabel("Profit of ESS")
plt.xlabel("Iterations")
plt.legend()
plt.grid()


#FOR NLP
results_NLP_ESS = DataFrame(ϵ1 = [], ϵ2 = [], iter = [], profit = [])
results_NLP_GEN = DataFrame(ϵ1 = [], ϵ2 = [], iter = [], profit = [])
results_NLP_WIND = DataFrame(ϵ1 = [], ϵ2 = [], iter = [], profit = [])

if flag_sequanctial_red == true
    global iter_limit = 50
    global iterations = 0
    global tol = 0.001
    global epsilon1_new = 100
    global epsilon2_new = 100
    #Rebuild the model, otherwise there is some bias
    NLP = build_NLP_dual_quantity_mod(t_steps, non_str_gens,
    non_str_res, str_agents, current, optimizer)

    JuMP.fix(NLP[:ϵ1], epsilon1_new; force = true)
    JuMP.fix(NLP[:ϵ2], epsilon2_new; force = true)
    optimize!(NLP)
    objective_value(NLP)
    profit_ESS, profit_GEN, profit_RES = compute_profits()
    println("PROFITS:
                    PROFIT_ESS: $profit_ESS,
                    PROFIT_GEN: $profit_GEN,
                    PROFIT_RES: $profit_RES")


    global iterations += 1
    push!(results_NLP_ESS, [epsilon1_new, epsilon2_new, iterations, profit_ESS])
    ############THE SMAL SCRIPT SNIPPET FOR UPDATING THE REGULARIZATION TERM
    while epsilon1_new >= tol || epsilon2_new >= tol
        #After reaching iteration limit exit:
        if iterations >= iter_limit
            break
        end
        #########SETTING THE WARM START POINT TO THE PREVIOUS SOLUTION#########
        #if terminition status is local optimal, otherwise should not use the solution as starting point
        term_stat = termination_status(NLP)
        # if term_stat == MOI.OPTIMAL || term_stat == MOI.LOCALLY_SOLVED
        set_start_value.(all_variables(NLP), value.(all_variables(NLP)))
        ###################FIXING EPRSILONS TO THE REDUCED PARAMS##############
        global epsilon1_new = value.(NLP[:ϵ1])*0.1
        if epsilon2_new >= tol
            global epsilon2_new = value.(NLP[:ϵ2])*0.1
        end
        # else
        #     set_start_value.(all_variables(MPPDC), prev_vals)
        #     global epsilon_new = value.(MPPDC[:ϵ])*0.99
        # end
        fix(NLP[:ϵ1], epsilon1_new; force = true)
        fix(NLP[:ϵ2], epsilon2_new; force = true)

        println("
        #########
        THE NEW RELAXATION TERMS:
                    # ϵ1: $epsilon1_new and
                    ϵ2: $epsilon2_new and
                    #ITERATIONS: $iterations
                    #TERMINATION: $term_stat
        #########")
        ################OPTIMIZING############
        optimize!(NLP)
        ##############CACLCULATING PROFIT############
        profit_ESS, profit_GEN, profit_RES = compute_profits()
        println("PROFITS:
                        PROFIT_ESS: $profit_ESS,
                        PROFIT_GEN: $profit_GEN,
                        PROFIT_RES: $profit_RES")
        global iterations += 1
        push!(results_NLP_ESS, [epsilon1_new, epsilon2_new, iterations, profit_ESS])
        ############PROGRESSING IN TIME##############
        global prev_vals = value.(all_variables(NLP))
    end
end

figure()
#plot here the line with the optimal
plt.plot(results_NLP_ESS.iter[1:5], results_NLP_ESS.profit[1:5]./10, label="ϵ¹_init=1, ϵ²_init=1", color ="k", linestyle="-", linewidth =1.5, marker="^", markersize = 12)
plt.plot(results_NLP_ESS.iter[6:11], results_NLP_ESS.profit[6:11]./10, label="ϵ¹_init=1, ϵ²_init=10", color ="k", linestyle="--", linewidth =1.5, marker="D", markersize = 10)
plt.plot(results_NLP_ESS.iter[12:18], results_NLP_ESS.profit[12:18]./10, label="ϵ¹_init=1, ϵ²_init=100", color ="k", linestyle=":", linewidth =1.5, marker="s", markersize = 10)
# plt.plot(results_NLP_ESS.iter[19:26], results_NLP_ESS.profit[19:26], label="ϵ_init = 1/1000", color ="k", linestyle="-.", linewidth =1.5, marker="p", markersize = 10)
plt.plot(results_NLP_ESS.iter[19:25], results_NLP_ESS.profit[19:25]./10, label="ϵ¹_init=10, ϵ²_init=100", color ="k", linestyle="-.", linewidth =1.5, marker="*", markersize = 10)
plt.plot(results_NLP_ESS.iter[26:end], results_NLP_ESS.profit[26:end]./10, label="ϵ¹_init=100, ϵ²_init=100", color ="k", linestyle="-", linewidth =1.5, marker="o", markersize = 10)
# plt.plot(results_NLP_ESS.iter[41:48], results_NLP_ESS.profit[41:48], label="ϵ_init = 1000/1000", color ="k", linestyle=":", linewidth =1.5, marker="2", markersize = 10)
plt.plot(range(1,8,step = 1), ones(8).*235.579/10, label = "MILP", color = "r", linewidth = 2.5)
# plt.xticks(range(0.02, 0.24, step = 0.02))
plt.ylabel("Profit of ESS (k€)")
plt.xlabel("Number of iterations")
plt.legend()
plt.grid()
# plt.margins(0,0)

lambda = value.(MPPDC[:λ][:])
#NON-STRATEGIC GENERATORS
gen_1 = value.(MPPDC[:g][:GEN1,:])
gen_2 = value.(MPPDC[:g][:GEN2,:])
gen_3 = value.(MPPDC[:g][:GEN3,:])
gen_5 = value.(MPPDC[:g][:GEN5,:])
gen_6 = value.(MPPDC[:g][:GEN6,:])
gen_7 = value.(MPPDC[:g][:GEN7,:])
gen_8 = value.(MPPDC[:g][:GEN8,:])

#DEMAND
d = value.(MPPDC[:d])
#NON-STRATEGIC RES
res_PV = value.(MPPDC[:w][:PV,:])
#DISPATCHED QUANTITIES FROMS STR AGENTS
dch = value.(MPPDC[:dch_st])
ch = value.(MPPDC[:ch_st])
gen = value.(MPPDC[:g_st])
res = value.(MPPDC[:w_st])


# lambda = value.(NLP[:λ][:])
# #NON-STRATEGIC GENERATORS
# gen_1 = value.(NLP[:g][:GEN1,:])
# gen_2 = value.(NLP[:g][:GEN2,:])
# gen_3 = value.(NLP[:g][:GEN3,:])
# gen_5 = value.(NLP[:g][:GEN5,:])
# gen_6 = value.(NLP[:g][:GEN6,:])
# gen_7 = value.(NLP[:g][:GEN7,:])
# gen_8 = value.(NLP[:g][:GEN8,:])
#
# #DEMAND
# d = value.(NLP[:d])
# #NON-STRATEGIC RES
# res_PV = value.(NLP[:w][:PV,:])
# #DISPATCHED QUANTITIES FROMS STR AGENTS
# dch = value.(NLP[:dch_st])
# ch = value.(NLP[:ch_st])
# gen = value.(NLP[:g_st])
# res = value.(NLP[:w_st])

#READING THE BIDS OF UL AGENT
# #QUANTITY-PRICE BIDS
# #ESS
# e_ch_bid_strat = value.(MPPDC[:CH_ST])
# p_ch_bid_strat = dat_str.ESS_STR.MC_ch
# e_dch_bid_strat = value.(MPPDC[:DCH_ST])
# p_dch_bid_strat = dat_str.ESS_STR.MC_dch
# #GENCO
# p_g_bid_strat = dat_str.GEN_STR.MC_g
# e_g_bid_strat = dat_str.GEN_STR.G
# #RES
# e_r_bid_strat = p_g_bid_strat = dat_str.WIND_STR.G
# p_r_bid_strat = p_g_bid_strat = dat_str.WIND_STR.MC_g
#CALCULATING PROFITS
earned_ESS = sum(dch[t] * lambda[t] for t in t_steps)
paid_ESS = sum(ch[t] * lambda[t] for t in t_steps)
profit_ESS =  earned_ESS - paid_ESS
println("PROFIT_ESS: ", profit_ESS)

earned_GEN = sum(gen[t] * lambda[t] for t in t_steps)
paid_GEN = sum(dat_str.GEN_STR.MC_g[t] * gen[t] for t in t_steps)
profit_GEN =  earned_GEN - paid_GEN
println("PROFIT_GEN: ", profit_GEN)

earned_RES = sum(res[t] * lambda[t] for t in t_steps)
paid_RES = sum(dat_str.WIND_STR.MC_g[t] * res[t] for t in t_steps)
profit_RES =  earned_RES - paid_RES
println("PROFIT_RES: ", profit_RES)


plot_flag_EPEC = true
if plot_flag_EPEC == true
     #PLOTS HERE
     using PyCall, PyPlot
     pygui(true)
     sns = pyimport("seaborn")
     sns.set()
     color_palette = sns.color_palette("Paired",12)
     sns.set_palette(color_palette)
     # sns.set_style("ticks")
     sns.set_style("darkgrid")
     sns.set_context("paper", font_scale=3.0)
     fig = figure()
     # ax1 = subplot(3,1,1)
     stackplot(t_steps,
     gen_1[:], gen_2[:], gen_3[:],
     gen_5[:], gen_6[:], gen_7[:], gen_8[:],
     res_PV[:], dch[:], gen[:], res[:],
     edgecolor="k", linewidth = 2.0,
     labels =["GEN1", "GEN2", "GEN3",
     "GEN5", "GEN6", "GEN7", "GEN8",
     "PV", "DCH-ST", "GEN-ST", "WIND-ST"],
     alpha=0.7
     )
     stackplot(t_steps,edgecolor="k", linewidth = 1.5, -ch[:], labels =["CH-ST"], alpha=0.7)
     plot(t_steps, load[:], label="LOAD",color ="k",linestyle=":", linewidth =3.5,marker="D", markersize = 8)
     # grid()
     legend(bbox_to_anchor=(0.5, 1.19), loc = "upper center",ncol=5)
     x_stk = 0:1:8
     plt.xticks(x_stk)
     plt.xlim((1,8))
     plt.ylabel("Generation (MW)")
     plt.xlabel("Time (h)")
     # savefig("myplot.png")



     #Bidding prices
     sns = pyimport("seaborn")
     sns.set()
     sns.set_style("darkgrid")
     sns.set_context("paper")
     fig = figure()
     # ax2 = subplot(3,1,2)
     # plot(t_steps, p_dch_bid_strat[:], label="BID_DCH_STR", color = "tab:brown", drawstyle = "steps-mid", linestyle="--", linewidth =2.5)
     # plot(t_steps, p_ch_bid_strat[:], label="BID_CH_STR", color = "tab:pink", drawstyle = "steps-mid", linestyle="--", linewidth =2.5)
     # plot(t_steps, p_g_bid_strat[:], label="BID_GEN_STR", color = "tab:olive", drawstyle = "steps-mid", linestyle="--", linewidth =2.5)
     # plot(t_steps, p_r_bid_strat * ones((length(t_steps))), label="BID_RES_STR", color = "tab:cyan", drawstyle = "steps-mid", linestyle="--", linewidth =2.5)
     # plot(t_steps, p_g["mid"] * ones((length(t_steps))), color = "tab:orange", label="BID_GEN_MID", drawstyle = "steps-mid", linewidth =1.5)
     # plot(t_steps, p_g["base"] * ones((length(t_steps))), color = "tab:blue", label="BID_GEN_BASE", drawstyle = "steps-mid", linewidth =1.5)
     # plot(t_steps, p_g["peak"] * ones((length(t_steps))), color = "tab:green", label="BID_GEN_PEAK", drawstyle = "steps-mid", linewidth =1.5)
     # plot(t_steps, p_r["PV"] * ones((length(t_steps))), color = "tab:purple", label="BID_PV", drawstyle = "steps-mid",linewidth =1.5)
     # plot(t_steps, p_r["wind"] * ones((length(t_steps))), color = "tab:red", label="BID_WIND", drawstyle = "steps-mid",linewidth =1.5)
     plot(t_steps, lambda[:], label="Market Price EPEC", drawstyle = "steps-mid", marker="D",color ="k", markersize = 8, linestyle=":",)
     legend()
     title("PRICE BIDS")
     # #Bidding quantities
     # ax3 = subplot(3,1,3)
     # plot(t_steps, e_dch_bid_strat[:], label="BID_DCH_STR", color = "tab:brown", drawstyle = "steps-mid", linestyle="--", linewidth =2.5)
     # plot(t_steps, e_ch_bid_strat[:], label="BID_CH_STR", color = "tab:pink", drawstyle = "steps-mid", linestyle="--", linewidth =2.5)
     # plot(t_steps, e_g_bid_strat[:], label="BID_GEN_STR", color = "tab:olive", drawstyle = "steps-mid", linestyle="--", linewidth =2.5)
     # plot(t_steps, e_r_bid_strat[:], label="BID_RES_STR", color = "tab:cyan", drawstyle = "steps-mid", linestyle="--", linewidth =2.5)
     # plot(t_steps, P_g["mid"] * ones((length(t_steps))), color = "tab:orange", label="BID_GEN_MID", drawstyle = "steps-mid", linewidth =1.5)
     # plot(t_steps, P_g["base"] * ones((length(t_steps))), color = "tab:blue", label="BID_GEN_BASE", drawstyle = "steps-mid", linewidth =1.5)
     # plot(t_steps, P_g["peak"] * ones((length(t_steps))), color = "tab:green", label="BID_GEN_PEAK", drawstyle = "steps-mid", linewidth =1.5)
     # plot(t_steps, P_r["PV"][:], label="BID_PV", color = "tab:purple", drawstyle = "steps-mid", linewidth =1.5)
     # plot(t_steps, P_r["wind"][:], label="BID_WIND", color = "tab:red", drawstyle = "steps-mid", linewidth =1.5)
     # title("QUANTITY BIDS")
     # grid()
end


flag_statistics = true
if flag_statistics == true
    using Statistics
    avg_mp = mean(lambda)
    max_mp = maximum(lambda)
    min_mp = minimum(lambda)
    var_mp = var(lambda)

    #CALCULATING SOCIAL WELFARE
    SW =
        (sum(
        d[t] .* p_cap
        for t in t_steps)
    +
        sum( # minus the producers surplus
            # -
            # dch[t] * p_dch_bid_strat[t]
            # +
            # ch[t] * p_ch_bid_strat[t]
            -
            res[t] * dat_str.WIND_STR.MC_g[t]
            -
            res_PV[t] * dat_nonstr.PV.c_g[t]
            -
            gen_1[t] * dat_nonstr[:GEN1].c_g[t]
            -
            gen_2[t] * dat_nonstr[:GEN2].c_g[t]
            -
            gen_3[t] * dat_nonstr[:GEN3].c_g[t]
            -
            gen_5[t] * dat_nonstr[:GEN5].c_g[t]
            -
            gen_6[t] * dat_nonstr[:GEN6].c_g[t]
            -
            gen_7[t] * dat_nonstr[:GEN7].c_g[t]
            -
            gen_8[t] * dat_nonstr[:GEN8].c_g[t]
            -
            gen[t] * dat_str.GEN_STR.MC_g[t]
        for t in t_steps)
    )

    println("STATISTICS:
    ....Avg market price: $avg_mp
    ....Maximum: $max_mp
    ....Minimum: $min_mp
    ....Social Welfare: $SW")
end



push!(results_validation, ("SCR-NCP-$current $bidding", profit_ESS, profit_GEN, profit_RES, SW, avg_mp, 1.0, 1.0, optimizer))

CSV.write("C:/Users/u0125813/Desktop/Leuven/KU_Leuven/EPECs/scr_nlp/RESULTS/CASE2/results_validation.csv", results_validation)
results_validation = CSV.read("C:/Users/u0125813/Desktop/Leuven/KU_Leuven/EPECs/scr_nlp/RESULTS/CASE2/results_validation.csv")

results_validation = DataFrame(model= String[], ESS = Float64[], GEN = Float64[], WIND_STR = Float64[],
                            SW= Float64[], Avg_MP = Float64[], ϵ1 = Float64[], ϵ2 = Float64[], solver = [])



# insert!.(eachcol(results_validation, false), 5, ["SCR-NCP-$current $bidding", profit_ESS, profit_GEN, profit_RES, SW, avg_mp, 1.0, 1.0, "IPopt"])
using ExcelFiles
results_validation |> save("C:/Users/u0125813/Desktop/Leuven/KU_Leuven/EPECs/scr_nlp/RESULTS/CASE2/results_validation.xlsx")
