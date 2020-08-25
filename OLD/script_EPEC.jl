using JuMP
using Ipopt
using PyPlot
using DataFrames
using StructArrays
using CSV
using LinearAlgebra
using KNITRO
using Gurobi
using Statistics


include("EPECs/EPEC_complex.jl") #EPEC with complex bids
include("EPECs/EPEC_PC_complex.jl") #Price consistent linear formulation, MILP based
include("EPECs/EPEC_PC_MPPDC_complex.jl") #Price consistent non-linear formulation, MPPDC
include("EPECs/EPEC_PC_MPPDC_quantity.jl") #Price consistent non-linear formulation, MPPDC
include("EPECs/EPEC_quantity_mod.jl")
include("EPECs/EPEC_quantity_mod_PC.jl")
include("MCPs/MCP_complex.jl") #MCP

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


function compute_profits()
    lambda = value.(EPEC[:λ][:])
    #NON-STRATEGIC GENERATORS
    gen_1 = value.(EPEC[:g][:GEN1,:])
    gen_2 = value.(EPEC[:g][:GEN2,:])
    gen_3 = value.(EPEC[:g][:GEN3,:])
    gen_5 = value.(EPEC[:g][:GEN5,:])
    gen_6 = value.(EPEC[:g][:GEN6,:])
    gen_7 = value.(EPEC[:g][:GEN7,:])
    gen_8 = value.(EPEC[:g][:GEN8,:])

    #DEMAND
    d = value.(EPEC[:d])
    #NON-STRATEGIC RES
    res_PV = value.(EPEC[:w][:PV,:])
    #DISPATCHED QUANTITIES FROMS STR AGENTS
    dch = value.(EPEC[:dch_st])
    ch = value.(EPEC[:ch_st])
    gen = value.(EPEC[:g_st])
    res = value.(EPEC[:w_st])
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

function compute_profits_MCP()
    lambda = value.(MCP[:λ][:])
    #NON-STRATEGIC GENERATORS
    gen_1 = value.(MCP[:g][:GEN1,:])
    gen_2 = value.(MCP[:g][:GEN2,:])
    gen_3 = value.(MCP[:g][:GEN3,:])
    gen_5 = value.(MCP[:g][:GEN5,:])
    gen_6 = value.(MCP[:g][:GEN6,:])
    gen_7 = value.(MCP[:g][:GEN7,:])
    gen_8 = value.(MCP[:g][:GEN8,:])

    #DEMAND
    d = value.(MCP[:d])
    #NON-STRATEGIC RES
    res_PV = value.(MCP[:w][:PV,:])
    #DISPATCHED QUANTITIES FROMS STR AGENTS
    dch = value.(MCP[:dch_st])
    ch = value.(MCP[:ch_st])
    gen = value.(MCP[:g_st])
    res = value.(MCP[:w_st])
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

###########FIRST SOLVE THE MCP EPEC############
optimizer = with_optimizer(Gurobi.Optimizer)
MCP = build_MCP(t_steps, non_str_gens, non_str_res, optimizer)
optimize!(MCP)
objective_value(MCP)

profit_ESS, profit_GEN, profit_RES = compute_profits_MCP()
flag_statistics = true
if flag_statistics == true
    lambda = value.(MCP[:λ][:])
    #NON-STRATEGIC GENERATORS
    gen_1 = value.(MCP[:g][:GEN1,:])
    gen_2 = value.(MCP[:g][:GEN2,:])
    gen_3 = value.(MCP[:g][:GEN3,:])
    gen_5 = value.(MCP[:g][:GEN5,:])
    gen_6 = value.(MCP[:g][:GEN6,:])
    gen_7 = value.(MCP[:g][:GEN7,:])
    gen_8 = value.(MCP[:g][:GEN8,:])

    #DEMAND
    d = value.(MCP[:d])
    #NON-STRATEGIC RES
    res_PV = value.(MCP[:w][:PV,:])
    #DISPATCHED QUANTITIES FROMS STR AGENTS
    dch = value.(MCP[:dch_st])
    ch = value.(MCP[:ch_st])
    gen = value.(MCP[:g_st])
    res = value.(MCP[:w_st])

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

# ###########SECOND SOLVE THE PRICE CONSISTENT EPEC############
# optimizer = with_optimizer(Gurobi.Optimizer)
# optimizer = (with_optimizer(Ipopt.Optimizer,
#                                 max_iter = 3000,
#                                 print_level = 0,
#                                 # linear_solver = "ma86"
#                                 # tol = 1.0e-12,
#                                 # acceptable_tol = 1.0e-9
#                                 ))

optimizer = with_optimizer(KNITRO.Optimizer,
                algorithm =3,
                hessian_no_f=1,
                datacheck=0,
                outlev=0,
                # ms_enable= 1,
                # ms_maxsolves = 5
                )

#Social Planner objective type####
flag_sequanctial_red = true
EPEC = build_EPEC_PC_MPPDC_quantity()
#
if flag_sequanctial_red == true
    global iter_limit = 100
    global iterations = 0
    global tol = 1e-4
    global epsilon_new = 100
    JuMP.fix(EPEC[:ϵ], epsilon_new; force = true)

    optimize!(EPEC)
    objective_value(EPEC)
    ############THE SMAL SCRIPT SNIPPET FOR UPDATING THE REGULARIZATION TERM
    while epsilon_new >= tol && iterations < iter_limit
        #########SETTING THE WARM START POINT TO THE PREVIOUS SOLUTION#########
        #if terminition status is local optimal, otherwise should not use the solution as starting point
        term_stat = termination_status(EPEC)
        # if term_stat == MOI.OPTIMAL || term_stat == MOI.LOCALLY_SOLVED
        set_start_value.(all_variables(EPEC), value.(all_variables(EPEC)))
        ###################FIXING EPRSILONS TO THE REDUCED PARAMS##############
        global epsilon_new = value.(EPEC[:ϵ])*0.5
        # else
        #     set_start_value.(all_variables(EPEC), prev_vals)
        #     global epsilon_new = value.(EPEC[:ϵ])*0.99
        # end
        fix(EPEC[:ϵ], epsilon_new; force = true)
        println("\\
        ######### \\
        THE NEW RELAXATION TERMS: \\
                    ϵ: $epsilon_new and \\
                    #ITERATIONS: $iterations \\
                    #TERMINATION: $term_stat
        #########")
        ################OPTIMIZING############
        optimize!(EPEC)
        profit_ESS, profit_GEN, profit_RES = compute_profits()
        println("PROFITS:
                        PROFIT_ESS: $profit_ESS,
                        PROFIT_GEN: $profit_GEN,
                        PROFIT_RES: $profit_RES")

        ##############CACLCULATING PROFIT############
        # profit = compute_profit()
        # println("PROFIT: ", profit)
        ############PROGRESSING IN TIME##############
        global iterations += 1
        global prev_vals = value.(all_variables(EPEC))
    end
end

# EPEC = build_EPEC_PC_complex()
# JuMP.fix(EPEC[:M1], 10000; force = true)
# JuMP.fix(EPEC[:M1], 10000; force = true)
# optimize!(EPEC)
# objective_value(EPEC)
#
profit_ESS, profit_GEN, profit_RES = compute_profits()

#########SOLVING SCR-NLPs##########
optimizer = with_optimizer(KNITRO.Optimizer,
                # algorithm =2,
                hessian_no_f=1,
                datacheck=0,
                outlev=0,
                # ms_enable= 1,
                # ms_maxsolves = 50
                )
#
# optimizer = (with_optimizer(Ipopt.Optimizer,
#                                 max_iter = 3000,
#                                 print_level = 0
#                                 # tol = 1.0e-16,
#                                 # acceptable_tol = 1.0e-10
#                                 ))
##Social Planner objective type####

# mod_dat_str = deepcopy(dat_str)
# mod_dat_str.GEN_STR.G[:] = value.(EPEC[:G_ST])
# mod_dat_str.WIND_STR.G[:] = value.(EPEC[:W_ST])
# mod_dat_str.ESS_STR.DCH[:] = value.(EPEC[:DCH_ST])
# mod_dat_str.ESS_STR.CH[:] = value.(EPEC[:CH_ST])
#
SP_obj = "empty"
bidding = "quantity"

if bidding ==  "complex"
    EPEC = build_EPEC_complex(t_steps, non_str_gens, non_str_res, str_agents, optimizer, SP_obj)
elseif bidding ==  "quantity"
    EPEC = build_EPEC_quantity_mod(t_steps, non_str_gens, non_str_res, str_agents, optimizer, SP_obj)
else
    println("ERROR: Pick 'complex' or 'quantity' for bidding.")
end

#For the price consistent SCR-NCP:
# EPEC = build_EPEC_quantity_mod_PC(t_steps, non_str_gens, non_str_res, str_agents, optimizer, SP_obj)

##########TO FIX THE EPEC RESULTS#########

# for t in t_steps
#     fix(EPEC[:G_ST][t], mod_dat_str.GEN_STR.G[t]; force = true)
#     fix(EPEC[:W_ST][t], mod_dat_str.WIND_STR.G[t]; force = true)
#     fix(EPEC[:DCH_ST][t], mod_dat_str.ESS_STR.DCH[t]; force = true)
#     fix(EPEC[:CH_ST][t], mod_dat_str.ESS_STR.CH[t]; force = true)
# end
#######OR WARMSTART#########
flag_sequanctial_red = true
#START TIMER
@time begin

if flag_sequanctial_red == true
    global iter_limit = 40
    global iterations = 0
    global tol = 0.01
    global epsilon1_new = 9.0
    global epsilon2_new = 9.0

    JuMP.fix(EPEC[:ϵ1], epsilon1_new; force = true)
    JuMP.fix(EPEC[:ϵ2], epsilon2_new; force = true)
    optimize!(EPEC)
    objective_value(EPEC)
    profit_ESS, profit_GEN, profit_RES = compute_profits()
    println("PROFITS:
                    PROFIT_ESS: $profit_ESS,
                    PROFIT_GEN: $profit_GEN,
                    PROFIT_RES: $profit_RES")



    ############THE SMALL SCRIPT SNIPPET FOR UPDATING THE REGULARIZATION TERM
    while epsilon1_new >= tol || epsilon2_new >= tol
        #After reaching iteration limit exit:
        if iterations >= iter_limit
            break
        end
        #########SETTING THE WARM START POINT TO THE PREVIOUS SOLUTION#########
        #if terminition status is local optimal, otherwise should not use the solution as starting point
        term_stat = termination_status(EPEC)
        # if term_stat == MOI.OPTIMAL || term_stat == MOI.LOCALLY_SOLVED
        set_start_value.(all_variables(EPEC), value.(all_variables(EPEC)))
        ###################FIXING EPRSILONS TO THE REDUCED PARAMS##############
        global epsilon1_new = value.(EPEC[:ϵ1])*0.5
        if epsilon2_new >= tol
            global epsilon2_new = value.(EPEC[:ϵ2])*0.5
        end
        # else
        #     set_start_value.(all_variables(EPEC), prev_vals)
        #     global epsilon_new = value.(EPEC[:ϵ])*0.99
        # end
        fix(EPEC[:ϵ1], epsilon1_new; force = true)
        fix(EPEC[:ϵ2], epsilon2_new; force = true)

        println("
        #########
        THE NEW RELAXATION TERMS:
                    # ϵ1: $epsilon1_new and
                    ϵ2: $epsilon2_new and
                    #ITERATIONS: $iterations
                    #TERMINATION: $term_stat
        #########")
        ################OPTIMIZING############
        optimize!(EPEC)
        ##############CACLCULATING PROFIT############
        profit_ESS, profit_GEN, profit_RES = compute_profits()
        println("PROFITS:
                        PROFIT_ESS: $profit_ESS,
                        PROFIT_GEN: $profit_GEN,
                        PROFIT_RES: $profit_RES")
        ############PROGRESSING IN TIME##############
        global iterations += 1
        global prev_vals = value.(all_variables(EPEC))
    end
end

# #END TIMER
end

lambda = value.(EPEC[:λ][:])
#NON-STRATEGIC GENERATORS
gen_1 = value.(EPEC[:g][:GEN1,:])
gen_2 = value.(EPEC[:g][:GEN2,:])
gen_3 = value.(EPEC[:g][:GEN3,:])
gen_5 = value.(EPEC[:g][:GEN5,:])
gen_6 = value.(EPEC[:g][:GEN6,:])
gen_7 = value.(EPEC[:g][:GEN7,:])
gen_8 = value.(EPEC[:g][:GEN8,:])

#DEMAND
d = value.(EPEC[:d])
#NON-STRATEGIC RES
res_PV = value.(EPEC[:w][:PV,:])
#DISPATCHED QUANTITIES FROMS STR AGENTS
dch = value.(EPEC[:dch_st])
ch = value.(EPEC[:ch_st])
gen = value.(EPEC[:g_st])
res = value.(EPEC[:w_st])

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


#############################################################

##########STORING RESULTS#########
#Define the df if was not existing
if @isdefined(results_EPEC) == false
    results_EPEC_quantity = DataFrame(bidding = String[], SC_obj = String[], ESS = Float64[], GEN = Float64[], WIND_STR = Float64[],
                                SW= Float64[], Avg_MP = Float64[], ϵ1 = Float64[], ϵ2 = Float64[], confirmedNE = Bool[], solver = [])
    # OR
    results_EPEC_OLD = CSV.read("C:/Users/u0125813/Desktop/Leuven/KU_Leuven/EPECs/scr_nlp/RESULTS/CASE2/results_EPEC_quantity.csv"; copycols = true)
end


profit_ESS, profit_GEN, profit_RES = compute_profits()

push!(results_EPEC_quantity, ("PC-MOPEC", "NA", profit_ESS, profit_GEN, profit_RES, SW, avg_mp, 1.0, 0.0, false, "KNITRO", profit_ESS + profit_GEN + profit_RES))
# results_EPEC_complex = vcat(results_EPEC_OLD, results_EPEC_complex)
results_EPEC_quantity = sort(results_EPEC_quantity, order(:SC_obj))
#SAVING THE RESULTS
CSV.write("C:/Users/u0125813/Desktop/Leuven/KU_Leuven/EPECs/scr_nlp/RESULTS/CASE2/results_EPEC_quantity.csv", results_EPEC_quantity)

using ExcelFiles
results_EPEC_complex |> save("C:/Users/u0125813/Desktop/Leuven/KU_Leuven/EPECs/scr_nlp/RESULTS/CASE2/output_complex.xlsx")

insert!.(eachcol(results_EPEC_complex, false), 13, ["SCR-NCP", "$SP_obj", profit_ESS, profit_GEN, profit_RES, SW, avg_mp, 10.0, 10.0, false, "KNITRO_ALG:1", profit_ESS + profit_GEN + profit_RES])

results_EPEC_quantity[:prod_surp] = [sum(results_EPEC_quantity[i, 3:5]) for i in 1:9]
round(results_EPEC_complex[:, 3:7], 2)


deleterows!(results_EPEC_quantity, 5)

results_EPEC_quantity.bidding[end] = "SCR-NCP"
