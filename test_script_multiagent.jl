using JuMP
using PyPlot
using DataFrames
using StructArrays
using CSV
using LinearAlgebra
using KNITRO
using Gurobi
using Statistics

#####TESTS FOR SCR-NCP QUANTITY, TESTED VIA DIAGONALIZATION COMPLEX#######
include("EPECs/EPEC_quantity_mod.jl")
include("diagonalization.jl")
include("MCPs/MCP_complex.jl")

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


#RAMPING LIMIT OF THE STRATEGIC GENERATOR:
Rᴳ = 1.0
res_dat = CSV.read("wind_PV.csv", header=true, delim= ";")
eight_res  = res_dat[1:3:end,:] #selecting every 3rd of the res_dat
renewable_nonstr = ["PV"]
renewable_str = ["WIND"]
capacity_re_nonstr = eight_res.PV/100
capacity_re_str =  (eight_res.Wind/100) * 1.5

#LOADING THE DEMAND
load_dat = CSV.read("load.csv",header=true, delim= ";")
load = load_dat.load[1:3:end,:]/100


eff_ch = eff_dch = generators_dat.ESS[end]
energy_max = (generators_dat.ESS[end-2]/100)*1.5
power_max = (generators_dat.ESS[end-4]/100)*1.5

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


function compute_profits(EPEC)
    lambda = value.(EPEC[:λ][:])
    #NON-STRATEGIC GENERATORS
    gen_1 = value.(EPEC[:g]["GEN1",:])
    gen_2 = value.(EPEC[:g]["GEN2",:])
    gen_3 = value.(EPEC[:g]["GEN3",:])
    gen_5 = value.(EPEC[:g]["GEN5",:])
    gen_6 = value.(EPEC[:g]["GEN6",:])
    gen_7 = value.(EPEC[:g]["GEN7",:])
    gen_8 = value.(EPEC[:g]["GEN8",:])

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
    gen_1 = value.(MCP[:g]["GEN1",:])
    gen_2 = value.(MCP[:g]["GEN2",:])
    gen_3 = value.(MCP[:g]["GEN3",:])
    gen_5 = value.(MCP[:g]["GEN5",:])
    gen_6 = value.(MCP[:g]["GEN6",:])
    gen_7 = value.(MCP[:g]["GEN7",:])
    gen_8 = value.(MCP[:g]["GEN8",:])

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
t1 = time()
MCP = build_MCP(t_steps, non_str_gens, non_str_res, optimizer)
optimize!(MCP)
objective_value(MCP)
t2 = time()
t_MCP = t2 - t1

profit_ESS, profit_GEN, profit_RES = compute_profits_MCP()
flag_statistics = true
if flag_statistics == true
    lambda = value.(MCP[:λ][:])
    #NON-STRATEGIC GENERATORS
    gen_1 = value.(MCP[:g]["GEN1",:])
    gen_2 = value.(MCP[:g]["GEN2",:])
    gen_3 = value.(MCP[:g]["GEN3",:])
    gen_5 = value.(MCP[:g]["GEN5",:])
    gen_6 = value.(MCP[:g]["GEN6",:])
    gen_7 = value.(MCP[:g]["GEN7",:])
    gen_8 = value.(MCP[:g]["GEN8",:])

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
    SW_MCP = SW =
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
            gen_1[t] * dat_nonstr["GEN1"].c_g[t]
            -
            gen_2[t] * dat_nonstr["GEN2"].c_g[t]
            -
            gen_3[t] * dat_nonstr["GEN3"].c_g[t]
            -
            gen_5[t] * dat_nonstr["GEN5"].c_g[t]
            -
            gen_6[t] * dat_nonstr["GEN6"].c_g[t]
            -
            gen_7[t] * dat_nonstr["GEN7"].c_g[t]
            -
            gen_8[t] * dat_nonstr["GEN8"].c_g[t]
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

simu_EPEC_quantity = DataFrame(time_build = [], time = [], bidding = String[], SC_obj = String[], ESS = Float64[], GEN = Float64[], WIND_STR = Float64[],
                            SW= Float64[], PS= Float64[], Avg_MP = Float64[], ϵ1 = Float64[], ϵ2 = Float64[],
                            confirmedNE = Bool[], solver = [], ESSᴰ = Float64[], GENᴰ = Float64[], WIND_STRᴰ = Float64[])

PS_MCP = profit_ESS + profit_GEN + profit_RES

push!(simu_EPEC_quantity, (0, t_MCP, "MCP", "NA", profit_ESS, profit_GEN, profit_RES,
        SW, profit_ESS + profit_GEN + profit_RES, avg_mp, 0.0, 0.0, true, "GUROBI", 0, 0, 0))


# obj_types = ["empty", "competitive", "collusive", "ESS_fav"]
obj_types = ["empty"]
param_adj = [i for i in range(1.0, stop = 2.0, step = 0.5)]
param_adj = [1,10,20,30,40,50,60,70,80,90,100]
flag_sequanctial_red = true
flag_statistics = true


for o in obj_types
    for n in param_adj

        global feasibility = true

        optimizer = with_optimizer(KNITRO.Optimizer,
                        # algorithm = 1,
                        hessian_no_f=1,
                        datacheck=0,
                        outlev=0,
                        # ms_enable= 1,
                        # ms_maxsolves = 50
                        )

        SP_obj = o
        bidding = "quantity"
        println("....Building EPEC model with objective: $o, and epsilon: $n ...")
        t_1_build = time()
        EPEC = build_EPEC_quantity_mod(t_steps, non_str_gens, non_str_res, str_agents, optimizer, SP_obj)
        t_2_build = time()
        t_build_EPEC = t_2_build - t_1_build
        println("....Solving EPEC model with objective: $o, and epsilon: $n ...")
        if flag_sequanctial_red == true
            global iter_limit = 40
            global iterations = 0
            global tol = 0.0001
            global epsilon1_new = n
            global epsilon2_new = n

            JuMP.fix(EPEC[:ϵ1], epsilon1_new; force = true)
            JuMP.fix(EPEC[:ϵ2], epsilon2_new; force = true)
            optimize!(EPEC)
            objective_value(EPEC)
            profit_ESS, profit_GEN, profit_RES = compute_profits(EPEC)
            println("PROFITS:
                            PROFIT_ESS: $profit_ESS,
                            PROFIT_GEN: $profit_GEN,
                            PROFIT_RES: $profit_RES")



            ############THE SMALL SCRIPT SNIPPET FOR UPDATING THE REGULARIZATION TERM
            t1 = time()
            while epsilon1_new >= tol || epsilon2_new >= tol
                #After reaching iteration limit exit:
                if iterations >= iter_limit
                    break
                end
                #########SETTING THE WARM START POINT TO THE PREVIOUS SOLUTION#########
                #if terminition status is local optimal, otherwise exit
                term_stat = termination_status(EPEC)
                if term_stat == MOI.OPTIMAL || term_stat == MOI.LOCALLY_SOLVED
                    set_start_value.(all_variables(EPEC), value.(all_variables(EPEC)))
                else
                    feasibility = false
                    break
                end
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
                profit_ESS, profit_GEN, profit_RES = compute_profits(EPEC)
                println("PROFITS:
                                PROFIT_ESS: $profit_ESS,
                                PROFIT_GEN: $profit_GEN,
                                PROFIT_RES: $profit_RES")
                ############PROGRESSING IN TIME##############
                global iterations += 1
                global prev_vals = value.(all_variables(EPEC))
            end
            t2 = time()
            t_EPEC = t2 - t1
        end

        if feasibility == true
            println("....Model was solved to a feasible point in $t_EPEC sec :) ...")
            lambda = value.(EPEC[:λ][:])
            #NON-STRATEGIC GENERATORS
            gen_1 = value.(EPEC[:g]["GEN1",:])
            gen_2 = value.(EPEC[:g]["GEN2",:])
            gen_3 = value.(EPEC[:g]["GEN3",:])
            gen_5 = value.(EPEC[:g]["GEN5",:])
            gen_6 = value.(EPEC[:g]["GEN6",:])
            gen_7 = value.(EPEC[:g]["GEN7",:])
            gen_8 = value.(EPEC[:g]["GEN8",:])

            #DEMAND
            d = value.(EPEC[:d])
            #NON-STRATEGIC RES
            res_PV = value.(EPEC[:w][:PV,:])
            #DISPATCHED QUANTITIES FROMS STR AGENTS
            dch = value.(EPEC[:dch_st])
            ch = value.(EPEC[:ch_st])
            gen = value.(EPEC[:g_st])
            res = value.(EPEC[:w_st])

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
                        gen_1[t] * dat_nonstr["GEN1"].c_g[t]
                        -
                        gen_2[t] * dat_nonstr["GEN2"].c_g[t]
                        -
                        gen_3[t] * dat_nonstr["GEN3"].c_g[t]
                        -
                        gen_5[t] * dat_nonstr["GEN5"].c_g[t]
                        -
                        gen_6[t] * dat_nonstr["GEN6"].c_g[t]
                        -
                        gen_7[t] * dat_nonstr["GEN7"].c_g[t]
                        -
                        gen_8[t] * dat_nonstr["GEN8"].c_g[t]
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

            global bidding = "complex"
            iter_diag, profit_ESS_MILP_ESS, profit_GEN_MILP_GEN, profit_RES_MILP_WIND = run_diagonaluization(EPEC)

            if iter_diag == 0
                confirmed_by_diag = true
            else
                confirmed_by_diag = false
            end

            push!(simu_EPEC_quantity, (t_build_EPEC, t_EPEC, "SCR-NC", "$SP_obj", profit_ESS, profit_GEN, profit_RES,
                    SW, profit_ESS + profit_GEN + profit_RES, avg_mp, n, n, confirmed_by_diag, "KNITRO",
                    profit_ESS_MILP_ESS, profit_GEN_MILP_GEN, profit_RES_MILP_WIND))

        else
            println("....Model converged to an infeasible point, after $iterations iterations :( ...")
        end

    end
end

deleterows!(simu_EPEC_quantity_10_M1, 5)

mean(simu_EPEC_quantity_10_M0.SW)
mean(simu_EPEC_quantity_10_M1.SW)
mean(simu_EPEC_quantity_10_M2.SW)
mean(simu_EPEC_quantity_10_M3.SW)

mean(simu_EPEC_quantity_10_M0.ESS)
mean(simu_EPEC_quantity_10_M1.ESS)
mean(simu_EPEC_quantity_10_M2.ESS)
mean(simu_EPEC_quantity_10_M3.ESS)

mean(simu_EPEC_quantity_10_M0.GEN)
mean(simu_EPEC_quantity_10_M1.GEN)
mean(simu_EPEC_quantity_10_M2.GEN)
mean(simu_EPEC_quantity_10_M3.GEN)

mean(simu_EPEC_quantity_10_M0.WIND_STR)
mean(simu_EPEC_quantity_10_M1.WIND_STR)
mean(simu_EPEC_quantity_10_M2.WIND_STR)
mean(simu_EPEC_quantity_10_M3.WIND_STR)

std(simu_EPEC_quantity_10_M0.ESS)
std(simu_EPEC_quantity_10_M1.ESS)
std(simu_EPEC_quantity_10_M2.ESS)
std(simu_EPEC_quantity_10_M3.ESS)

std(simu_EPEC_quantity_10_M0.SW)
std(simu_EPEC_quantity_10_M1.SW)
std(simu_EPEC_quantity_10_M2.SW)
std(simu_EPEC_quantity_10_M3.SW)

std(simu_EPEC_quantity_10_M0.GEN)
std(simu_EPEC_quantity_10_M1.GEN)
std(simu_EPEC_quantity_10_M2.GEN)
std(simu_EPEC_quantity_10_M3.GEN)

std(simu_EPEC_quantity_10_M0.WIND_STR)
std(simu_EPEC_quantity_10_M1.WIND_STR)
std(simu_EPEC_quantity_10_M2.WIND_STR)
std(simu_EPEC_quantity_10_M3.WIND_STR)



deleterows!(simu_EPEC_quantity_10_M0, 5)
simu_EPEC_quantity_10_M3 = simu_EPEC_quantity
simu_EPEC_quantity
simu_EPEC_quantity_ESS_fav_withtime = deepcopy(simu_EPEC_quantity)
# deleterows!(simu_EPEC_quantity_empty_withtime, 1)
Statistics.mean(simu_EPEC_quantity_ESS_fav_withtime.time[:])
simu_EPEC_quantity_ESS_fav_withtime |> @filter(_.confirmedNE == true) |> DataFrame

#
CSV.write("C:/Users/u0125813/Desktop/Leuven/KU_Leuven/EPECs/scr_nlp/RESULTS/CASE2/simu_EPEC_quantity_ESS_fav_withtime.csv", simu_EPEC_quantity_ESS_fav_withtime)

CSV.write("C:/Users/u0125813/Desktop/Leuven/KU_Leuven/EPECs/scr_nlp/RESULTS/CASE2/simu_EPEC_quantity_collusive_withtime.csv", simu_EPEC_quantity_collusive_withtime)

CSV.write("C:/Users/u0125813/Desktop/Leuven/KU_Leuven/EPECs/scr_nlp/RESULTS/CASE2/simu_EPEC_quantity_empty_withtime.csv", simu_EPEC_quantity_empty_withtime)

CSV.write("C:/Users/u0125813/Desktop/Leuven/KU_Leuven/EPECs/scr_nlp/RESULTS/CASE2/simu_EPEC_quantity_competitive_withtime.csv", simu_EPEC_quantity_competitive_withtime)


CSV.write("C:/Users/u0125813/Desktop/Leuven/KU_Leuven/EPECs/scr_nlp/RESULTS/CASE2/simu_EPEC_quantity_10_M0.csv", simu_EPEC_quantity_10_M0)

CSV.write("C:/Users/u0125813/Desktop/Leuven/KU_Leuven/EPECs/scr_nlp/RESULTS/CASE2/simu_EPEC_quantity_10_M1.csv", simu_EPEC_quantity_10_M1)

CSV.write("C:/Users/u0125813/Desktop/Leuven/KU_Leuven/EPECs/scr_nlp/RESULTS/CASE2/simu_EPEC_quantity_10_M2.csv", simu_EPEC_quantity_10_M2)

CSV.write("C:/Users/u0125813/Desktop/Leuven/KU_Leuven/EPECs/scr_nlp/RESULTS/CASE2/simu_EPEC_quantity_10_M3.csv", simu_EPEC_quantity_10_M3)


simu_EPEC_quantity_10_M1 = CSV.read("C:/Users/u0125813/Desktop/Leuven/KU_Leuven/EPECs/scr_nlp/RESULTS/CASE2/simu_EPEC_quantity_10_M1.csv"; copycols = true)


using PyCall, PyPlot
pygui(true)
figure()
#plot here the line with the optimal
plt.scatter(simu_EPEC_quantity_10_M0.ϵ1, (simu_EPEC_quantity_10_M2.ESS - simu_EPEC_quantity_10_M0.ESS), color = "tab:purple", marker="D")
plt.scatter(simu_EPEC_quantity_10_M0.ϵ1, (simu_EPEC_quantity_10_M1.ESS - simu_EPEC_quantity_10_M0.ESS), color = "tab:cyan", marker="s")
# plt.plot(results_MPPDC_ESS.iter[6:11], results_MPPDC_ESS.profit[6:11], label="ϵ_init = 10", color ="k", linestyle="--", linewidth =1.5, marker="D", markersize = 10)
# plt.plot(results_MPPDC_ESS.iter[12:18], results_MPPDC_ESS.profit[12:18], label="ϵ_init = 100", color ="k", linestyle=":", linewidth =1.5, marker="s", markersize = 10)
# plt.plot(results_MPPDC_ESS.iter[19:26], results_MPPDC_ESS.profit[19:26], label="ϵ_init = 1000", color ="k", linestyle="-.", linewidth =1.5, marker="p", markersize = 10)
fill_between(simu_EPEC_quantity_10_M0.ϵ1, (mean(simu_EPEC_quantity_10_M2.ESS) - mean(simu_EPEC_quantity_10_M0.ESS)), (simu_EPEC_quantity_10_M2.ESS - simu_EPEC_quantity_10_M0.ESS) , color="tab:purple", alpha=0.4, label = "∇₂",)
fill_between(simu_EPEC_quantity_10_M0.ϵ1, (mean(simu_EPEC_quantity_10_M1.ESS) - mean(simu_EPEC_quantity_10_M0.ESS)), (simu_EPEC_quantity_10_M1.ESS - simu_EPEC_quantity_10_M0.ESS) , color="tab:cyan", alpha=0.4, label = "∇₁",)
plt.plot(simu_EPEC_quantity_10_M0.ϵ1, ones(10).*(mean(simu_EPEC_quantity_10_M1.ESS) - mean(simu_EPEC_quantity_10_M0.ESS)), linestyle="--", label = "Avg. ∇₁", color = "k", linewidth = 2.5)
plt.plot(simu_EPEC_quantity_10_M0.ϵ1, ones(10).*(mean(simu_EPEC_quantity_10_M2.ESS) - mean(simu_EPEC_quantity_10_M0.ESS)), linestyle=":", label = "Avg. ∇₂", color = "k", linewidth = 2.5)
plt.xticks(simu_EPEC_quantity_10_M0.ϵ1)
plt.ylabel("∇ Profit ESS")
plt.xlabel("ϵ")
plt.legend()
plt.grid()
plt.savefig("errors",bbox_inches="tight",dpi=300)



figure()
#plot here the line with the optimal
plt.scatter(simu_EPEC_quantity_10_M0.ϵ1, (simu_EPEC_quantity_10_M2.ESS/10 - simu_EPEC_quantity_10_M0.ESS/10), color = "tab:cyan", marker="D")
plt.scatter(simu_EPEC_quantity_10_M0.ϵ1, (simu_EPEC_quantity_10_M1.ESS/10 - simu_EPEC_quantity_10_M0.ESS/10), color = "tab:purple", marker="s")
# plt.plot(results_MPPDC_ESS.iter[6:11], results_MPPDC_ESS.profit[6:11], label="ϵ_init = 10", color ="k", linestyle="--", linewidth =1.5, marker="D", markersize = 10)
# plt.plot(results_MPPDC_ESS.iter[12:18], results_MPPDC_ESS.profit[12:18], label="ϵ_init = 100", color ="k", linestyle=":", linewidth =1.5, marker="s", markersize = 10)
# plt.plot(results_MPPDC_ESS.iter[19:26], results_MPPDC_ESS.profit[19:26], label="ϵ_init = 1000", color ="k", linestyle="-.", linewidth =1.5, marker="p", markersize = 10)
fill_between(simu_EPEC_quantity_10_M0.ϵ1, (mean(simu_EPEC_quantity_10_M2.ESS/10) - mean(simu_EPEC_quantity_10_M0.ESS/10)), (simu_EPEC_quantity_10_M2.ESS/10 - simu_EPEC_quantity_10_M0.ESS/10) , color="tab:cyan", alpha=0.4, label = "ESS ∇₂",)
fill_between(simu_EPEC_quantity_10_M0.ϵ1, (mean(simu_EPEC_quantity_10_M1.ESS/10) - mean(simu_EPEC_quantity_10_M0.ESS/10)), (simu_EPEC_quantity_10_M1.ESS/10 - simu_EPEC_quantity_10_M0.ESS/10) , color="tab:purple", alpha=0.4, label = "ESS ∇₁",)
plt.plot(simu_EPEC_quantity_10_M0.ϵ1, ones(10).*(mean(simu_EPEC_quantity_10_M1.ESS/10) - mean(simu_EPEC_quantity_10_M0.ESS/10)), linestyle="--", label = "ESS Avg. ∇₁", color = "tab:purple", linewidth = 2.5)
plt.plot(simu_EPEC_quantity_10_M0.ϵ1, ones(10).*(mean(simu_EPEC_quantity_10_M2.ESS/10) - mean(simu_EPEC_quantity_10_M0.ESS/10)), linestyle=":", label = "ESS Avg. ∇₂", color = "tab:cyan", linewidth = 2.5)
plt.xticks(simu_EPEC_quantity_10_M0.ϵ1)


plt.scatter(simu_EPEC_quantity_10_M0.ϵ1, (simu_EPEC_quantity_10_M3.WIND_STR/10 - simu_EPEC_quantity_10_M0.WIND_STR/10), color = "tab:red", marker="p")
plt.scatter(simu_EPEC_quantity_10_M0.ϵ1, (simu_EPEC_quantity_10_M1.WIND_STR/10 - simu_EPEC_quantity_10_M0.WIND_STR/10), color = "tab:olive", marker="^")
# plt.plot(results_MPPDC_ESS.iter[6:11], results_MPPDC_ESS.profit[6:11], label="ϵ_init = 10", color ="k", linestyle="--", linewidth =1.5, marker="D", markersize = 10)
# plt.plot(results_MPPDC_ESS.iter[12:18], results_MPPDC_ESS.profit[12:18], label="ϵ_init = 100", color ="k", linestyle=":", linewidth =1.5, marker="s", markersize = 10)
# plt.plot(results_MPPDC_ESS.iter[19:26], results_MPPDC_ESS.profit[19:26], label="ϵ_init = 1000", color ="k", linestyle="-.", linewidth =1.5, marker="p", markersize = 10)
fill_between(simu_EPEC_quantity_10_M0.ϵ1, (mean(simu_EPEC_quantity_10_M3.WIND_STR/10) - mean(simu_EPEC_quantity_10_M0.WIND_STR/10)), (simu_EPEC_quantity_10_M3.WIND_STR/10 - simu_EPEC_quantity_10_M0.WIND_STR/10) , color="tab:red", alpha=0.4, label = "RG ∇₃",)
fill_between(simu_EPEC_quantity_10_M0.ϵ1, (mean(simu_EPEC_quantity_10_M1.WIND_STR/10) - mean(simu_EPEC_quantity_10_M0.WIND_STR/10)), (simu_EPEC_quantity_10_M1.WIND_STR/10 - simu_EPEC_quantity_10_M0.WIND_STR/10) , color="tab:olive", alpha=0.4, label = "RG ∇₁",)
plt.plot(simu_EPEC_quantity_10_M0.ϵ1, ones(10).*(mean(simu_EPEC_quantity_10_M1.WIND_STR/10) - mean(simu_EPEC_quantity_10_M0.WIND_STR/10)), linestyle="--", label = "RG Avg. ∇₁", color = "tab:olive", linewidth = 2.5)
plt.plot(simu_EPEC_quantity_10_M0.ϵ1, ones(10).*(mean(simu_EPEC_quantity_10_M3.WIND_STR/10) - mean(simu_EPEC_quantity_10_M0.WIND_STR/10)), linestyle=":", label = "RG Avg. ∇₃", color = "tab:red", linewidth = 2.5)
plt.xticks(simu_EPEC_quantity_10_M0.ϵ1)
plt.ylabel("∇ Profit")
plt.xlabel("ϵ")
plt.legend(ncol=3)
plt.grid()
plt.savefig("errors",bbox_inches="tight",dpi=300)
