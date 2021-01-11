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
include("MILPs/MILP_quantity.jl") #The MILP model with quantity bidding
include("MPPDCs/NLPs/NLP_dual_MPPDC_quantity.jl") #The NLP model with quantity bidding


###############DECLARE PARTICIPATING AGENTS###########
rel_path = "DATA"
filepath = joinpath(@__DIR__, rel_path)
cd(filepath)
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


#SETTING THE CURRENT AGENT
#NOTE: ACCEPTED STRINGS: "ESS_STR", "GEN_STR", "WIND_STR"
current = "ESS_STR"
bidding = "quantity"
model_type = "MPEC"
model_type = "MPPDC"
model_type = "MILP"
objective_type = "dual" #Dual refers to the linearized objective functon, whereas primal is the original including bi-linear terms
flag_regularize = true

MPPDC = build_MPPDC_quantity(t_steps, non_str_gens,
non_str_res, current, optimizer, objective_type)

println("SOLUTION METHOD:
                ....Agent: $current
                ....Objective_type: $objective_type")


flag_sequanctial_red = true
#NOTE: MIND THE STEP SIZE
results_MPPDC_ESS = DataFrame(ϵ = [], iter = [], profit = [])
results_MPPDC_GEN = DataFrame(ϵ = [], iter = [], profit = [])
results_MPPDC_WIND = DataFrame(ϵ = [], iter = [], profit = [])

#FOR MPPDC THE SEQUANTIAL ALGORITHM
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


@show res_MPPDC = results_MPPDC_ESS.profit[end]

#FOR NLP
NLP = build_NLP_dual_quantity_mod(t_steps, non_str_gens,
non_str_res, str_agents, current, optimizer)

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

@show res_NLP = results_NLP_ESS.profit[end]

#OR BUILD WITH KKT CONDITIONS INSTEAD OF MPPDC, NOTE: With the linearized objective, LP solver may be used
optimizer = with_optimizer(Gurobi.Optimizer)
MILP = build_MILP_quantity(t_steps, non_str_gens,
non_str_res, current, optimizer, objective_type)
optimize!(MILP)

@show res_MILP = objective_value(MILP)
