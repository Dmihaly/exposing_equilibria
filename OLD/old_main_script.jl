using JuMP
using Ipopt
using PyPlot
using DataFrames
using StructArrays
using CSV
using LinearAlgebra
using KNITRO


###############DECLARE PARTICIPATING AGENTS###########
cd("C:/Users/u0125813/Desktop/Leuven/KU_Leuven/JULIA/Projects/SCR_NLP_PAPER/DATA/")

generators_dat = CSV.read("generators_data.csv", header=true, delim= ";")
MC_g = generators_dat[1,2:8]
MC_g_str = generators_dat[1,9]
#CAPACITIES OF GENERATORS
capacity_g = generators_dat[2,2:9]./100
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
load = load_dat.load[1:3:end,:]/150


eff_ch = eff_dch = generators_dat.ESS[end]
energy_max = generators_dat.ESS[end-2]/100
power_max = generators_dat.ESS[end-4]/100

MC_re_str = MC_re_nonstr = 0 #Zero marginal cost

str_stor = ["ESS_STR"]
str_gens = ["GEN_STR", "WIND_STR"]

generators = names(MC_g)
renewables = [:PV]

non_str_gens = vcat(generators)
non_str_res = vcat(renewables)
non_str_agents = vcat(generators, renewables)

#NOTE:ATM There is no non-str ESS involved in the pool
#MARKET PRICE LIMITS
p_cap = 100
p_floor = 0

#############INIT TIME STEPS##########
function t_interval(lenght, resolution)
    t_steps = [i for i in range(1, stop = lenght, step = resolution)]
    return t_steps
end

t_steps = t_interval(8,1)

###########INIT BID, TECH LIMIT DFs#########
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
    MC_dch = zeros(length(t_steps)) * p_floor,
    c_dch = zeros(length(t_steps)),
    CH = zeros(length(t_steps)),
    c_ch = zeros(length(t_steps)),
    MC_ch = ones(length(t_steps)) * p_cap)
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

#######THE ABSTRACT MODEL TYPE FOR THE MPPDC############
function init_MPPDC()
    ############INIT THE MODEL INSTANCE############
    MPPDC =  Model()

    #NOTE:OTHER STRATEGIC DEC VARS ARE INITIALIZED BUT FIXED
    @variables(MPPDC, begin
    ############UL DECISION VARS################
    #########PRICE BIDS#############
    c_g_st[t in t_steps] #Generation price bid str
    c_w_st[t in t_steps] #Wind price bid str
    c_dch_st[t in t_steps] #ESS discharge price bid str
    c_ch_st[t in t_steps] #ESS charge price bid str
    c_g[i in non_str_gens, t in t_steps] #Generation price bid nonstr
    c_w[k in non_str_res, t in t_steps] #Wind price bid nonstr
    c_d[t in t_steps] #WTP of the load
    #########QUANTITY BIDS############
    0 <= G_ST[t in t_steps] #Generation quantity bid
    0 <= W_ST[t in t_steps] #Wind quantity bid
    0 <= DCH_ST[t in t_steps] #ESS discharge quantity bid
    0 <= CH_ST[t in t_steps] #ESS charge quantity bid
    0 <= SOC[t in t_steps] #ESS State of Charge var
    E_max[t in t_steps] == 0 #max capacity of the ESS
    #############LL DECISION VARS###########
    #####PARAMETERS#########
    #NOTE: PARAMETERS NEED TO BE FIXED FROM OUTSIDE
    0 <= G[i in non_str_gens, t in t_steps] #Generation capacity of non-str generator
    0 <= W[k in non_str_res, t in t_steps] #Wind capacity of non-str re producer
    0 <= D[t in t_steps] #Max load
    ######PRIMALS#######
    0 <= g_st[t in t_steps] #Dispatched generation of the strategic agent
    0 <= w_st[t in t_steps] #Dispatched wind of the strategic agent
    0 <= g[i in non_str_gens, t in t_steps] #Dispatched generation
    0 <= w[k in non_str_res, t in t_steps] #Dispatched wind
    0 <= dch_st[t in t_steps] #Dispatched strategic discharge
    0 <= ch_st[t in t_steps] #Dispatched strategic charge
    0 <= d[t in t_steps] #Dispatched demand
    #######DUALS#########
    0 <= β⁺g[i in non_str_gens, t in t_steps] #Dual of the upper bound of the generation
    0 <= β⁻g[i in non_str_gens, t in t_steps] #Dual of the lower bound of the generation
    0 <= β⁺w[k in non_str_res, t in t_steps] #Dual of the upper bound of the RE generation
    0 <= β⁻w[k in non_str_res, t in t_steps] #Dual of the lower bound of the RE generation
    0 <= β⁺d[t in t_steps] #Dual of the upper bound of the demand
    0 <= β⁻d[t in t_steps] #Dual of the lower bound of the demand

    0 <= β⁺g_st[t in t_steps] #Dual of the upper bound of the str generation
    0 <= β⁻g_st[t in t_steps] #Dual of the lower bound of the str generation
    0 <= β⁺w_st[t in t_steps] #Dual of the upper bound of the str RE generation
    0 <= β⁻w_st[t in t_steps] #Dual of the lower bound of the str RE generation
    0 <= β⁺ch_st[t in t_steps] #Dual of the upper bound of the str charging
    0 <= β⁻ch_st[t in t_steps] #Dual of the lower bound of the str charging
    0 <= β⁺dch_st[t in t_steps] #Dual of the upper bound of the str discharging
    0 <= β⁻dch_st[t in t_steps] #Dual of the lower bound of the str discharging
    λ[t in t_steps] #Market price, dual of the energy balance cstr of LL
    ########CONSTANT TERMS FOR RAMPING###########
    ramp_limit[t in t_steps] == 1.0 #Set to 100% as default
    ηch == 1.0 #Set to 100% as default
    ηdch == 1.0 #Set to 100% as default
    ϵ == 0.0001
    end)

    ##########ENERGY BALANCE##########
    balance_eq = @expression(MPPDC, [t in t_steps],
    - g_st[t] - w_st[t] -dch_st[t] + ch_st[t]
    - sum(g[i,t] for i in non_str_gens) - sum(w[k,t] for k in non_str_res) + d[t])

    #######STRONG DUALITY EQ##########
    #PRIMAL OBJ - DUAL OBJ
    SDE_eq = @expression(MPPDC,
                        sum(
                        g_st[t] * c_g_st[t] +
                        w_st[t] * c_w_st[t] +
                        dch_st[t] * c_dch_st[t] +
                        - ch_st[t] * c_ch_st[t] +
                        sum(w[k,t] * c_w[k,t] for k in non_str_res) +
                        sum(g[i,t] * c_g[i,t] for i in non_str_gens) -
                        d[t] * c_d[t] +
                        sum(β⁺g[i,t] * G[i,t] for i in non_str_gens) +
                        sum(β⁺w[k,t] * W[k,t] for k in non_str_res) +
                        β⁺g_st[t] * G_ST[t] +
                        β⁺w_st[t] * W_ST[t] +
                        β⁺ch_st[t] * CH_ST[t] +
                        β⁺dch_st[t] * DCH_ST[t] +
                        β⁺d[t] * D[t]
                        for t in t_steps)
    )

    ##########OBJECTIVES OF THE VARIOUS UL AGENTS#########
    #ESS-STR
    obj_ESS = @expression(MPPDC, sum(- dch_st[t] * λ[t]
    + ch_st[t] * λ[t] for t in t_steps)
    )
    #GEN-STR
    obj_GEN = @expression(MPPDC, sum(c_g_st[t] * g_st[t] -
    g_st[t] * λ[t] for t in t_steps)
    )
    #WIND-STR
    obj_WIND = @expression(MPPDC, sum(c_w_st[t] * w_st[t]
    - w_st[t] * λ[t] for t in t_steps)
    )

    #################ADDING CONSTRAINTS#############
    @constraints(MPPDC, begin
    #ENFORCE THAT NOBODY GETS FORCE TO NEGATIVE PROFITS BY OTHER STR AGENTS
    - obj_ESS >= 0
    # - obj_GEN >= 0
    # - obj_WIND >= 0
    #ENERGY BALANCE CSTR
    en_balance[t in t_steps], balance_eq[t] == 0
    #ENFORCING STRONG DUALITY
    SDE, SDE_eq <= ϵ
    #LIMITS OF THE DISPATCHED QUANTITIES OF THE STRATEGIC AGENTS, THEY CANT BE SET FROM OUTSIDE THE MODEL
    [t in t_steps], g_st[t] <= G_ST[t]
    [t in t_steps], w_st[t] <= W_ST[t]
    [t in t_steps], dch_st[t] <= DCH_ST[t]
    [t in t_steps], ch_st[t] <= CH_ST[t]
    [t in t_steps], SOC[t] <= E_max[t]
    #FOR THE NON STR AGENTS
    [i in non_str_gens, t in t_steps], g[i,t] <= G[i,t]
    [k in non_str_res, t in t_steps], w[k,t] <= W[k,t]
    [t in t_steps], d[t] <= D[t]
    #LAGRANGIAN STATIONARITY
    #w.r.t g_{i,t}:
    [i in non_str_gens, t in t_steps], c_g[i,t] - λ[t] - β⁻g[i,t] + β⁺g[i,t] == 0
    #w.r.t g_st_{i,t}:
    [t in t_steps], c_g_st[t] - λ[t] - β⁻g_st[t] + β⁺g_st[t] == 0
    #w.r.t w_{k,t}:
    [k in non_str_res, t in t_steps], c_w[k,t] - λ[t] - β⁻w[k,t] + β⁺w[k,t] == 0
    #w.r.t w_st_{t}:
    [t in t_steps], c_w_st[t] - λ[t] - β⁻w_st[t] + β⁺w_st[t] == 0
    #w.r.t dch_st_{t}:
    [t in t_steps], c_dch_st[t] - λ[t] - β⁻dch_st[t] + β⁺dch_st[t] == 0
    #w.r.t ch_st_{t}:
    [t in t_steps], -c_ch_st[t] + λ[t] - β⁻ch_st[t] + β⁺ch_st[t] == 0
    #w.r.t d_{t}:
    [t in t_steps], -c_d[t] + λ[t] - β⁻d[t] + β⁺d[t] == 0
    ############UL CONSTRAINTS############
    #NOTE:AGENT SPECIFIC CONSTRAINTS NEED TO BE FIXED FOR EACH AGENT INDIVIDUALLY
    #RAMPING CSTRS
    #Ramping up constraint at time step 1
    ramp_up_0[t = 1], g_st[t] - g_st[end] <= ramp_limit[t]
    #Ramping up constraint
    ramp_up[t = 2:t_steps[end]], g_st[t] - g_st[t-1] <= ramp_limit[t]
    #Ramping down constraint at time step 1
    ramp_down_0[t = 1], g_st[t] - g_st[end] >= - ramp_limit[t]
    #Ramping down constraint
    ramp_down[t = 2:t_steps[end]], g_st[t] - g_st[t-1] >= - ramp_limit[t]
    end)

    @NLconstraints(MPPDC, begin
    #EXCLUSIVE CH/DCH
    # excl[t in t_steps], dch_st[t] * ch_st[t] == 0
    #SOC RULES
    soc_rule_0[t = 1], SOC[end] - SOC[t] + ch_st[t] * ηch - (dch_st[t]/ηdch) == 0
    soc_rule[t = 2:t_steps[end]], SOC[t-1] - SOC[t] + ch_st[t] * ηch - (dch_st[t]/ηdch) == 0
    end)

    return MPPDC
end

MPPDC = init_MPPDC()

###############SETTING VARIABLE BOUNDS#############
function set_bounds_UL_vars()
    for t in t_steps
        #RAMPING LIMIT OF STR GENERATOR
        JuMP.fix(MPPDC[:ramp_limit][t], Rᴳ * dat_str.GEN_STR.P_max[t]; force = true)
        #ESS EFFICIENCIES
        JuMP.fix(MPPDC[:ηch], eff_ch; force = true)
        JuMP.fix(MPPDC[:ηdch], eff_dch; force = true)
        #SETTING BOUNDS OF THE UPPER LEVEL VARS
        JuMP.set_upper_bound(MPPDC[:c_g_st][t], p_cap)
        JuMP.set_lower_bound(MPPDC[:c_g_st][t], p_floor)
        JuMP.set_upper_bound(MPPDC[:c_w_st][t], p_cap)
        JuMP.set_lower_bound(MPPDC[:c_w_st][t], p_floor)
        JuMP.set_upper_bound(MPPDC[:c_ch_st][t], p_cap)
        JuMP.set_lower_bound(MPPDC[:c_ch_st][t], p_floor)
        JuMP.set_upper_bound(MPPDC[:c_dch_st][t], p_cap)
        JuMP.set_lower_bound(MPPDC[:c_dch_st][t], p_floor)

        JuMP.set_upper_bound(MPPDC[:G_ST][t], dat_str.GEN_STR.P_max[t])
        JuMP.set_upper_bound(MPPDC[:W_ST][t], dat_str.WIND_STR.P_max[t])
        JuMP.set_upper_bound(MPPDC[:DCH_ST][t], dat_str.ESS_STR.P_max[t])
        JuMP.set_upper_bound(MPPDC[:CH_ST][t], dat_str.ESS_STR.P_max[t])
        JuMP.fix(MPPDC[:E_max][t], dat_str.ESS_STR.E_max[t]; force = true) #By this the upper bound var. on SoC is set
    end
end

#LL VARIABLE BOUNDS OF THE NON-STR AGENTS
#NOTE: FOR STR AGENTS THE BOUNDS ARE THE DECISION VARIABLES
function set_bounds_LL_vars()
    for t in t_steps
        for i in non_str_gens
            JuMP.fix(MPPDC[:c_g][i,t], dat_nonstr[i].c_g[t]; force = true) #MC OF GEN
            JuMP.fix(MPPDC[:G][i,t], dat_nonstr[i].G[t]; force = true) #MAX CAPACITY OF GEN
        end
        for k in non_str_res
            JuMP.fix(MPPDC[:c_w][k,t], dat_nonstr[k].c_g[t]; force = true) #COST OF WIND
            JuMP.fix(MPPDC[:W][k,t], dat_nonstr[k].G[t]; force = true) #MAX WIND CAPACITY
        end
        JuMP.fix(MPPDC[:c_d][t], p_cap; force = true) #WTP SET TO PRICE CAP
        JuMP.fix(MPPDC[:D][t], load[t]; force = true) #MAX LOAD
    end
end

#FIXING THE OTHER STRATEGIC AGENT VARIABLES
#NOTE: This is not a neccesity but may ease the solvers job
current = "ESS"
function fix_other_agent_vars()
    #THIS IS JUST UGLY....
    #ALLOW ONLY QUANTITY BIDDING FOR EACH AGENT
    if current == "ESS"
        for t in t_steps
            #QUANTITIES ARE EQUAL TO MAX CAPACITY
            JuMP.fix(MPPDC[:G_ST][t], dat_str.GEN_STR.P_max[t]; force = true)
            JuMP.fix(MPPDC[:W_ST][t], dat_str.WIND_STR.P_max[t]; force = true)
            #PRICES ARE EQUAL TO THE MARGINAL PRICES
            JuMP.fix(MPPDC[:c_g_st][t], dat_str.GEN_STR.MC_g[t]; force = true)
            JuMP.fix(MPPDC[:c_w_st][t], dat_str.WIND_STR.MC_g[t]; force = true)
            JuMP.fix(MPPDC[:c_dch_st][t], dat_str.ESS_STR.MC_dch[t]; force = true)
            JuMP.fix(MPPDC[:c_ch_st][t], dat_str.ESS_STR.MC_ch[t]; force = true)
        end
    elseif current == "GEN_STR"
        for t in t_steps
            #QUANTITIES ARE EQUAL TO MAX CAPACITY
            JuMP.fix(MPPDC[:W_ST][t], dat_str.WIND_STR.P_max[t]; force = true)
            JuMP.fix(MPPDC[:DCH_ST][t], dat_str.ESS_STR.P_max[t]; force = true)
            JuMP.fix(MPPDC[:CH_ST][t], dat_str.ESS_STR.P_max[t]; force = true)
            #PRICES ARE EQUAL TO THE MARGINAL PRICES
            JuMP.fix(MPPDC[:c_g_st][t], dat_str.GEN_STR.MC_g[t]; force = true)
            JuMP.fix(MPPDC[:c_w_st][t], dat_str.WIND_STR.MC_g[t]; force = true)
            # JuMP.fix(MPPDC[:c_dch_st][t], dat_str.ESS_STR.MC_dch[t]; force = true)
            # JuMP.fix(MPPDC[:c_ch_st][t], dat_str.ESS_STR.MC_ch[t]; force = true)
        end
    elseif current == "WIND_STR"
        for t in t_steps
            #QUANTITIES ARE EQUAL TO MAX CAPACITY
            JuMP.fix(MPPDC[:G_ST][t], dat_str.GEN_STR.P_max[t]; force = true)
            JuMP.fix(MPPDC[:DCH_ST][t], dat_str.ESS_STR.P_max[t]; force = true)
            JuMP.fix(MPPDC[:CH_ST][t], dat_str.ESS_STR.P_max[t]; force = true)
            #PRICES ARE EQUAL TO THE MARGINAL PRICES
            JuMP.fix(MPPDC[:c_w_st][t], dat_str.WIND_STR.MC_g[t]; force = true)
            JuMP.fix(MPPDC[:c_g_st][t], dat_str.GEN_STR.MC_g[t]; force = true)
            # JuMP.fix(MPPDC[:c_dch_st][t], dat_str.ESS_STR.MC_dch[t]; force = true)
            # JuMP.fix(MPPDC[:c_ch_st][t], dat_str.ESS_STR.MC_ch[t]; force = true)

        end
    else
        println("ERROR: The current agent is not part of the list of str agents!!!")
    end
end

setting_UL_bounds = set_bounds_UL_vars()
setting_LL_bounds = set_bounds_LL_vars()
fixing_other_strategic_vars = fix_other_agent_vars()

#NOTE: USE THE SAME MPPDC FOR ALL UL AGENTS AND FIX THE CORRESPONDING VARS AND SET THE RIGHT OBJECTIVE

function return_objective()
    #DECLARATION OF THE CURRENT OBJECTIVE
    if current == "ESS"
        obj = sum(- MPPDC[:dch_st][t] * MPPDC[:λ][t]
        + MPPDC[:ch_st][t] * MPPDC[:λ][t] for t in t_steps)

    elseif current == "GEN_STR"
        obj = sum(dat_str.GEN_STR.MC_g[t] * MPPDC[:g_st][t] -
        MPPDC[:g_st][t] * MPPDC[:λ][t] for t in t_steps)

    elseif current == "WIND_STR"
        obj = sum(dat_str.WIND_STR.MC_g[t] * MPPDC[:w_st][t]
        - MPPDC[:w_st][t] * MPPDC[:λ][t] for t in t_steps)
    else
        println("ERROR: The current agent is not part of the list of str agents!!!")
    end

    return obj
end
obj = return_objective()

#NOW SETTING THE RESPECTIVE OBJECTIVE
JuMP.set_objective_function(MPPDC, obj)
JuMP.set_objective_sense(MPPDC, MOI.MIN_SENSE)

###########SETTING THE OPTIMIZER FOR THE MODEL###############
# set_optimizer(MPPDC, with_optimizer(KNITRO.Optimizer,
#                                                 algorithm =3,
#                                                 hessian_no_f=1,
#                                                 datacheck=0,
#                                                 ms_enable= 1,
#                                                 ms_maxsolves = 100
#                                                  ))
set_optimizer(MPPDC, with_optimizer(Ipopt.Optimizer))

optimize!(MPPDC)
objective_value(MPPDC)

#########THE WHILE LOOP TO DECREASE EPSILON#########
flag_sequanctial_red = false
if flag_sequanctial_red == true
    global iter_limit = 10
    global iterations = 0
    global tol = 1e-4
    global epsilon_new = 1
    ############THE SMAL SCRIPT SNIPPET FOR UPDATING THE REGULARIZATION TERM
    while epsilon_new >= tol && iterations < iter_limit
        #########SETTING THE WARM START POINT TO THE PREVIOUS SOLUTION#########
        set_start_value.(all_variables(MPPDC), value.(all_variables(MPPDC)))
        ###################FIXING EPRSILONS TO THE REDUCED PARAMS##############
        global epsilon_new = value.(MPPDC[:ϵ])*0.1
        fix(MPPDC[:ϵ], epsilon_new; force = true)
        println("\\
        ######### \\
        THE NEW RELAXATION TERMS: \\
                    ϵ: $epsilon_new and \\
                    #ITERATIONS: $iterations \\
        #########")
        ################OPTIMIZING############
        optimize!(MPPDC)
        ##############CACLCULATING PROFIT############
        # profit = compute_profit()
        # println("PROFIT: ", profit)
        ############PROGRESSING IN TIME##############
        global iterations += 1
    end
end



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
#READING THE BIDS OF UL AGENT
#QUANTITY-PRICE BIDS
#ESS
e_ch_bid_strat = value.(MPPDC[:DCH_ST])
p_ch_bid_strat = value.(MPPDC[:c_dch_st])
e_dch_bid_strat = value.(MPPDC[:CH_ST])
p_dch_bid_strat = value.(MPPDC[:c_ch_st])
#GENCO
p_g_bid_strat = value.(MPPDC[:c_g_st])
e_g_bid_strat = value.(MPPDC[:G_ST])
#RES
e_r_bid_strat = value.(MPPDC[:W_ST])
p_r_bid_strat = value.(MPPDC[:c_w_st])

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
     sns.set_style("ticks")
     sns.set_context("paper", font_scale=3.0)
     fig = figure()
     # ax1 = subplot(3,1,1)
     stackplot(t_steps,
     gen_1[:], gen_2[:], gen_3[:], gen_5[:],
     gen_6[:], gen_7[:], gen_8[:],
     res_PV[:], dch[:], gen[:], res[:],
     edgecolor="k", linewidth = 1.5,
     labels =["GEN1", "GEN2", "GEN3","GEN5", "GEN6", "GEN7", "GEN8","PV", "DCH-ST", "GEN-ST", "WIND-ST"],
     )
     stackplot(t_steps,edgecolor="k", linewidth = 1.5, -ch[:], labels =["CH-ST"])
     plot(t_steps, load[:], label="LOAD",color ="k",linestyle=":", linewidth =3.5,marker="D", markersize = 8)
     grid()
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


# objective_function(MPPDC)
# objective_sense(MPPDC)
