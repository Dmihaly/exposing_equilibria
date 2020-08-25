"This module is created to build MPPDC with complex bids, for all UL agents,
Definitation of which UL problem is being built is based on the 'current' string passed in"
function build_MPPDC_complex(t_steps, non_str_gens, non_str_res, current, optimizer)

    non_str_agents = vcat(generators, renewables)
    ηch = 0.95 #Needs to be passed in
    ηdch = 0.95 #Needs to be passed in
    #############FIRST INIT THE ABSTRACT MODEL##########
    function init_MPPDC_complex()
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
        E_max[t in t_steps] == 0 #max capacity of the ESS
        0 <= SOC[t in t_steps] #ESS State of Charge var
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
        0 <= κ⁺soc[t in t_steps] #Dual of the upper bound of the SoC cstr
        0 <= κ⁻soc[t in t_steps] #Dual of the lower bound of the SoC cstr
        μsoc[t in t_steps] #Dual assiciated with the SoC rules
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
                            κ⁺soc[t] * E_max[t] +
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
        #w.r.t. SoC{l,t}
        [t = t_steps[end]], + μsoc[1] - μsoc[t] - κ⁻soc[t] + κ⁺soc[t] == 0
        [t in 1:t_steps[end-1]], + μsoc[t+1] - μsoc[t] - κ⁻soc[t] + κ⁺soc[t] == 0
        #w.r.t d_{t}:
        [t in t_steps], -c_d[t] + λ[t] - β⁻d[t] + β⁺d[t] == 0
        ############UL CONSTRAINTS############
        #SOC RULES
        soc_rule_0[t = 1], SOC[end] - SOC[t] + ch_st[t] * ηch - (dch_st[t]/ηdch) == 0
        soc_rule[t = 2:t_steps[end]], SOC[t-1] - SOC[t] + ch_st[t] * ηch - (dch_st[t]/ηdch) == 0
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
        excl[t in t_steps], dch_st[t] * ch_st[t] == 0
        #LS
        #w.r.t dch_st_{t}:
        [t in t_steps], c_dch_st[t] - λ[t] - β⁻dch_st[t] + β⁺dch_st[t] - μsoc[t]/ηdch == 0
        #w.r.t ch_st_{t}:
        [t in t_steps], -c_ch_st[t] + λ[t] - β⁻ch_st[t] + β⁺ch_st[t] + μsoc[t] * ηch == 0
        end)

        return MPPDC
    end


    MPPDC = init_MPPDC_complex()

    ###############SETTING VARIABLE BOUNDS#############
    function set_bounds_UL_vars()
        #ESS EFFICIENCIES
        # JuMP.fix(MPPDC[:ηch], eff_ch; force = true)
        # JuMP.fix(MPPDC[:ηdch], eff_dch; force = true)
        for t in t_steps
            #RAMPING LIMIT OF STR GENERATOR
            JuMP.fix(MPPDC[:ramp_limit][t], Rᴳ * dat_str.GEN_STR.P_max[t]; force = true)
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
    function fix_other_agent_vars()
        #THIS IS JUST UGLY....
        #ALLOW ONLY QUANTITY BIDDING FOR EACH AGENT
        if current == "ESS_STR"
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

    ##########SET THE OBJECTIVE ACCORDING TO THE UL AGENT##########
    function return_objective()
        #DECLARATION OF THE CURRENT OBJECTIVE
        if current == "ESS_STR"
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

    JuMP.set_objective_function(MPPDC, obj)
    JuMP.set_objective_sense(MPPDC, MOI.MIN_SENSE)

    ##########SETTING THE OPTIMIZER#########
    set_optimizer(MPPDC, with_optimizer(optimizer))


    return MPPDC
end
