function build_EPEC_PC_MPPDC_quantity(optimizer)
    ηch = 0.95 #Needs to be passed in
    ηdch = 0.95 #Needs to be passed in
    soc_init = dat_str.ESS_STR.E_max[1] * 0.5
    non_str_agents = vcat(generators, renewables)
    function init_EPEC_PC_quantity()
        ############INIT THE MODEL INSTANCE############
        EPEC =  Model()
        EPEC.ext[:param] = Dict()
        param = EPEC.ext[:param]
        EPEC.ext[:eq] = Dict()
        eq = EPEC.ext[:eq]
        #PRICE-TAKING AGENTS QUANTITY BIDS
        param[:G] = @variable(EPEC, G[i in non_str_gens, t in t_steps] == 0) #Generation capacity of non-str generator
        param[:W] = @variable(EPEC, W[k in non_str_res, t in t_steps]  == 0) #Generation capacity of non-str RES
        param[:D] = @variable(EPEC, D[t in t_steps]  == 0) #Max load
        #########MAX POWER OF THE UL AGENTS###########
        # param[:P_G_ST] = @variable(EPEC, P_G_ST[t in t_steps] == 0)
        # param[:P_W_ST] = @variable(EPEC, P_W_ST[t in t_steps] == 0)
        # param[:P_DCH_ST] = @variable(EPEC, P_DCH_ST[t in t_steps] == 0)
        # param[:P_CH_ST] = @variable(EPEC, P_CH_ST[t in t_steps] == 0)
        param[:E_max] = @variable(EPEC, E_max[t in t_steps] == 0) #max capacity of the ESS
        #########PRICE BIDS#############
        #NOTE:ESS HAS BOTH PRCE AND QUANTITY AS DEC VARS!
        param[:c_g_st] = @variable(EPEC, c_g_st[t in t_steps] == 0) #Generation price bid str
        param[:c_w_st] = @variable(EPEC, c_w_st[t in t_steps]  == 0)#Wind price bid str
        param[:c_g] = @variable(EPEC, c_g[i in non_str_gens, t in t_steps] == 0) #Generation price bid nonstr
        param[:c_w] = @variable(EPEC, c_w[k in non_str_res, t in t_steps] == 0) #Wind price bid nonstr
        param[:c_d] = @variable(EPEC, c_d[t in t_steps] == 0) #WTP of the load

        for t in t_steps
            #MAXIMUM ENERGY CONTENT OF ESS
            JuMP.fix(E_max[t], dat_str.ESS_STR.E_max[t]; force = true) #By this the upper bound var. on SoC is set

            for i in non_str_gens
                JuMP.fix(c_g[i,t], dat_nonstr[i].c_g[t]; force = true) #MC OF GEN
                JuMP.fix(G[i,t], dat_nonstr[i].G[t]; force = true) #MAX CAPACITY OF GEN
            end
            for k in non_str_res
                JuMP.fix(c_w[k,t], dat_nonstr[k].c_g[t]; force = true) #COST OF WIND
                JuMP.fix(W[k,t], dat_nonstr[k].G[t]; force = true) #MAX WIND CAPACITY
            end
            JuMP.fix(c_d[t], p_cap; force = true) #WTP SET TO PRICE CAP
            JuMP.fix(D[t], load[t]; force = true) #MAX LOAD
            #W.R.T OTHER STRATEGIC AGENTS
            #PRICES ARE EQUAL TO THE MARGINAL PRICES
            JuMP.fix(c_g_st[t], dat_str.GEN_STR.MC_g[t]; force = true)
            JuMP.fix(c_w_st[t], dat_str.WIND_STR.MC_g[t]; force = true)
            # JuMP.set_value(c_ch_st[t], dat_str.ESS_STR.MC_ch[t])
            # JuMP.set_value(c_dch_st[t], dat_str.ESS_STR.MC_dch[t])
        end

        @variables(EPEC, begin
            ############UL DECISION VARS################
            #REGULARIZATION A PARAMETER
            ϵ == 0.0001
            r_dch_st[t in t_steps]
            r_ch_st[t in t_steps]
            r_g_st[t in t_steps]
            r_w_st[t in t_steps]
            ############PRICE BIDS###########
            c_dch_st[t in t_steps]
            c_ch_st[t in t_steps]
            #########QUANTITY BIDS############
            0 <= G_ST[t in t_steps] #Generator str quantity bid
            0 <= W_ST[t in t_steps] #Wind str quantity bid
            0 <= DCH_ST[t in t_steps] #ESS discharge quantity bid
            0 <= CH_ST[t in t_steps] #ESS charge quantity bid
            0 <= SOC[t in t_steps] #ESS State of Charge var
            #############LL DECISION VARS###########
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
        end)

        ##########ENERGY BALANCE##########
        eq[:en_balance] = balance_eq = @expression(EPEC, [t in t_steps],
        - g_st[t] - w_st[t] -dch_st[t] + ch_st[t]
        - sum(g[i,t] for i in non_str_gens) - sum(w[k,t] for k in non_str_res) + d[t])

        #######STRONG DUALITY EQ##########
        #PRIMAL OBJ - DUAL OBJ
        eq[:SDE_eq] = SDE_eq = @expression(EPEC,
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
                            r_g_st[t] +
                            r_w_st[t] +
                            r_ch_st[t] +
                            r_dch_st[t] +
                            β⁺d[t] * D[t]
                            for t in t_steps)

        )

        @constraints(EPEC, begin
            # ENERGY BALANCE CSTR
            en_balance[t in t_steps], balance_eq[t] == 0
            #SOC RULES
            soc_rule_0[t = 1], SOC[end] - SOC[t] + ch_st[t] * ηch - (dch_st[t]/ηdch) == 0
            soc_rule[t = 2:t_steps[end]], SOC[t-1] - SOC[t] + ch_st[t] * ηch - (dch_st[t]/ηdch) == 0
            soc_0, SOC[end] == soc_init
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

            #ENFORCING STRONG DUALITY
            SDE, SDE_eq == 0
            [t in t_steps], β⁺ch_st[t] * CH_ST[t] - r_ch_st[t] <= ϵ
            [t in t_steps], -β⁺ch_st[t] * CH_ST[t] + r_ch_st[t] <= ϵ
            [t in t_steps], β⁺dch_st[t] * DCH_ST[t] - r_dch_st[t]  <= ϵ
            [t in t_steps], -β⁺dch_st[t] * DCH_ST[t] + r_dch_st[t]  <= ϵ
            [t in t_steps], β⁺g_st[t] * G_ST[t] - r_g_st[t] <= ϵ
            [t in t_steps], -β⁺g_st[t] * G_ST[t] + r_g_st[t] <= ϵ
            [t in t_steps], β⁺w_st[t] * W_ST[t] - r_w_st[t] <= ϵ
            [t in t_steps], -β⁺w_st[t] * W_ST[t] + r_w_st[t] <= ϵ

        end)

        return EPEC
    end
    EPEC = init_EPEC_PC_quantity()
    ###############SETTING VARIABLE BOUNDS#############
    function set_UL_bounds()
        for t in t_steps
            JuMP.set_upper_bound(EPEC[:DCH_ST][t], dat_str.ESS_STR.P_max[t])
            JuMP.set_upper_bound(EPEC[:CH_ST][t], dat_str.ESS_STR.P_max[t])
            JuMP.set_upper_bound(EPEC[:G_ST][t], dat_str.GEN_STR.P_max[t])
            JuMP.set_upper_bound(EPEC[:W_ST][t], dat_str.WIND_STR.P_max[t])
        end
    end

    setting_UL_bounds = set_UL_bounds()

    ##########SET THE OBJECTIVE ACCORDING TO THE UL AGENT##########
    function return_objective()
        #DECLARATION OF THE CURRENT OBJECTIVE
        # obj = (
        #     #COMMON PART
        #
        #
        #     #ESS
        #     - (
        #         - sum(EPEC[:g_st][t] * dat_str.GEN_STR.MC_g[t] for t in t_steps) #GEN_STR
        #         - sum(EPEC[:w_st][t] * dat_str.WIND_STR.MC_g[t] for t in t_steps) #WIND_STR
        #         - sum(sum(EPEC[:g][i,t] * dat_nonstr[i].c_g[t] for t in t_steps) for i in non_str_gens) #GENs
        #         - sum(sum(EPEC[:w][k,t] * dat_nonstr[k].c_g[t] for t in t_steps) for k in non_str_res) #GENs
        #         + sum(EPEC[:d][t] .* p_cap for t in t_steps)
        #         - sum(EPEC[:β⁺g_st][t] * EPEC[:G_ST][t] for t in t_steps) #GEN_STR
        #         - sum(EPEC[:β⁺w_st][t] * EPEC[:W_ST][t] for t in t_steps) #WIND_STR
        #         - sum(sum(EPEC[:β⁺g][i,t] * dat_nonstr[i].G[t] for t in t_steps) for i in non_str_gens) #GENs
        #         - sum(sum(EPEC[:β⁺w][k,t] * dat_nonstr[k].G[t] for t in t_steps) for k in non_str_res) #GENs
        #         - sum(EPEC[:β⁺d][t] * load[t] for t in t_steps)
        #         )
        #     #GEN
        #     - (
        #         - sum(EPEC[:w_st][t] * dat_str.WIND_STR.MC_g[t] for t in t_steps) #WIND_STR
        #         - sum(EPEC[:dch_st][t] * EPEC[:c_dch_st][t] for t in t_steps) #DCH_STR
        #         + sum(EPEC[:ch_st][t] * EPEC[:c_ch_st][t] for t in t_steps) #DCH_STR
        #         - sum(sum(EPEC[:g][i,t] * dat_nonstr[i].c_g[t] for t in t_steps) for i in non_str_gens) #GENs
        #         - sum(sum(EPEC[:w][k,t] * dat_nonstr[k].c_g[t] for t in t_steps) for k in non_str_res) #GENs
        #         + sum(EPEC[:d][t] .* p_cap for t in t_steps)
        #         - sum(EPEC[:β⁺dch_st][t] * EPEC[:DCH_ST][t] for t in t_steps) #WIND_STR
        #         - sum(EPEC[:β⁺ch_st][t] * EPEC[:CH_ST][t] for t in t_steps) #WIND_STR
        #         - sum(EPEC[:β⁺w_st][t] * EPEC[:W_ST][t] for t in t_steps) #WIND_STR
        #         - sum(sum(EPEC[:β⁺g][i,t] * dat_nonstr[i].G[t] for t in t_steps) for i in non_str_gens) #GENs
        #         - sum(sum(EPEC[:β⁺w][k,t] * dat_nonstr[k].G[t] for t in t_steps) for k in non_str_res) #GENs
        #         - sum(EPEC[:β⁺d][t] * load[t] for t in t_steps)
        #     )
        #     #WIND
        #     - (
        #         - sum(EPEC[:g_st][t] * dat_str.GEN_STR.MC_g[t] for t in t_steps) #GEN_STR
        #         - sum(EPEC[:dch_st][t] * EPEC[:c_dch_st][t] for t in t_steps) #DCH_STR
        #         + sum(EPEC[:ch_st][t] * EPEC[:c_ch_st][t] for t in t_steps) #DCH_STR
        #         - sum(sum(EPEC[:g][i,t] * dat_nonstr[i].c_g[t] for t in t_steps) for i in non_str_gens) #GENs
        #         - sum(sum(EPEC[:w][k,t] * dat_nonstr[k].c_g[t] for t in t_steps) for k in non_str_res) #GENs
        #         + sum(EPEC[:d][t] .* p_cap for t in t_steps)
        #         - sum(EPEC[:β⁺g_st][t] * EPEC[:G_ST][t] for t in t_steps) #GEN_STR
        #         - sum(EPEC[:β⁺dch_st][t] * EPEC[:DCH_ST][t] for t in t_steps) #WIND_STR
        #         - sum(EPEC[:β⁺ch_st][t] * EPEC[:CH_ST][t] for t in t_steps) #WIND_STR
        #         # - sum(EPEC[:β⁺w_st][t] * dat_str.WIND_STR.P_max[t] for t in t_steps) #WIND_STR
        #         - sum(sum(EPEC[:β⁺g][i,t] * dat_nonstr[i].G[t] for t in t_steps) for i in non_str_gens) #GENs
        #         - sum(sum(EPEC[:β⁺w][k,t] * dat_nonstr[k].G[t] for t in t_steps) for k in non_str_res) #GENs
        #         - sum(EPEC[:β⁺d][t] * load[t] for t in t_steps)
        #     )
        # )

        obj = (
                (
                    sum(- EPEC[:dch_st][t] * EPEC[:λ][t]
                    + EPEC[:ch_st][t] * EPEC[:λ][t] for t in t_steps)
                )
                +
                (
                    sum(dat_str.GEN_STR.MC_g[t] * EPEC[:g_st][t] -
                    EPEC[:g_st][t] * EPEC[:λ][t] for t in t_steps)
                )
                +
                (
                    obj = sum(dat_str.WIND_STR.MC_g[t] * EPEC[:w_st][t]
                    - EPEC[:w_st][t] * EPEC[:λ][t] for t in t_steps)
                )
            )

        return obj
    end

    obj = return_objective()

    JuMP.set_objective_function(EPEC, obj)
    JuMP.set_objective_sense(EPEC, MOI.MIN_SENSE)

    ##########SETTING THE OPTIMIZER#########
    set_optimizer(EPEC, optimizer)


    return EPEC
end
