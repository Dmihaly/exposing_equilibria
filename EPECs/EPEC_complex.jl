"This module is created to build the SCR-NLP version of the MPPDC with complex bids, with dual objective function, for all UL agents,
Selection of social planner's objective is based on the string SP_obj"
function build_EPEC_complex(t_steps, non_str_gens, non_str_res, str_agents, optimizer, SP_obj)
    ηch = 0.95 #Needs to be passed in
    ηdch = 0.95 #Needs to be passed in
    soc_init = dat_str.ESS_STR.E_max[1] * 0.5
    non_str_agents = vcat(generators, renewables)
    ####INIT ABSTRACT MODEL######
    println("...Initializing EPEC complex abstract model...")
    function init_EPEC_complex()
        ############INIT THE MODEL INSTANCE############
        EPEC =  Model()
        EPEC.ext[:param] = Dict()
        param = EPEC.ext[:param]
        EPEC.ext[:eq] = Dict()
        eq = EPEC.ext[:eq]
        #PRICE-TAKING AGENTS QUANTITY BIDS
        param[:G] = @NLparameter(EPEC, G[i in non_str_gens, t in t_steps] == 0) #Generation capacity of non-str generator
        param[:W] = @NLparameter(EPEC, W[k in non_str_res, t in t_steps]  == 0) #Generation capacity of non-str RES
        param[:D] = @NLparameter(EPEC, D[t in t_steps]  == 0) #Max load
        #########MAX POWER OF THE UL AGENTS###########
        param[:P_G_ST] = @NLparameter(EPEC, P_G_ST[t in t_steps] == 0)
        param[:P_W_ST] = @NLparameter(EPEC, P_W_ST[t in t_steps] == 0)
        param[:P_DCH_ST] = @NLparameter(EPEC, P_DCH_ST[t in t_steps] == 0)
        param[:P_CH_ST] = @NLparameter(EPEC, P_CH_ST[t in t_steps] == 0)
        param[:E_max] = @NLparameter(EPEC, E_max[t in t_steps] == 0) #max capacity of the ESS
        #########PRICE BIDS#############
        #NOTE:ESS HAS BOTH PRCE AND QUANTITY AS DEC VARS!
        param[:c_g_st] = @NLparameter(EPEC, c_g_st[t in t_steps] == 0) #Generation price bid str
        param[:c_w_st] = @NLparameter(EPEC, c_w_st[t in t_steps]  == 0)#Wind price bid str
        param[:c_g] = @NLparameter(EPEC, c_g[i in non_str_gens, t in t_steps] == 0) #Generation price bid nonstr
        param[:c_w] = @NLparameter(EPEC, c_w[k in non_str_res, t in t_steps] == 0) #Wind price bid nonstr
        param[:c_d] = @NLparameter(EPEC, c_d[t in t_steps] == 0) #WTP of the load
        #############UPDATING PARAMS#############
        for t in t_steps
            #MAXIMUM ENERGY CONTENT OF ESS
            JuMP.set_value(E_max[t], dat_str.ESS_STR.E_max[t]) #By this the upper bound var. on SoC is set
            #Updating max powers of UL agents
            JuMP.set_value(P_G_ST[t], dat_str.GEN_STR.P_max[t])
            JuMP.set_value(P_W_ST[t], dat_str.WIND_STR.P_max[t])
            JuMP.set_value(P_DCH_ST[t], dat_str.ESS_STR.P_max[t])
            JuMP.set_value(P_CH_ST[t], dat_str.ESS_STR.P_max[t])

            for i in non_str_gens
                JuMP.set_value(c_g[i,t], dat_nonstr[i].c_g[t]) #MC OF GEN
                JuMP.set_value(G[i,t], dat_nonstr[i].G[t]) #MAX CAPACITY OF GEN
            end
            for k in non_str_res
                JuMP.set_value(c_w[k,t], dat_nonstr[k].c_g[t]) #COST OF WIND
                JuMP.set_value(W[k,t], dat_nonstr[k].G[t]) #MAX WIND CAPACITY
            end
            JuMP.set_value(c_d[t], p_cap) #WTP SET TO PRICE CAP
            JuMP.set_value(D[t], load[t]) #MAX LOAD
            #W.R.T OTHER STRATEGIC AGENTS
            #PRICES ARE EQUAL TO THE MARGINAL PRICES
            JuMP.set_value(c_g_st[t], dat_str.GEN_STR.MC_g[t])
            JuMP.set_value(c_w_st[t], dat_str.WIND_STR.MC_g[t])
        end
        #########ADDING VARS###########
        @variables(EPEC, begin
            #REGULARIZATION A PARAMETER
            ϵ1 == 0.0001
            ϵ2 == 0.0001
            ############UL DECISION VARS################
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
            0 <= κ⁺soc[t in t_steps] #Dual of the upper bound of the SoC cstr
            0 <= κ⁻soc[t in t_steps] #Dual of the lower bound of the SoC cstr
            λ[t in t_steps] #Market price, dual of the energy balance cstr of LL
            ########SECOND STAGE DUALS#############
            0 <= Κ⁺soc[j in str_agents, t in t_steps] #Dual of the upper bound of the SoC cstr
            0 <= Κ⁻soc[j in str_agents, t in t_steps] #Dual of the lower bound of the SoC cstr
            Μsoc[j in str_agents, t in t_steps] #Dual assiciated with the SoC rules
            Μ_0[j in str_agents] #Dual assiciated with the SOC init

            0 <= ν⁺g[j in str_agents, i in non_str_gens, t in t_steps]
            0 <= ν⁻g[j in str_agents, i in non_str_gens, t in t_steps]
            0 <= ν⁺w[j in str_agents, k in non_str_res, t in t_steps]
            0 <= ν⁻w[j in str_agents, k in non_str_res, t in t_steps]
            0 <= ν⁺d[j in str_agents, t in t_steps]
            0 <= ν⁻d[j in str_agents, t in t_steps]
            0 <= ν⁺g_st[j in str_agents, t in t_steps]
            0 <= ν⁻g_st[j in str_agents, t in t_steps]
            0 <= ν⁺w_st[j in str_agents, t in t_steps]
            0 <= ν⁻w_st[j in str_agents, t in t_steps]
            0 <= ν⁺ch_st[j in str_agents, t in t_steps]
            0 <= ν⁻ch_st[j in str_agents, t in t_steps]
            0 <= ν⁺dch_st[j in str_agents, t in t_steps]
            0 <= ν⁻dch_st[j in str_agents, t in t_steps]
            #DUAL OF THE STRONG DUALITY EQ DUAL
            0 <= σ_PD[j in str_agents]
            #DUAL FOR THE MARKET PRICE (2nd stage)
            π[j in str_agents, t in t_steps]
            #DUAL OF THE SoC rules
            μsoc[t in t_steps]
            μ_0 #Dual assiciated with the SOC init
            #Duals of the LL equalities for KKTs
            ψg[j in str_agents, i in non_str_gens, t in t_steps]
            ψw[j in str_agents, k in non_str_res, t in t_steps]
            ψd[j in str_agents, t in t_steps]
            #For strategic agents:
            ψch_st[j in str_agents, t in t_steps]
            ψdch_st[j in str_agents, t in t_steps]
            ψsoc[j in str_agents, t in t_steps]
            ψg_st[j in str_agents, t in t_steps]
            ψw_st[j in str_agents, t in t_steps]
            0 <= δ⁺g[j in str_agents, i in non_str_gens, t in t_steps]
            0 <= δ⁻g[j in str_agents, i in non_str_gens, t in t_steps]
            0 <= δ⁺w[j in str_agents, k in non_str_res, t in t_steps]
            0 <= δ⁻w[j in str_agents, k in non_str_res, t in t_steps]
            0 <= δ⁺d[j in str_agents, t in t_steps]
            0 <= δ⁻d[j in str_agents, t in t_steps]
            0 <= δ⁺g_st[j in str_agents, t in t_steps]
            0 <= δ⁻g_st[j in str_agents, t in t_steps]
            0 <= δ⁺w_st[j in str_agents, t in t_steps]
            0 <= δ⁻w_st[j in str_agents, t in t_steps]
            0 <= δ⁺ch_st[j in str_agents, t in t_steps]
            0 <= δ⁻ch_st[j in str_agents, t in t_steps]
            0 <= δ⁺dch_st[j in str_agents, t in t_steps]
            0 <= δ⁻dch_st[j in str_agents, t in t_steps]
            0 <= δ⁺soc[j in str_agents, t in t_steps]
            0 <= δ⁻soc[j in str_agents, t in t_steps]
            #Duals of the UL bounds, specific to each agent
            #FOR THE ESS
            0 <= ϵ⁺dch[t in t_steps] #Dual of the upper bound of the dch limit cstr
            0 <= ϵ⁻dch[t in t_steps] #Dual of the lower bound of the dch limit cstr
            0 <= ϵ⁺ch[t in t_steps] #Dual of the upper bound of the ch limit cstr
            0 <= ϵ⁻ch[t in t_steps] #Dual of the lower bound of the ch limit cstr
            #FOR THE GEN
            0 <= ϵ⁺g[t in t_steps] #Dual of the upper bound of the gen limit cstr
            0 <= ϵ⁻g[t in t_steps] #Dual of the lower bound of the gen limit cstr
            #FOR THE GEN
            0 <= ϵ⁺w[t in t_steps] #Dual of the upper bound of the res limit cstr
            0 <= ϵ⁻w[t in t_steps] #Dual of the lower bound of the res limit cstr
            #############SLACKS#############
            #LL Slacks
            0 <= s_g[i in non_str_gens, t in t_steps]
            0 <= s_w[k in non_str_res, t in t_steps]
            0 <= s_d[t in t_steps]
            0 <= s_ch_st[t in t_steps]
            0 <= s_dch_st[t in t_steps]
            0 <= s_g_st[t in t_steps]
            0 <= s_w_st[t in t_steps]
            #UL slacks
            #Quantity bounds
            0 <= t_g_bid[t in t_steps]
            0 <= t_w_bid[t in t_steps]
            0 <= t_ch_bid[t in t_steps]
            0 <= t_dch_bid[t in t_steps]
            #SoC
            0 <= t_soc[t in t_steps]
            #THe param for regularization
        end)

        ##########ENERGY BALANCE##########
        eq[:en_balance] = balance_eq = @expression(EPEC, [t in t_steps],
        - g_st[t] - w_st[t] -dch_st[t] + ch_st[t]
        - sum(g[i,t] for i in non_str_gens) - sum(w[k,t] for k in non_str_res) + d[t])

        #######STRONG DUALITY EQ##########
        #PRIMAL OBJ - DUAL OBJ
        eq[:SDE_eq] = SDE_eq = @NLexpression(EPEC,
                            sum(
                            g_st[t] * c_g_st[t] +
                            w_st[t] * c_w_st[t] +
                            # dch_st[t] * c_dch_st[t] +
                            # - ch_st[t] * c_ch_st[t] +
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
                            + μ_0 * soc_init

        )


        ###########COMMON CONSTRAINTS#########
        #NOTE: MOVE HERE THE SHARED CSTRS
        ############PRIMAL FEASIBILITY OF THE MPPDC#####################
        @constraints(EPEC, begin
            # ENERGY BALANCE CSTR
            en_balance[t in t_steps], balance_eq[t] == 0
            #SOC RULES
            #NOTE: SOC RULES ARE NOT COMMON TO ALL AGENTS BUT THERE IS NO DIFFERENTIATION IN J VARS
            soc_rule_0[t = 1], SOC[end] - SOC[t] + ch_st[t] * ηch - (dch_st[t]/ηdch) == 0
            soc_rule[t = 2:t_steps[end]], SOC[t-1] - SOC[t] + ch_st[t] * ηch - (dch_st[t]/ηdch) == 0
            soc_0, SOC[end] == soc_init
        end)

        @NLconstraints(EPEC, begin
            #LIMITS OF THE DISPATCHED QUANTITIES OF THE STRATEGIC AGENTS, THEY CANT BE SET FROM OUTSIDE THE MODEL
            [t in t_steps], g_st[t] <= G_ST[t]
            [t in t_steps], w_st[t] <= W_ST[t]
            [t in t_steps], dch_st[t] <= DCH_ST[t]
            [t in t_steps], ch_st[t] <= CH_ST[t]
            [t in t_steps], SOC[t] <= E_max[t]
            #NOTE: THIS IS STILL NEEDED TO MAKE SURE CS WORKS AS A MULTIPLICATION
            [t in t_steps], CH_ST[t] <= P_CH_ST[t]
            [t in t_steps], DCH_ST[t] <= P_DCH_ST[t]
            [t in t_steps], G_ST[t] <= P_G_ST[t]
            [t in t_steps], W_ST[t] <= P_W_ST[t]
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
            [t in t_steps], - λ[t] - β⁻dch_st[t] + β⁺dch_st[t] - μsoc[t]/ηdch == 0
            #w.r.t ch_st_{t}:
            [t in t_steps], + λ[t] - β⁻ch_st[t] + β⁺ch_st[t] + μsoc[t] * ηch == 0
            #w.r.t. SoC{l,t}
            [t = t_steps[end]], + μsoc[1] - μsoc[t] - κ⁻soc[t] + κ⁺soc[t] + μ_0 == 0
            [t in 1:t_steps[end-1]], + μsoc[t+1] - μsoc[t] - κ⁻soc[t] + κ⁺soc[t] == 0
            #w.r.t d_{t}:
            [t in t_steps], -c_d[t] + λ[t] - β⁻d[t] + β⁺d[t] == 0

            #INHERITING FROM LL
            [i in non_str_gens, t in t_steps], G[i,t] -  g[i,t] - s_g[i,t] == 0
            [t in t_steps], G_ST[t] -  g_st[t] - s_g_st[t] == 0
            [k in non_str_res, t in t_steps], W[k,t] -  w[k,t] - s_w[k,t] == 0
            [t in t_steps], W_ST[t] -  w_st[t] - s_w_st[t] == 0
            [t in t_steps], DCH_ST[t] -  dch_st[t] - s_dch_st[t] == 0
            [t in t_steps], CH_ST[t] -  ch_st[t] - s_ch_st[t] == 0
            [t in t_steps], D[t] -  d[t] - s_d[t] == 0
            #########ADDING THE SLACKS TO THE CS CONSTRAINTS#######
            # #ENFORCING STRONG DUALITY
            (ϵ1 - SDE_eq) >= 0 #PRIMAL FEASIBILITY PART
            #NOTE: ESS CAN RESPECT IT BY USING ECONOMIC SENSE, FOR OTHERS NEEDED
            #EXCLUSIVE CH/DCH
            # excl[t in t_steps], dch_st[t] * ch_st[t] == 0
        end)

        #################FOR THE ESS##########
            j = "ESS_STR"

            ##############CONSTRAINTS INHERITING FROM SECOND-STAGE#########
            @NLconstraints(EPEC, begin
                ############LAGRANGIAN STATIONARITY OF THE SECOND STAGE KKTs#########
                #w.r.t DCH_st{t}:
                [t in t_steps],  - ϵ⁻dch[t] + ϵ⁺dch[t] - ν⁺dch_st[j,t] + σ_PD[j] * β⁺dch_st[t] == 0
                #w.r.t CH_st{t}:
                [t in t_steps],  - ϵ⁻ch[t] + ϵ⁺ch[t] - ν⁺ch_st[j,t] + σ_PD[j] * β⁺ch_st[t] == 0
                #w.r.t dch_st{t}:
                [t in t_steps],  - π[j,t] - ν⁻dch_st[j,t] + ν⁺dch_st[j,t]  - (Μsoc[j,t]/ηdch)  == 0
                #w.r.t ch_st{t}:
                [t in t_steps],  + π[j,t] - ν⁻ch_st[j,t] + ν⁺ch_st[j,t] + (Μsoc[j,t]*ηch) == 0
                #w.r.t. SoC{l,t}
                [t = t_steps[end]], + Μsoc[j,1] - Μsoc[j,t] - Κ⁻soc[j,t] + Κ⁺soc[j,t] + Μ_0[j] == 0
                [t in 1:t_steps[end-1]], + Μsoc[j,t+1] - Μsoc[j,t] - Κ⁻soc[j,t] + Κ⁺soc[j,t] == 0
                #w.r.t g_st{t}:
                [t in t_steps], + c_g_st[t] - π[j,t] - ν⁻g_st[j,t] + ν⁺g_st[j,t] + σ_PD[j] * c_g_st[t] == 0
                #w.r.t g{t}:
                [i in non_str_gens, t in t_steps], + c_g[i,t] - π[j,t] - ν⁻g[j,i,t] + ν⁺g[j,i,t] + σ_PD[j] * c_g[i,t] == 0
                #w.r.t w_st{t}:
                [t in t_steps], + c_w_st[t] - π[j,t] - ν⁻w_st[j,t] + ν⁺w_st[j,t] + σ_PD[j] * c_w_st[t] == 0
                #w.r.t w{t}:
                [k in non_str_res, t in t_steps], + c_w[k,t] - π[j,t] - ν⁻w[j,k,t] + ν⁺w[j,k,t] + σ_PD[j] * c_w[k,t] == 0
                #w.r.t d{t}:
                [t in t_steps], - c_d[t] + π[j,t] - ν⁻d[j,t] + ν⁺d[j,t] - σ_PD[j] * c_d[t] == 0

                #W.R.T. DUALS
                #w.r.t λ{t}:
                [t in t_steps], + ψd[j,t] - ψg_st[j,t]  - ψw_st[j,t] - sum(ψg[j,i,t] for i in non_str_gens) - sum(ψw[j,k,t] for k in non_str_res) + ψch_st[j,t] - ψdch_st[j,t] == 0
                #w.r.t β⁺g_st{t}:
                [t in t_steps], + G_ST[t] + σ_PD[j] * G_ST[t] + ψg_st[j,t] - δ⁺g_st[j,t] == 0
                #w.r.t β⁻g_st{t}:
                [t in t_steps], - ψg_st[j,t] - δ⁻g_st[j,t] == 0
                #w.r.t β⁺g{i,t}:
                [i in non_str_gens, t in t_steps], + G[i,t] + σ_PD[j] * G[i,t] + ψg[j,i,t] - δ⁺g[j,i,t] == 0
                #w.r.t β⁻g{i,t}:
                [i in non_str_gens, t in t_steps], - ψg[j,i,t] - δ⁻g[j,i,t] == 0
                #w.r.t β⁺w_st{t}:
                [t in t_steps], + W_ST[t] + σ_PD[j] * W_ST[t] + ψw_st[j,t] - δ⁺w_st[j,t] == 0
                #w.r.t β⁻w_st{t}:
                [t in t_steps], - ψw_st[j,t] - δ⁻w_st[j,t] == 0
                #w.r.t β⁺w{k,t}:
                [k in non_str_res, t in t_steps], + W[k,t] + σ_PD[j] * W[k,t] + ψw[j,k,t] - δ⁺w[j,k,t] == 0
                #w.r.t β⁻w{k,t}:
                [k in non_str_res, t in t_steps], - ψw[j,k,t] - δ⁻w[j,k,t] == 0
                #w.r.t β⁺dch_st{t}:
                [t in t_steps], σ_PD[j] * DCH_ST[t] + ψdch_st[j,t] - δ⁺dch_st[j,t] == 0
                #w.r.t β⁻dch_st{t}:
                [t in t_steps], - ψdch_st[j,t] - δ⁻dch_st[j,t] == 0
                #w.r.t β⁺ch_st{t}:
                [t in t_steps], σ_PD[j] * CH_ST[t] + ψch_st[j,t] - δ⁺ch_st[j,t] == 0
                #w.r.t β⁻ch_st{t}:
                [t in t_steps], - ψch_st[j,t] - δ⁻ch_st[j,t] == 0
                #w.r.t κ⁺soc{t}:
                [t in t_steps], + E_max[t] + σ_PD[j] * E_max[t] + ψsoc[j,t] - δ⁺soc[j,t] == 0
                #w.r.t κ⁻soc{t}:
                [t in t_steps], - ψsoc[j,t] - δ⁻soc[j,t] == 0
                #w.r.t β⁺d{t}:
                [t in t_steps], + D[t] + σ_PD[j] * D[t] + ψd[j,t] - δ⁺d[j,t] == 0
                #w.r.t β⁻d{t}:
                [t in t_steps], - ψd[j,t] - δ⁻d[j,t] == 0
                #w.r.t μsoc{t}:
                [t = 1], - (ψdch_st[j,t]/ηdch) + ψch_st[j,t] * ηch + ψsoc[j,end] - ψsoc[j,t] == 0
                [t = 2:t_steps[end]], - (ψdch_st[j,t]/ηdch) + ψch_st[j,t] * ηch + ψsoc[j,t-1] - ψsoc[j,t] == 0
                #w.r.t μ0:
                + soc_init + σ_PD[j] * soc_init + ψsoc[j,end] == 0
                #########ADDING THE SLACKS TO THE CS CONSTRAINTS#######
                #INHERITING FROM UL
                #PRIMAL FEASIBILITY FROM THE UL
                #NOTE: THIS IS STILL NEEDED TO MAKE SURE CS WORKS AS A MULTIPLICATION
                [t in t_steps], CH_ST[t] <= P_CH_ST[t]
                [t in t_steps], DCH_ST[t] <= P_DCH_ST[t]
                #FROM QUANTITY BOUNDS
                [t in t_steps], P_CH_ST[t] - CH_ST[t] - t_ch_bid[t] == 0
                [t in t_steps], P_DCH_ST[t] - DCH_ST[t] - t_dch_bid[t] == 0
                [t in t_steps], E_max[t] - SOC[t] - t_soc[t] == 0


            end)
            ############COMPLEMENTARITY SLACKNESS##########
            @NLconstraints(EPEC, begin
                #RESULTING FROM THE SDE INEQUALITY
                (ϵ1 - SDE_eq) * σ_PD[j] <= ϵ2
                #LL related
                [i in non_str_gens, t in t_steps], g[i,t] * ν⁻g[j,i,t]  <= ϵ2
                [i in non_str_gens, t in t_steps], s_g[i,t] * ν⁺g[j,i,t]  <= ϵ2
                [i in non_str_gens, t in t_steps], β⁻g[i,t] *  δ⁻g[j,i,t]  <= ϵ2
                [i in non_str_gens, t in t_steps], β⁺g[i,t] *  δ⁺g[j,i,t] <= ϵ2
                [k in non_str_res, t in t_steps], w[k,t] * ν⁻w[j,k,t]  <= ϵ2
                [k in non_str_res, t in t_steps], s_w[k,t] * ν⁺w[j,k,t]  <= ϵ2
                [k in non_str_res, t in t_steps], β⁻w[k,t] *  δ⁻w[j,k,t]  <= ϵ2
                [k in non_str_res, t in t_steps], β⁺w[k,t] *  δ⁺w[j,k,t] <= ϵ2
                [t in t_steps], g_st[t] * ν⁻g_st[j,t]  <= ϵ2
                [t in t_steps], s_g_st[t] * ν⁺g_st[j,t]  <= ϵ2
                [t in t_steps], β⁻g_st[t] *  δ⁻g_st[j,t]  <= ϵ2
                [t in t_steps], β⁺g_st[t] *  δ⁺g_st[j,t] <= ϵ2
                [t in t_steps], w_st[t] * ν⁻w_st[j,t]  <= ϵ2
                [t in t_steps], s_w_st[t] * ν⁺w_st[j,t]  <= ϵ2
                [t in t_steps], β⁻w_st[t] *  δ⁻w_st[j,t]  <= ϵ2
                [t in t_steps], β⁺w_st[t] *  δ⁺w_st[j,t] <= ϵ2
                [t in t_steps], dch_st[t] * ν⁻dch_st[j,t]  <= ϵ2
                [t in t_steps], s_dch_st[t] * ν⁺dch_st[j,t]  <= ϵ2
                [t in t_steps], β⁻dch_st[t] *  δ⁻dch_st[j,t]  <= ϵ2
                [t in t_steps], β⁺dch_st[t] *  δ⁺dch_st[j,t] <= ϵ2
                [t in t_steps], ch_st[t] * ν⁻ch_st[j,t]  <= ϵ2
                [t in t_steps], s_ch_st[t] * ν⁺ch_st[j,t]  <= ϵ2
                [t in t_steps], β⁻ch_st[t] *  δ⁻ch_st[j,t]  <= ϵ2
                [t in t_steps], β⁺ch_st[t] *  δ⁺ch_st[j,t] <= ϵ2
                [t in t_steps], d[t] * ν⁻d[j,t]  <= ϵ2
                [t in t_steps], s_d[t] * ν⁺d[j,t]  <= ϵ2
                [t in t_steps], β⁻d[t] *  δ⁻d[j,t]  <= ϵ2
                [t in t_steps], β⁺d[t] *  δ⁺d[j,t] <= ϵ2
                [t in t_steps], κ⁺soc[t] * δ⁺soc[j,t] <= ϵ2
                [t in t_steps], κ⁻soc[t] * δ⁻soc[j,t] <= ϵ2
                #UL related
                #Quantity related
                [t in t_steps], t_ch_bid[t] * ϵ⁺ch[t]  <= ϵ2
                [t in t_steps], CH_ST[t] * ϵ⁻ch[t]  <= ϵ2
                [t in t_steps], t_dch_bid[t] * ϵ⁺dch[t]  <= ϵ2
                [t in t_steps], DCH_ST[t] * ϵ⁻dch[t]  <= ϵ2
                [t in t_steps], t_soc[t] * Κ⁺soc[j,t]  <= ϵ2
                [t in t_steps], SOC[t] * Κ⁻soc[j,t] <= ϵ2
            end)

        #################FOR THE GENERATOR##########
            j = "GEN_STR"

            ##############CONSTRAINTS INHERITING FROM SECOND-STAGE#########
            @NLconstraints(EPEC, begin
                ############LAGRANGIAN STATIONARITY OF THE SECOND STAGE KKTs#########
                #W.R.T. PRIMALS
                # #w.r.t G_st{t}:
                [t in t_steps],  - ϵ⁻g[t] + ϵ⁺g[t] - ν⁺g_st[j,t] + σ_PD[j] * β⁺g_st[t] == 0
                #w.r.t dch_st{t}:
                [t in t_steps], - π[j,t] - ν⁻dch_st[j,t] + ν⁺dch_st[j,t]  - (Μsoc[j,t]/ηdch)  == 0
                #w.r.t ch_st{t}:
                [t in t_steps], + π[j,t] - ν⁻ch_st[j,t] + ν⁺ch_st[j,t] + (Μsoc[j,t]*ηch) == 0
                #w.r.t. SoC{l,t}
                [t = t_steps[end]], + Μsoc[j,1] - Μsoc[j,t] - Κ⁻soc[j,t] + Κ⁺soc[j,t] + Μ_0[j] == 0
                [t in 1:t_steps[end-1]], + Μsoc[j,t+1] - Μsoc[j,t] - Κ⁻soc[j,t] + Κ⁺soc[j,t] == 0
                #w.r.t g_st{t}:
                [t in t_steps], + dat_str.GEN_STR.MC_g[t] - π[j,t] - ν⁻g_st[j,t] + ν⁺g_st[j,t] + σ_PD[j] * c_g_st[t] == 0
                #w.r.t g{t}:
                [i in non_str_gens, t in t_steps], + c_g[i,t] - π[j,t] - ν⁻g[j,i,t] + ν⁺g[j,i,t] + σ_PD[j] * c_g[i,t] == 0
                #w.r.t w_st{t}:
                [t in t_steps], + c_w_st[t] - π[j,t] - ν⁻w_st[j,t] + ν⁺w_st[j,t] + σ_PD[j] * c_w_st[t] == 0
                #w.r.t w{t}:
                [k in non_str_res, t in t_steps], + c_w[k,t] - π[j,t] - ν⁻w[j,k,t] + ν⁺w[j,k,t] + σ_PD[j] * c_w[k,t] == 0
                #w.r.t d{t}:
                [t in t_steps], - c_d[t] + π[j,t] - ν⁻d[j,t] + ν⁺d[j,t] - σ_PD[j] * c_d[t] == 0

                #W.R.T. DUALS
                #w.r.t λ{t}:
                [t in t_steps], + ψd[j,t] - ψg_st[j,t]  - ψw_st[j,t] - sum(ψg[j,i,t] for i in non_str_gens) - sum(ψw[j,k,t] for k in non_str_res) + ψch_st[j,t] - ψdch_st[j,t] == 0
                #w.r.t β⁺g_st{t}:
                [t in t_steps], σ_PD[j] * G_ST[t] + ψg_st[j,t] - δ⁺g_st[j,t] == 0
                #w.r.t β⁻g_st{t}:
                [t in t_steps], - ψg_st[j,t] - δ⁻g_st[j,t] == 0
                #w.r.t β⁺g{i,t}:
                [i in non_str_gens, t in t_steps], + G[i,t] + σ_PD[j] * G[i,t] + ψg[j,i,t] - δ⁺g[j,i,t] == 0
                #w.r.t β⁻g{i,t}:
                [i in non_str_gens, t in t_steps], - ψg[j,i,t] - δ⁻g[j,i,t] == 0
                #w.r.t β⁺w_st{t}:
                [t in t_steps], + W_ST[t] + σ_PD[j] * W_ST[t] + ψw_st[j,t] - δ⁺w_st[j,t] == 0
                #w.r.t β⁻w_st{t}:
                [t in t_steps], - ψw_st[j,t] - δ⁻w_st[j,t] == 0
                #w.r.t β⁺w{k,t}:
                [k in non_str_res, t in t_steps], + W[k,t] +  σ_PD[j] * W[k,t] + ψw[j,k,t] - δ⁺w[j,k,t] == 0
                #w.r.t β⁻w{k,t}:
                [k in non_str_res, t in t_steps], - ψw[j,k,t] - δ⁻w[j,k,t] == 0
                #w.r.t β⁺dch_st{t}:
                [t in t_steps], + DCH_ST[t] + σ_PD[j] * DCH_ST[t] + ψdch_st[j,t] - δ⁺dch_st[j,t] == 0
                #w.r.t β⁻dch_st{t}:
                [t in t_steps], - ψdch_st[j,t] - δ⁻dch_st[j,t] == 0
                #w.r.t β⁺ch_st{t}:
                [t in t_steps], + CH_ST[t] + σ_PD[j] * CH_ST[t] + ψch_st[j,t] - δ⁺ch_st[j,t] == 0
                #w.r.t β⁻ch_st{t}:
                [t in t_steps], - ψch_st[j,t] - δ⁻ch_st[j,t] == 0
                #w.r.t κ⁺soc{t}:
                [t in t_steps], + E_max[t] + σ_PD[j] * E_max[t] + ψsoc[j,t] - δ⁺soc[j,t] == 0
                #w.r.t κ⁻soc{t}:
                [t in t_steps], - ψsoc[j,t] - δ⁻soc[j,t] == 0
                #w.r.t β⁺d{t}:
                [t in t_steps], + D[t] + σ_PD[j] * D[t] + ψd[j,t] - δ⁺d[j,t] == 0
                #w.r.t β⁻d{t}:
                [t in t_steps], - ψd[j,t] - δ⁻d[j,t] == 0
                #w.r.t μsoc{t}:
                [t = 1], - (ψdch_st[j,t]/ηdch) + ψch_st[j,t] * ηch + ψsoc[j,end] - ψsoc[j,t] == 0
                [t = 2:t_steps[end]], - (ψdch_st[j,t]/ηdch) + ψch_st[j,t] * ηch + ψsoc[j,t-1] - ψsoc[j,t] == 0
                #w.r.t μ0:
                + soc_init + σ_PD[j] * soc_init + ψsoc[j,end] == 0
                #########ADDING THE SLACKS TO THE CS CONSTRAINTS#######
                #INHERITING FROM UL
                #PRIMAL FEASIBILITY FROM THE UL
                #NOTE: THIS IS STILL NEEDED TO MAKE SURE CS WORKS AS A MULTIPLICATION
                [t in t_steps], G_ST[t] <= P_G_ST[t]
                #FROM QUANTITY BOUNDS
                [t in t_steps], P_G_ST[t] - G_ST[t] - t_g_bid[t] == 0

            end)
            ############COMPLEMENTARITY SLACKNESS##########
            @NLconstraints(EPEC, begin
                #RESULTING FROM THE SDE INEQUALITY
                (ϵ1 - SDE_eq) * σ_PD[j] <= ϵ2
                #LL related
                [i in non_str_gens, t in t_steps], g[i,t] * ν⁻g[j,i,t]  <= ϵ2
                [i in non_str_gens, t in t_steps], s_g[i,t] * ν⁺g[j,i,t]  <= ϵ2
                [i in non_str_gens, t in t_steps], β⁻g[i,t] *  δ⁻g[j,i,t]  <= ϵ2
                [i in non_str_gens, t in t_steps], β⁺g[i,t] *  δ⁺g[j,i,t] <= ϵ2
                [k in non_str_res, t in t_steps], w[k,t] * ν⁻w[j,k,t]  <= ϵ2
                [k in non_str_res, t in t_steps], s_w[k,t] * ν⁺w[j,k,t]  <= ϵ2
                [k in non_str_res, t in t_steps], β⁻w[k,t] *  δ⁻w[j,k,t]  <= ϵ2
                [k in non_str_res, t in t_steps], β⁺w[k,t] *  δ⁺w[j,k,t] <= ϵ2
                [t in t_steps], g_st[t] * ν⁻g_st[j,t]  <= ϵ2
                [t in t_steps], s_g_st[t] * ν⁺g_st[j,t]  <= ϵ2
                [t in t_steps], β⁻g_st[t] *  δ⁻g_st[j,t]  <= ϵ2
                [t in t_steps], β⁺g_st[t] *  δ⁺g_st[j,t] <= ϵ2
                [t in t_steps], w_st[t] * ν⁻w_st[j,t]  <= ϵ2
                [t in t_steps], s_w_st[t] * ν⁺w_st[j,t]  <= ϵ2
                [t in t_steps], β⁻w_st[t] *  δ⁻w_st[j,t]  <= ϵ2
                [t in t_steps], β⁺w_st[t] *  δ⁺w_st[j,t] <= ϵ2
                [t in t_steps], dch_st[t] * ν⁻dch_st[j,t]  <= ϵ2
                [t in t_steps], s_dch_st[t] * ν⁺dch_st[j,t]  <= ϵ2
                [t in t_steps], β⁻dch_st[t] *  δ⁻dch_st[j,t]  <= ϵ2
                [t in t_steps], β⁺dch_st[t] *  δ⁺dch_st[j,t] <= ϵ2
                [t in t_steps], ch_st[t] * ν⁻ch_st[j,t]  <= ϵ2
                [t in t_steps], s_ch_st[t] * ν⁺ch_st[j,t]  <= ϵ2
                [t in t_steps], β⁻ch_st[t] *  δ⁻ch_st[j,t]  <= ϵ2
                [t in t_steps], β⁺ch_st[t] *  δ⁺ch_st[j,t] <= ϵ2
                [t in t_steps], d[t] * ν⁻d[j,t]  <= ϵ2
                [t in t_steps], s_d[t] * ν⁺d[j,t]  <= ϵ2
                [t in t_steps], β⁻d[t] *  δ⁻d[j,t]  <= ϵ2
                [t in t_steps], β⁺d[t] *  δ⁺d[j,t] <= ϵ2
                [t in t_steps], κ⁺soc[t] * δ⁺soc[j,t] <= ϵ2
                [t in t_steps], κ⁻soc[t] * δ⁻soc[j,t] <= ϵ2
                #UL related
                #Quantity related
                [t in t_steps], t_g_bid[t] * ϵ⁺g[t]  <= ϵ2
                [t in t_steps], G_ST[t] * ϵ⁻g[t]  <= ϵ2
                [t in t_steps], t_soc[t] * Κ⁺soc[j,t]  <= ϵ2
                [t in t_steps], SOC[t] * Κ⁻soc[j,t] <= ϵ2
            end)


        #################FOR THE WIND##########
            j = "WIND_STR"
            ##############CONSTRAINTS INHERITING FROM SECOND-STAGE#########
            @NLconstraints(EPEC, begin
                ############LAGRANGIAN STATIONARITY OF THE SECOND STAGE KKTs#########
                #W.R.T. PRIMALS
                # #w.r.t W_st{t}:
                [t in t_steps],  - ϵ⁻w[t] + ϵ⁺w[t] - ν⁺w_st[j,t] + σ_PD[j] * β⁺w_st[t] == 0
                #w.r.t dch_st{t}:
                [t in t_steps], - π[j,t] - ν⁻dch_st[j,t] + ν⁺dch_st[j,t]  - (Μsoc[j,t]/ηdch)  == 0
                #w.r.t ch_st{t}:
                [t in t_steps], + π[j,t] - ν⁻ch_st[j,t] + ν⁺ch_st[j,t] + (Μsoc[j,t]*ηch) == 0
                #w.r.t. SoC{l,t}
                [t = t_steps[end]], + Μsoc[j,1] - Μsoc[j,t] - Κ⁻soc[j,t] + Κ⁺soc[j,t] + Μ_0[j] == 0
                [t in 1:t_steps[end-1]], + Μsoc[j,t+1] - Μsoc[j,t] - Κ⁻soc[j,t] + Κ⁺soc[j,t] == 0
                #w.r.t g_st{t}:
                [t in t_steps], + c_g_st[t] - π[j,t] - ν⁻g_st[j,t] + ν⁺g_st[j,t] + σ_PD[j] * c_g_st[t] == 0
                #w.r.t g{t}:
                [i in non_str_gens, t in t_steps], + c_g[i,t] - π[j,t] - ν⁻g[j,i,t] + ν⁺g[j,i,t] + σ_PD[j] * c_g[i,t] == 0
                #w.r.t w_st{t}:
                [t in t_steps], + dat_str.WIND_STR.MC_g[t] - π[j,t] - ν⁻w_st[j,t] + ν⁺w_st[j,t] + σ_PD[j] * c_w_st[t] == 0
                #w.r.t w{t}:
                [k in non_str_res, t in t_steps], + c_w[k,t] - π[j,t] - ν⁻w[j,k,t] + ν⁺w[j,k,t] + σ_PD[j] * c_w[k,t] == 0
                #w.r.t d{t}:
                [t in t_steps], - c_d[t] + π[j,t] - ν⁻d[j,t] + ν⁺d[j,t] - σ_PD[j] * c_d[t] == 0

                #W.R.T. DUALS
                #w.r.t λ{t}:
                [t in t_steps], + ψd[j,t] - ψg_st[j,t]  - ψw_st[j,t] - sum(ψg[j,i,t] for i in non_str_gens) - sum(ψw[j,k,t] for k in non_str_res) + ψch_st[j,t] - ψdch_st[j,t] == 0
                #w.r.t β⁺g_st{t}:
                [t in t_steps], +G_ST[t] + σ_PD[j] * G_ST[t] + ψg_st[j,t] - δ⁺g_st[j,t] == 0
                #w.r.t β⁻g_st{t}:
                [t in t_steps], - ψg_st[j,t] - δ⁻g_st[j,t] == 0
                #w.r.t β⁺g{i,t}:
                [i in non_str_gens, t in t_steps], + G[i,t] + σ_PD[j] * G[i,t] + ψg[j,i,t] - δ⁺g[j,i,t] == 0
                #w.r.t β⁻g{i,t}:
                [i in non_str_gens, t in t_steps], - ψg[j,i,t] - δ⁻g[j,i,t] == 0
                #w.r.t β⁺w_st{t}:
                [t in t_steps], + σ_PD[j] * W_ST[t] + ψw_st[j,t] - δ⁺w_st[j,t] == 0
                #w.r.t β⁻w_st{t}:
                [t in t_steps], - ψw_st[j,t] - δ⁻w_st[j,t] == 0
                #w.r.t β⁺w{k,t}:
                [k in non_str_res, t in t_steps], + W[k,t] +  σ_PD[j] * W[k,t] + ψw[j,k,t] - δ⁺w[j,k,t] == 0
                #w.r.t β⁻w{k,t}:
                [k in non_str_res, t in t_steps], - ψw[j,k,t] - δ⁻w[j,k,t] == 0
                #w.r.t β⁺dch_st{t}:
                [t in t_steps], + DCH_ST[t] + σ_PD[j] * DCH_ST[t] + ψdch_st[j,t] - δ⁺dch_st[j,t] == 0
                #w.r.t β⁻dch_st{t}:
                [t in t_steps], - ψdch_st[j,t] - δ⁻dch_st[j,t] == 0
                #w.r.t β⁺ch_st{t}:
                [t in t_steps], + CH_ST[t] + σ_PD[j] * CH_ST[t] + ψch_st[j,t] - δ⁺ch_st[j,t] == 0
                #w.r.t β⁻ch_st{t}:
                [t in t_steps], - ψch_st[j,t] - δ⁻ch_st[j,t] == 0
                #w.r.t κ⁺soc{t}:
                [t in t_steps], + E_max[t] + σ_PD[j] * E_max[t] + ψsoc[j,t] - δ⁺soc[j,t] == 0
                #w.r.t κ⁻soc{t}:
                [t in t_steps], - ψsoc[j,t] - δ⁻soc[j,t] == 0
                #w.r.t β⁺d{t}:
                [t in t_steps], + D[t] + σ_PD[j] * D[t] + ψd[j,t] - δ⁺d[j,t] == 0
                #w.r.t β⁻d{t}:
                [t in t_steps], - ψd[j,t] - δ⁻d[j,t] == 0
                #w.r.t μsoc{t}:
                [t = 1], - (ψdch_st[j,t]/ηdch) + ψch_st[j,t] * ηch + ψsoc[j,end] - ψsoc[j,t] == 0
                [t = 2:t_steps[end]], - (ψdch_st[j,t]/ηdch) + ψch_st[j,t] * ηch + ψsoc[j,t-1] - ψsoc[j,t] == 0
                #w.r.t μ0:
                + soc_init + σ_PD[j] * soc_init + ψsoc[j,end] == 0
                #########ADDING THE SLACKS TO THE CS CONSTRAINTS#######
                #INHERITING FROM UL
                #PRIMAL FEASIBILITY FROM THE UL
                #NOTE: THIS IS STILL NEEDED TO MAKE SURE CS WORKS AS A MULTIPLICATION
                [t in t_steps], W_ST[t] <= P_W_ST[t]
                #FROM QUANTITY BOUNDS
                [t in t_steps], P_W_ST[t] - W_ST[t] - t_w_bid[t] == 0

            end)
            ############COMPLEMENTARITY SLACKNESS##########
            @NLconstraints(EPEC, begin
                #RESULTING FROM THE SDE INEQUALITY
                (ϵ1 - SDE_eq) * σ_PD[j] <= ϵ2
                #LL related
                [i in non_str_gens, t in t_steps], g[i,t] * ν⁻g[j,i,t]  <= ϵ2
                [i in non_str_gens, t in t_steps], s_g[i,t] * ν⁺g[j,i,t]  <= ϵ2
                [i in non_str_gens, t in t_steps], β⁻g[i,t] *  δ⁻g[j,i,t]  <= ϵ2
                [i in non_str_gens, t in t_steps], β⁺g[i,t] *  δ⁺g[j,i,t] <= ϵ2
                [k in non_str_res, t in t_steps], w[k,t] * ν⁻w[j,k,t]  <= ϵ2
                [k in non_str_res, t in t_steps], s_w[k,t] * ν⁺w[j,k,t]  <= ϵ2
                [k in non_str_res, t in t_steps], β⁻w[k,t] *  δ⁻w[j,k,t]  <= ϵ2
                [k in non_str_res, t in t_steps], β⁺w[k,t] *  δ⁺w[j,k,t] <= ϵ2
                [t in t_steps], g_st[t] * ν⁻g_st[j,t]  <= ϵ2
                [t in t_steps], s_g_st[t] * ν⁺g_st[j,t]  <= ϵ2
                [t in t_steps], β⁻g_st[t] *  δ⁻g_st[j,t]  <= ϵ2
                [t in t_steps], β⁺g_st[t] *  δ⁺g_st[j,t] <= ϵ2
                [t in t_steps], w_st[t] * ν⁻w_st[j,t]  <= ϵ2
                [t in t_steps], s_w_st[t] * ν⁺w_st[j,t]  <= ϵ2
                [t in t_steps], β⁻w_st[t] *  δ⁻w_st[j,t]  <= ϵ2
                [t in t_steps], β⁺w_st[t] *  δ⁺w_st[j,t] <= ϵ2
                [t in t_steps], dch_st[t] * ν⁻dch_st[j,t]  <= ϵ2
                [t in t_steps], s_dch_st[t] * ν⁺dch_st[j,t]  <= ϵ2
                [t in t_steps], β⁻dch_st[t] *  δ⁻dch_st[j,t]  <= ϵ2
                [t in t_steps], β⁺dch_st[t] *  δ⁺dch_st[j,t] <= ϵ2
                [t in t_steps], ch_st[t] * ν⁻ch_st[j,t]  <= ϵ2
                [t in t_steps], s_ch_st[t] * ν⁺ch_st[j,t]  <= ϵ2
                [t in t_steps], β⁻ch_st[t] *  δ⁻ch_st[j,t]  <= ϵ2
                [t in t_steps], β⁺ch_st[t] *  δ⁺ch_st[j,t] <= ϵ2
                [t in t_steps], d[t] * ν⁻d[j,t]  <= ϵ2
                [t in t_steps], s_d[t] * ν⁺d[j,t]  <= ϵ2
                [t in t_steps], β⁻d[t] *  δ⁻d[j,t]  <= ϵ2
                [t in t_steps], β⁺d[t] *  δ⁺d[j,t] <= ϵ2
                [t in t_steps], κ⁺soc[t] * δ⁺soc[j,t] <= ϵ2
                [t in t_steps], κ⁻soc[t] * δ⁻soc[j,t] <= ϵ2
                #UL related
                #Quantity related
                [t in t_steps], t_w_bid[t] * ϵ⁺w[t]  <= ϵ2
                [t in t_steps], W_ST[t] * ϵ⁻w[t]  <= ϵ2
                [t in t_steps], t_soc[t] * Κ⁺soc[j,t]  <= ϵ2
                [t in t_steps], SOC[t] * Κ⁻soc[j,t] <= ϵ2
            end)

        return EPEC
    end

    EPEC = init_EPEC_complex()
    ##########GETTING THE SP OBJECTIVE

    function get_EPEC_obj()
        if SP_obj == "empty"
            obj = 0
        elseif SP_obj == "competitive"
            #COMPETITIVE OBJECTIVE CORRESPONDS TO THE TRUE SW MAXIMUM
            #CALCULATING SOCIAL WELFARE
            obj =(
                    #Maximize consumers surplus
                    - sum(EPEC[:d][t] .* p_cap for t in t_steps)
                    #Minimize generation cost
                    + sum(EPEC[:g_st][t] * dat_str.GEN_STR.MC_g[t] for t in t_steps) #GEN_STR
                    + sum(EPEC[:w_st][t] * dat_str.WIND_STR.MC_g[t] for t in t_steps) #WIND_STR
                    + sum(sum(EPEC[:g][i,t] * dat_nonstr[i].c_g[t] for t in t_steps) for i in non_str_gens) #GENs
                    + sum(sum(EPEC[:w][k,t] * dat_nonstr[k].c_g[t] for t in t_steps) for k in non_str_res) #GENs
                )

        elseif SP_obj == "fav_demand"
            obj =(
                    #Minimize payment
                    #NOTE: As demand is always fully satisfied this can be linearized by using the maximum demand
                    sum(EPEC[:λ][t] .* load[t] for t in t_steps)
                )

        elseif SP_obj == "collusive"
            # #COLLUSIVE EQUILIBRIUM BY MAXIMIZING PRODUCERS PROFITS
            # obj = (
            # sum(- EPEC[:dch_st][t] * EPEC[:λ][t]
            # + EPEC[:ch_st][t] * EPEC[:λ][t] for t in t_steps)
            # +
            # sum(dat_str.GEN_STR.MC_g[t] * EPEC[:g_st][t] -
            # EPEC[:g_st][t] * EPEC[:λ][t] for t in t_steps)
            # +
            # sum(dat_str.GEN_STR.MC_g[t] * EPEC[:w_st][t] -
            # EPEC[:w_st][t] * EPEC[:λ][t] for t in t_steps)
            # )

            obj = (
                #ESS
                - (
                    - sum(EPEC[:g_st][t] * dat_str.GEN_STR.MC_g[t] for t in t_steps) #GEN_STR
                    - sum(EPEC[:w_st][t] * dat_str.WIND_STR.MC_g[t] for t in t_steps) #WIND_STR
                    - sum(sum(EPEC[:g][i,t] * dat_nonstr[i].c_g[t] for t in t_steps) for i in non_str_gens) #GENs
                    - sum(sum(EPEC[:w][k,t] * dat_nonstr[k].c_g[t] for t in t_steps) for k in non_str_res) #GENs
                    + sum(EPEC[:d][t] .* p_cap for t in t_steps)
                    - sum(EPEC[:β⁺g_st][t] * EPEC[:G_ST][t] for t in t_steps) #GEN_STR
                    - sum(EPEC[:β⁺w_st][t] * EPEC[:W_ST][t] for t in t_steps) #WIND_STR
                    - sum(sum(EPEC[:β⁺g][i,t] * dat_nonstr[i].G[t] for t in t_steps) for i in non_str_gens) #GENs
                    - sum(sum(EPEC[:β⁺w][k,t] * dat_nonstr[k].G[t] for t in t_steps) for k in non_str_res) #GENs
                    - sum(EPEC[:β⁺d][t] * load[t] for t in t_steps)
                    - sum(EPEC[:κ⁺soc][t] * dat_str.ESS_STR.E_max[t] for t in t_steps)
                    )
                #GEN
                - (
                    - sum(EPEC[:w_st][t] * dat_str.WIND_STR.MC_g[t] for t in t_steps) #WIND_STR
                    - sum(sum(EPEC[:g][i,t] * dat_nonstr[i].c_g[t] for t in t_steps) for i in non_str_gens) #GENs
                    - sum(sum(EPEC[:w][k,t] * dat_nonstr[k].c_g[t] for t in t_steps) for k in non_str_res) #GENs
                    + sum(EPEC[:d][t] .* p_cap for t in t_steps)
                    - sum(EPEC[:β⁺dch_st][t] * EPEC[:DCH_ST][t] for t in t_steps) #WIND_STR
                    - sum(EPEC[:β⁺ch_st][t] * EPEC[:CH_ST][t] for t in t_steps) #WIND_STR
                    - sum(EPEC[:β⁺w_st][t] * EPEC[:W_ST][t] for t in t_steps) #WIND_STR
                    - sum(sum(EPEC[:β⁺g][i,t] * dat_nonstr[i].G[t] for t in t_steps) for i in non_str_gens) #GENs
                    - sum(sum(EPEC[:β⁺w][k,t] * dat_nonstr[k].G[t] for t in t_steps) for k in non_str_res) #GENs
                    - sum(EPEC[:β⁺d][t] * load[t] for t in t_steps)
                    - sum(EPEC[:κ⁺soc][t] * dat_str.ESS_STR.E_max[t] for t in t_steps)
                )
                #WIND
                - (
                    - sum(EPEC[:g_st][t] * dat_str.GEN_STR.MC_g[t] for t in t_steps) #GEN_STR
                    - sum(sum(EPEC[:g][i,t] * dat_nonstr[i].c_g[t] for t in t_steps) for i in non_str_gens) #GENs
                    - sum(sum(EPEC[:w][k,t] * dat_nonstr[k].c_g[t] for t in t_steps) for k in non_str_res) #GENs
                    + sum(EPEC[:d][t] .* p_cap for t in t_steps)
                    - sum(EPEC[:β⁺g_st][t] * EPEC[:G_ST][t] for t in t_steps) #GEN_STR
                    - sum(EPEC[:β⁺dch_st][t] * EPEC[:DCH_ST][t] for t in t_steps) #WIND_STR
                    - sum(EPEC[:β⁺ch_st][t] * EPEC[:CH_ST][t] for t in t_steps) #WIND_STR
                    # - sum(EPEC[:β⁺w_st][t] * dat_str.WIND_STR.P_max[t] for t in t_steps) #WIND_STR
                    - sum(sum(EPEC[:β⁺g][i,t] * dat_nonstr[i].G[t] for t in t_steps) for i in non_str_gens) #GENs
                    - sum(sum(EPEC[:β⁺w][k,t] * dat_nonstr[k].G[t] for t in t_steps) for k in non_str_res) #GENs
                    - sum(EPEC[:β⁺d][t] * load[t] for t in t_steps)
                    - sum(EPEC[:κ⁺soc][t] * dat_str.ESS_STR.E_max[t] for t in t_steps)
                )
            )

        elseif SP_obj == "ESS_fav"
            #THE OBJECTIVE FUNCTION THAT MAXIMIZES THE ESS PROFIT
            obj = sum(- EPEC[:dch_st][t] * EPEC[:λ][t]
            + EPEC[:ch_st][t] * EPEC[:λ][t] for t in t_steps)
            #NOTE: THIS IS THE LINEARIZED VERSION
            # obj =
            #     - (- sum(EPEC[:g_st][t] * dat_str.GEN_STR.MC_g[t] for t in t_steps) #GEN_STR
            #     - sum(EPEC[:w_st][t] * dat_str.WIND_STR.MC_g[t] for t in t_steps) #WIND_STR
            #     - sum(sum(EPEC[:g][i,t] * dat_nonstr[i].c_g[t] for t in t_steps) for i in non_str_gens) #GENs
            #     - sum(sum(EPEC[:w][k,t] * dat_nonstr[k].c_g[t] for t in t_steps) for k in non_str_res) #GENs
            #     + sum(EPEC[:d][t] .* p_cap for t in t_steps)
            #     - sum(EPEC[:β⁺g_st][t] * EPEC[:G_ST][t] for t in t_steps) #GEN_STR
            #     - sum(EPEC[:β⁺w_st][t] * EPEC[:W_ST][t] for t in t_steps) #WIND_STR
            #     - sum(sum(EPEC[:β⁺g][i,t] * dat_nonstr[i].G[t] for t in t_steps) for i in non_str_gens) #GENs
            #     - sum(sum(EPEC[:β⁺w][k,t] * dat_nonstr[k].G[t] for t in t_steps) for k in non_str_res) #GENs
            #     - sum(EPEC[:β⁺d][t] * load[t] for t in t_steps)
            #     - sum(EPEC[:κ⁺soc][t] * dat_str.ESS_STR.E_max[t] for t in t_steps)
            # )

        elseif SP_obj == "GEN_fav"
            #FAVORING THE CONVENTIONAL GENERATOR
            # obj = (
            #     sum(dat_str.GEN_STR.MC_g[t] * EPEC[:g_st][t] -
            #     EPEC[:g_st][t] * EPEC[:λ][t] for t in t_steps)
            # )

            obj = (
                #GEN
                - (
                    - sum(EPEC[:w_st][t] * dat_str.WIND_STR.MC_g[t] for t in t_steps) #WIND_STR
                    - sum(sum(EPEC[:g][i,t] * dat_nonstr[i].c_g[t] for t in t_steps) for i in non_str_gens) #GENs
                    - sum(sum(EPEC[:w][k,t] * dat_nonstr[k].c_g[t] for t in t_steps) for k in non_str_res) #GENs
                    + sum(EPEC[:d][t] .* p_cap for t in t_steps)
                    - sum(EPEC[:β⁺dch_st][t] * dat_str.ESS_STR.P_max[t] for t in t_steps) #WIND_STR
                    - sum(EPEC[:β⁺ch_st][t] * dat_str.ESS_STR.P_max[t] for t in t_steps) #WIND_STR
                    - sum(EPEC[:β⁺w_st][t] * dat_str.WIND_STR.P_max[t] for t in t_steps) #WIND_STR
                    - sum(sum(EPEC[:β⁺g][i,t] * dat_nonstr[i].G[t] for t in t_steps) for i in non_str_gens) #GENs
                    - sum(sum(EPEC[:β⁺w][k,t] * dat_nonstr[k].G[t] for t in t_steps) for k in non_str_res) #GENs
                    - sum(EPEC[:β⁺d][t] * load[t] for t in t_steps)
                    - sum(EPEC[:κ⁺soc][t] * dat_str.ESS_STR.E_max[t] for t in t_steps)
                )
            )

        else
            println("ERROR:The chosen Social Planner objective is not part of the list!!!")
        end

        return obj
    end

    println("...Setting social planner's objective: $SP_obj...")
    obj_EPEC = get_EPEC_obj()
    JuMP.set_objective_function(EPEC, obj_EPEC)
    JuMP.set_objective_sense(EPEC, MOI.MIN_SENSE)

    ##########SETTING THE OPTIMIZER#########
    set_optimizer(EPEC, optimizer)

    return EPEC
end
