"This module is created to build MPPDC with complex bids, for all UL agents,
Definitation of which UL problem is being built is based on the 'current' string passed in"
function build_MPPDC_complex(t_steps, non_str_gens, non_str_res, current, optimizer, objective_type)
    ηch = 0.95 #Needs to be passed in
    ηdch = 0.95 #Needs to be passed in
    soc_init = dat_str.ESS_STR.E_max[1] * 0.5
    non_str_agents = vcat(generators, renewables)
    #############FIRST INIT THE ABSTRACT MODEL##########
    if current == "ESS_STR"
        ####INIT ABSTRACT MODEL######
        function init_MPPDC_ESS()
            ############INIT THE MODEL INSTANCE############
            MPPDC =  Model()
            MPPDC.ext[:param] = Dict()
            param = MPPDC.ext[:param]
            MPPDC.ext[:eq] = Dict()
            eq = MPPDC.ext[:eq]
            ##############PARAMETERS IN ESS MPPDC###############
            #OTHER AGENT QUANTITY BIDS
            param[:G_ST] = @NLparameter(MPPDC, G_ST[t in t_steps] == 0) #Generation quantity bid
            param[:W_ST] = @NLparameter(MPPDC, W_ST[t in t_steps] == 0) #Wind quantity bid
            #MAX ENERGY CONTENT
            param[:E_max] = @NLparameter(MPPDC, E_max[t in t_steps] == 0) #max capacity of the ESS
            #PRICE-TAKING AGENTS QUANTITY BIDS
            param[:G] = @NLparameter(MPPDC, G[i in non_str_gens, t in t_steps] == 0) #Generation capacity of non-str generator
            param[:W] = @NLparameter(MPPDC, W[k in non_str_res, t in t_steps]  == 0) #Generation capacity of non-str RES
            param[:D] = @NLparameter(MPPDC, D[t in t_steps]  == 0) #Max load

            #PRICE BIDS, AS QUANTITY BIDDING ONLY
            #########PRICE BIDS#############
            param[:c_g_st] = @NLparameter(MPPDC, c_g_st[t in t_steps] == 0) #Generation price bid str
            param[:c_w_st] = @NLparameter(MPPDC, c_w_st[t in t_steps]  == 0)#Wind price bid str
            # @NLparameter(MPPDC, c_dch_st[t in t_steps] == 0) #ESS discharge price bid str
            # @NLparameter(MPPDC, c_ch_st[t in t_steps] == 0) #ESS charge price bid str
            param[:c_g] = @NLparameter(MPPDC, c_g[i in non_str_gens, t in t_steps] == 0) #Generation price bid nonstr
            param[:c_w] = @NLparameter(MPPDC, c_w[k in non_str_res, t in t_steps] == 0) #Wind price bid nonstr
            param[:c_d] = @NLparameter(MPPDC, c_d[t in t_steps] == 0) #WTP of the load

            ########CONSTANT TERMS FOR RAMPING###########
            param[:ramp_limit] = @NLparameter(MPPDC, ramp_limit[t in t_steps] == 1.0) #Set to 100% as default

            for t in t_steps
                #RAMPING LIMIT OF STR GENERATOR
                JuMP.set_value(ramp_limit[t], Rᴳ * dat_str.GEN_STR.P_max[t])
                #MAXIMUM ENERGY CONTENT OF ESS
                JuMP.set_value(E_max[t], dat_str.ESS_STR.E_max[t]) #By this the upper bound var. on SoC is set

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
                #QUANTITIES ARE EQUAL TO MAX CAPACITY
                JuMP.set_value(G_ST[t], dat_str.GEN_STR.P_max[t])
                JuMP.set_value(W_ST[t], dat_str.WIND_STR.P_max[t])
                #PRICES ARE EQUAL TO THE MARGINAL PRICES
                JuMP.set_value(c_g_st[t], dat_str.GEN_STR.MC_g[t])
                JuMP.set_value(c_w_st[t], dat_str.WIND_STR.MC_g[t])
                # JuMP.set_value(c_ch_st[t], dat_str.ESS_STR.MC_ch[t])
                # JuMP.set_value(c_dch_st[t], dat_str.ESS_STR.MC_dch[t])
            end

            @variables(MPPDC, begin
                #REGULARIZATION A PARAMETER
                ϵ == 0.0001
                ############UL DECISION VARS################
                ############PRICE BIDS###########
                # c_dch_st[t in t_steps]
                # c_ch_st[t in t_steps]
                #########QUANTITY BIDS############
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

                μsoc[t in t_steps] #Dual assiciated with the SoC rules
                μ_0 #Dual assiciated with the SOC init
                λ[t in t_steps] #Market price, dual of the energy balance cstr of LL
            end)


            ##########ENERGY BALANCE##########
            eq[:en_balance] = balance_eq = @expression(MPPDC, [t in t_steps],
            - g_st[t] - w_st[t] -dch_st[t] + ch_st[t]
            - sum(g[i,t] for i in non_str_gens) - sum(w[k,t] for k in non_str_res) + d[t])

            #######STRONG DUALITY EQ##########
            #PRIMAL OBJ - DUAL OBJ
            eq[:SDE_eq] = SDE_eq = @NLexpression(MPPDC,
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

            ##########OBJECTIVES OF THE VARIOUS UL AGENTS#########
            # #ESS-STR
            # obj_ESS = @NLexpression(MPPDC, sum(- dch_st[t] * λ[t]
            # + ch_st[t] * λ[t] for t in t_steps)
            # )
            #GEN-STR
            # obj_GEN = @NLexpression(MPPDC, sum(c_g_st[t] * g_st[t] -
            # g_st[t] * λ[t] for t in t_steps)
            # )
            # #WIND-STR
            # obj_WIND = @NLexpression(MPPDC, sum(c_w_st[t] * w_st[t]
            # - w_st[t] * λ[t] for t in t_steps)
            # )

            #################ADDING CONSTRAINTS#############
            @constraints(MPPDC, begin
                # ENERGY BALANCE CSTR
                en_balance[t in t_steps], balance_eq[t] == 0
                #SOC RULES
                soc_rule_0[t = 1], SOC[end] - SOC[t] + ch_st[t] * ηch - (dch_st[t]/ηdch) == 0
                soc_rule[t = 2:t_steps[end]], SOC[t-1] - SOC[t] + ch_st[t] * ηch - (dch_st[t]/ηdch) == 0
                soc_init, SOC[end] == soc_init
            end)

            @NLconstraints(MPPDC, begin
                #NOTE:DUE TO THE NL PARAMETERS MOST CONSTRAINTS HAS TO BE HANDELED AS NON-LINEAR
                ############UL CONSTRAINTS############
                #NOTE:AGENT SPECIFIC CONSTRAINTS NEED TO BE FIXED FOR EACH AGENT INDIVIDUALLY
                #RAMPING CSTRS
                # #Ramping up constraint at time step 1
                # ramp_up_0[t = 1], g_st[t] - g_st[end] <= ramp_limit[t]
                # #Ramping up constraint
                # ramp_up[t = 2:t_steps[end]], g_st[t] - g_st[t-1] <= ramp_limit[t]
                # #Ramping down constraint at time step 1
                # ramp_down_0[t = 1], g_st[t] - g_st[end] >= - ramp_limit[t]
                # #Ramping down constraint
                # ramp_down[t = 2:t_steps[end]], g_st[t] - g_st[t-1] >= - ramp_limit[t]
                #ENFORCE THAT NOBODY GETS FORCE TO NEGATIVE PROFITS BY OTHER STR AGENTS
                # - obj_ESS >= 0
                # - obj_GEN >= 0
                # - obj_WIND >= 0

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
                [t in t_steps], - λ[t] - β⁻dch_st[t] + β⁺dch_st[t] - μsoc[t]/ηdch == 0
                #w.r.t ch_st_{t}:
                [t in t_steps], + λ[t] - β⁻ch_st[t] + β⁺ch_st[t] + μsoc[t] * ηch == 0
                #w.r.t. SoC{l,t}
                [t = t_steps[end]], + μsoc[1] - μsoc[t] - κ⁻soc[t] + κ⁺soc[t] + μ_0 == 0
                [t in 1:t_steps[end-1]], + μsoc[t+1] - μsoc[t] - κ⁻soc[t] + κ⁺soc[t]== 0
                #w.r.t d_{t}:
                [t in t_steps], -c_d[t] + λ[t] - β⁻d[t] + β⁺d[t] == 0
                #ENFORCING STRONG DUALITY
                SDE, SDE_eq <= ϵ
                #NOTE: ESS CAN RESPECT IT BY USING ECONOMIC SENSE, FOR OTHERS NEEDED
                #EXCLUSIVE CH/DCH
                # excl[t in t_steps], dch_st[t] * ch_st[t] == 0
            end)

            return MPPDC
        end
        ####BUILD ABSTRACT MODEL#####
        MPPDC = init_MPPDC_ESS()

    elseif current == "GEN_STR"
        ####INIT ABSTRACT MODEL######
        function init_MPPDC_GEN()
            ############INIT THE MODEL INSTANCE############
            MPPDC =  Model()
            MPPDC.ext[:param] = Dict()
            param = MPPDC.ext[:param]
            MPPDC.ext[:eq] = Dict()
            eq = MPPDC.ext[:eq]
            ##############PARAMETERS IN ESS MPPDC###############
            #OTHER AGENT QUANTITY BIDS
            # @NLparameter(MPPDC, G_ST[t in t_steps] == 0) #Generation quantity bid
            param[:W_ST] = @NLparameter(MPPDC, W_ST[t in t_steps] == 0) #Wind quantity bid
            param[:DCH_ST] = @NLparameter(MPPDC, DCH_ST[t in t_steps] == 0) #ESS discharge quantity bid
            param[:CH_ST] = @NLparameter(MPPDC, CH_ST[t in t_steps] == 0) #ESS charge quantity bid
            #MAX ENERGY CONTENT
            param[:E_max] = @NLparameter(MPPDC, E_max[t in t_steps] == 0) #max capacity of the ESS
            #PRICE TAKING AGENTS QUANTITY BIDS
            param[:G] = @NLparameter(MPPDC, G[i in non_str_gens, t in t_steps] == 0) #Generation capacity of non-str generator
            param[:W] = @NLparameter(MPPDC, W[k in non_str_res, t in t_steps]  == 0) #Generation capacity of non-str RES
            param[:D] = @NLparameter(MPPDC, D[t in t_steps]  == 0) #Max load

            #PRICE BIDS, AS QUANTITY BIDDING ONLY
            #########PRICE BIDS#############
            param[:c_g_st] = @NLparameter(MPPDC, c_g_st[t in t_steps] == 0) #Generation price bid str
            param[:c_w_st] = @NLparameter(MPPDC, c_w_st[t in t_steps]  == 0)#Wind price bid str
            # param[:c_dch_st] = @NLparameter(MPPDC, c_dch_st[t in t_steps] == 0) #ESS discharge price bid str
            # param[:c_ch_st] = @NLparameter(MPPDC, c_ch_st[t in t_steps] == 0) #ESS charge price bid str
            param[:c_g] = @NLparameter(MPPDC, c_g[i in non_str_gens, t in t_steps] == 0) #Generation price bid nonstr
            param[:c_w] = @NLparameter(MPPDC, c_w[k in non_str_res, t in t_steps] == 0) #Wind price bid nonstr
            param[:c_d] = @NLparameter(MPPDC, c_d[t in t_steps] == 0) #WTP of the load

            ########CONSTANT TERMS FOR RAMPING###########
            param[:ramp_limit] = @NLparameter(MPPDC, ramp_limit[t in t_steps] == 1.0) #Set to 100% as default

            for t in t_steps
                #RAMPING LIMIT OF STR GENERATOR
                JuMP.set_value(ramp_limit[t], Rᴳ * dat_str.GEN_STR.P_max[t])
                #MAXIMUM ENERGY CONTENT OF ESS
                JuMP.set_value(E_max[t], dat_str.ESS_STR.E_max[t]) #By this the upper bound var. on SoC is set

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
                #QUANTITIES ARE EQUAL TO MAX CAPACITY
                # JuMP.set_value(G_ST[t], dat_str.GEN_STR.P_max[t])
                JuMP.set_value(W_ST[t], dat_str.WIND_STR.P_max[t])
                JuMP.set_value(DCH_ST[t], dat_str.ESS_STR.P_max[t])
                JuMP.set_value(CH_ST[t], dat_str.ESS_STR.P_max[t])
                # JuMP.set_value(DCH_ST[t], 0)
                # JuMP.set_value(CH_ST[t], 0)
                #PRICES ARE EQUAL TO THE MARGINAL PRICES
                JuMP.set_value(c_g_st[t], dat_str.GEN_STR.MC_g[t])
                JuMP.set_value(c_w_st[t], dat_str.WIND_STR.MC_g[t])
                # #NOTE: THIS MAY HAVE TO BE REMOVED AND LET THE SOLVER CHOSE
                # JuMP.set_value(c_dch_st[t], dat_str.ESS_STR.MC_dch[t])
                # JuMP.set_value(c_ch_st[t], dat_str.ESS_STR.MC_ch[t])
            end


            @variables(MPPDC, begin
                #REGULARIZATION A PARAMETER
                ϵ == 0.0001
                ############UL DECISION VARS################
                #########QUANTITY BIDS############
                0 <= G_ST[t in t_steps]
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

                μsoc[t in t_steps] #Dual assiciated with the SoC rules
                μ_0 #Dual assiciated with the SOC init
                λ[t in t_steps] #Market price, dual of the energy balance cstr of LL
            end)


            ##########ENERGY BALANCE##########
            eq[:en_balance] = balance_eq = @expression(MPPDC, [t in t_steps],
            - g_st[t] - w_st[t] -dch_st[t] + ch_st[t]
            - sum(g[i,t] for i in non_str_gens) - sum(w[k,t] for k in non_str_res) + d[t])

            #######STRONG DUALITY EQ##########
            #PRIMAL OBJ - DUAL OBJ
            eq[:SDE_eq] = SDE_eq = @NLexpression(MPPDC,
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

            ##########OBJECTIVES OF THE VARIOUS UL AGENTS#########
            #ESS-STR
            # obj_ESS = @NLexpression(MPPDC, sum(- dch_st[t] * λ[t]
            # + ch_st[t] * λ[t] for t in t_steps)
            # )
            # #GEN-STR
            # obj_GEN = @NLexpression(MPPDC, sum(c_g_st[t] * g_st[t] -
            # g_st[t] * λ[t] for t in t_steps)
            # )
            # #WIND-STR
            # obj_WIND = @NLexpression(MPPDC, sum(c_w_st[t] * w_st[t]
            # - w_st[t] * λ[t] for t in t_steps)
            # )

            #################ADDING CONSTRAINTS#############
            @constraints(MPPDC, begin
                # ENERGY BALANCE CSTR
                en_balance[t in t_steps], balance_eq[t] == 0
                #SOC RULES
                soc_rule_0[t = 1], SOC[end] - SOC[t] + ch_st[t] * ηch - (dch_st[t]/ηdch) == 0
                soc_rule[t = 2:t_steps[end]], SOC[t-1] - SOC[t] + ch_st[t] * ηch - (dch_st[t]/ηdch) == 0
                soc_init, SOC[end] == soc_init
            end)

            @NLconstraints(MPPDC, begin
                #NOTE:DUE TO THE NL PARAMETERS MOST CONSTRAINTS HAS TO BE HANDELED AS NON-LINEAR
                ############UL CONSTRAINTS############
                #NOTE:AGENT SPECIFIC CONSTRAINTS NEED TO BE FIXED FOR EACH AGENT INDIVIDUALLY
                #RAMPING CSTRS
                #Ramping up constraint at time step 1
                # ramp_up_0[t = 1], g_st[t] - g_st[end] <= ramp_limit[t]
                # #Ramping up constraint
                # ramp_up[t = 2:t_steps[end]], g_st[t] - g_st[t-1] <= ramp_limit[t]
                # #Ramping down constraint at time step 1
                # ramp_down_0[t = 1], g_st[t] - g_st[end] >= - ramp_limit[t]
                # #Ramping down constraint
                # ramp_down[t = 2:t_steps[end]], g_st[t] - g_st[t-1] >= - ramp_limit[t]
                #ENFORCE THAT NOBODY GETS FORCE TO NEGATIVE PROFITS BY OTHER STR AGENTS
                # - obj_ESS >= 0
                # - obj_GEN >= 0
                # - obj_WIND >= 0

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
                [t in t_steps], - λ[t] - β⁻dch_st[t] + β⁺dch_st[t] - μsoc[t]/ηdch == 0
                #w.r.t ch_st_{t}:
                [t in t_steps], + λ[t] - β⁻ch_st[t] + β⁺ch_st[t] + μsoc[t] * ηch == 0
                #w.r.t. SoC{l,t}
                [t = t_steps[end]], + μsoc[1] - μsoc[t] - κ⁻soc[t] + κ⁺soc[t] + μ_0 == 0
                [t in 1:t_steps[end-1]], + μsoc[t+1] - μsoc[t] - κ⁻soc[t] + κ⁺soc[t] == 0

                #w.r.t d_{t}:
                [t in t_steps], -c_d[t] + λ[t] - β⁻d[t] + β⁺d[t] == 0

                #ENFORCING STRONG DUALITY
                SDE, SDE_eq <= ϵ
                #NOTE: ESS CAN RESPECT IT BY USING ECONOMIC SENSE, FOR OTHERS NEEDED
                #EXCLUSIVE CH/DCH
                #NOTE: THIS MAY HAVE TO BE ADDED
                # excl[t in t_steps], dch_st[t] * ch_st[t] == 0
            end)

            return MPPDC
        end
        ####BUILD ABSTRACT MODEL#####
        MPPDC = init_MPPDC_GEN()

    elseif current == "WIND_STR"
        ####INIT ABSTRACT MODEL######
        function init_MPPDC_RES()
            ############INIT THE MODEL INSTANCE############
            MPPDC =  Model()
            MPPDC.ext[:param] = Dict()
            param = MPPDC.ext[:param]
            MPPDC.ext[:eq] = Dict()
            eq = MPPDC.ext[:eq]
            ##############PARAMETERS IN ESS MPPDC###############
            #OTHER AGENT QUANTITY BIDS
            # @NLparameter(MPPDC, G_ST[t in t_steps] == 0) #Generation quantity bid
            param[:G_ST] = @NLparameter(MPPDC, G_ST[t in t_steps] == 0) #Wind quantity bid
            param[:DCH_ST] = @NLparameter(MPPDC, DCH_ST[t in t_steps] == 0) #ESS discharge quantity bid
            param[:CH_ST] = @NLparameter(MPPDC, CH_ST[t in t_steps] == 0) #ESS charge quantity bid
            #MAX ENERGY CONTENT
            param[:E_max] = @NLparameter(MPPDC, E_max[t in t_steps] == 0) #max capacity of the ESS
            #PRICE TAKING AGENTS QUANTITY BIDS
            param[:G] = @NLparameter(MPPDC, G[i in non_str_gens, t in t_steps] == 0) #Generation capacity of non-str generator
            param[:W] = @NLparameter(MPPDC, W[k in non_str_res, t in t_steps]  == 0) #Generation capacity of non-str RES
            param[:D] = @NLparameter(MPPDC, D[t in t_steps]  == 0) #Max load

            #PRICE BIDS, AS QUANTITY BIDDING ONLY
            #########PRICE BIDS#############
            param[:c_g_st] = @NLparameter(MPPDC, c_g_st[t in t_steps] == 0) #Generation price bid str
            param[:c_w_st] = @NLparameter(MPPDC, c_w_st[t in t_steps]  == 0)#Wind price bid str
            # param[:c_dch_st] = @NLparameter(MPPDC, c_dch_st[t in t_steps] == 0) #ESS discharge price bid str
            # param[:c_ch_st] = @NLparameter(MPPDC, c_ch_st[t in t_steps] == 0) #ESS charge price bid str
            param[:c_g] = @NLparameter(MPPDC, c_g[i in non_str_gens, t in t_steps] == 0) #Generation price bid nonstr
            param[:c_w] = @NLparameter(MPPDC, c_w[k in non_str_res, t in t_steps] == 0) #Wind price bid nonstr
            param[:c_d] = @NLparameter(MPPDC, c_d[t in t_steps] == 0) #WTP of the load

            ########CONSTANT TERMS FOR RAMPING###########
            param[:ramp_limit] = @NLparameter(MPPDC, ramp_limit[t in t_steps] == 1.0) #Set to 100% as default

            for t in t_steps
                #RAMPING LIMIT OF STR GENERATOR
                JuMP.set_value(ramp_limit[t], Rᴳ * dat_str.GEN_STR.P_max[t])
                #MAXIMUM ENERGY CONTENT OF ESS
                JuMP.set_value(E_max[t], dat_str.ESS_STR.E_max[t]) #By this the upper bound var. on SoC is set

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
                #QUANTITIES ARE EQUAL TO MAX CAPACITY
                # JuMP.set_value(G_ST[t], dat_str.GEN_STR.P_max[t])
                JuMP.set_value(G_ST[t], dat_str.GEN_STR.P_max[t])
                JuMP.set_value(DCH_ST[t], dat_str.ESS_STR.P_max[t])
                JuMP.set_value(CH_ST[t], dat_str.ESS_STR.P_max[t])
                # JuMP.set_value(DCH_ST[t], 0)
                # JuMP.set_value(CH_ST[t], 0)
                #PRICES ARE EQUAL TO THE MARGINAL PRICES
                JuMP.set_value(c_g_st[t], dat_str.GEN_STR.MC_g[t])
                JuMP.set_value(c_w_st[t], dat_str.WIND_STR.MC_g[t])
                #NOTE: THIS MAY HAVE TO BE REMOVED AND LET THE SOLVER CHOSE
                # JuMP.set_value(c_dch_st[t], dat_str.ESS_STR.MC_dch[t])
                # JuMP.set_value(c_ch_st[t], dat_str.ESS_STR.MC_ch[t])
            end


            @variables(MPPDC, begin
                #REGULARIZATION A PARAMETER
                ϵ == 0.0001
                ############UL DECISION VARS################
                #########QUANTITY BIDS############
                0 <= W_ST[t in t_steps]
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
                0 <= κ⁺soc[t in t_steps] #Dual of the upper bound of the SoC cstr
                0 <= κ⁻soc[t in t_steps] #Dual of the lower bound of the SoC cstr

                μsoc[t in t_steps] #Dual assiciated with the SoC rules
                μ_0 #Dual assiciated with the SOC init

            end)


            ##########ENERGY BALANCE##########
            eq[:en_balance] = balance_eq = @expression(MPPDC, [t in t_steps],
            - g_st[t] - w_st[t] -dch_st[t] + ch_st[t]
            - sum(g[i,t] for i in non_str_gens) - sum(w[k,t] for k in non_str_res) + d[t])

            #######STRONG DUALITY EQ##########
            #PRIMAL OBJ - DUAL OBJ
            eq[:SDE_eq] = SDE_eq = @NLexpression(MPPDC,
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

            ##########OBJECTIVES OF THE VARIOUS UL AGENTS#########
            #ESS-STR
            # obj_ESS = @NLexpression(MPPDC, sum(- dch_st[t] * λ[t]
            # + ch_st[t] * λ[t] for t in t_steps)
            # )
            # #GEN-STR
            # obj_GEN = @NLexpression(MPPDC, sum(c_g_st[t] * g_st[t] -
            # g_st[t] * λ[t] for t in t_steps)
            # )
            # #WIND-STR
            # obj_WIND = @NLexpression(MPPDC, sum(c_w_st[t] * w_st[t]
            # - w_st[t] * λ[t] for t in t_steps)
            # )

            #################ADDING CONSTRAINTS#############
            @constraints(MPPDC, begin
                # ENERGY BALANCE CSTR
                en_balance[t in t_steps], balance_eq[t] == 0
                #SOC RULES
                soc_rule_0[t = 1], SOC[end] - SOC[t] + ch_st[t] * ηch - (dch_st[t]/ηdch) == 0
                soc_rule[t = 2:t_steps[end]], SOC[t-1] - SOC[t] + ch_st[t] * ηch - (dch_st[t]/ηdch) == 0
                soc_init, SOC[end] == soc_init
            end)

            @NLconstraints(MPPDC, begin
                #NOTE:DUE TO THE NL PARAMETERS MOST CONSTRAINTS HAS TO BE HANDELED AS NON-LINEAR
                ############UL CONSTRAINTS############
                #NOTE:AGENT SPECIFIC CONSTRAINTS NEED TO BE FIXED FOR EACH AGENT INDIVIDUALLY
                #RAMPING CSTRS
                # #Ramping up constraint at time step 1
                # ramp_up_0[t = 1], g_st[t] - g_st[end] <= ramp_limit[t]
                # #Ramping up constraint
                # ramp_up[t = 2:t_steps[end]], g_st[t] - g_st[t-1] <= ramp_limit[t]
                # #Ramping down constraint at time step 1
                # ramp_down_0[t = 1], g_st[t] - g_st[end] >= - ramp_limit[t]
                # #Ramping down constraint
                # ramp_down[t = 2:t_steps[end]], g_st[t] - g_st[t-1] >= - ramp_limit[t]
                #ENFORCE THAT NOBODY GETS FORCE TO NEGATIVE PROFITS BY OTHER STR AGENTS
                # - obj_ESS >= 0
                # - obj_GEN >= 0
                # - obj_WIND >= 0

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
                [t in t_steps], - λ[t] - β⁻dch_st[t] + β⁺dch_st[t] - μsoc[t]/ηdch == 0
                #w.r.t ch_st_{t}:
                [t in t_steps], + λ[t] - β⁻ch_st[t] + β⁺ch_st[t] + μsoc[t] * ηch == 0
                #w.r.t. SoC{l,t}
                [t = t_steps[end]], + μsoc[1] - μsoc[t] - κ⁻soc[t] + κ⁺soc[t] + μ_0 == 0
                [t in 1:t_steps[end-1]], + μsoc[t+1] - μsoc[t] - κ⁻soc[t] + κ⁺soc[t] == 0

                #w.r.t d_{t}:
                [t in t_steps], -c_d[t] + λ[t] - β⁻d[t] + β⁺d[t] == 0

                #ENFORCING STRONG DUALITY
                SDE, SDE_eq <= ϵ
                #NOTE: ESS CAN RESPECT IT BY USING ECONOMIC SENSE, FOR OTHERS NEEDED
                #EXCLUSIVE CH/DCH
                #NOTE: THIS MAY HAVE TO BE ADDED
                # excl[t in t_steps], dch_st[t] * ch_st[t] == 0
            end)

            return MPPDC
        end
        ####BUILD ABSTRACT MODEL#####
        MPPDC = init_MPPDC_RES()
    else
        println("ERROR: The current agent is not part of the list of str agents!!!")
    end



    ###############SETTING VARIABLE BOUNDS#############
    function set_UL_bounds()
        if current == "ESS_STR"
            for t in t_steps
                JuMP.set_upper_bound(MPPDC[:DCH_ST][t], dat_str.ESS_STR.P_max[t])
                JuMP.set_upper_bound(MPPDC[:CH_ST][t], dat_str.ESS_STR.P_max[t])
                # JuMP.set_upper_bound(MPPDC[:c_ch_st][t], p_cap)
                # JuMP.set_lower_bound(MPPDC[:c_ch_st][t], p_floor)
                # JuMP.set_upper_bound(MPPDC[:c_dch_st][t], p_cap)
                # JuMP.set_lower_bound(MPPDC[:c_dch_st][t], p_floor)
            end

        elseif current == "GEN_STR"
            for t in t_steps
                JuMP.set_upper_bound(MPPDC[:G_ST][t], dat_str.GEN_STR.P_max[t])
            end

        elseif current == "WIND_STR"
            for t in t_steps
                JuMP.set_upper_bound(MPPDC[:W_ST][t], dat_str.WIND_STR.P_max[t])
            end
        else
            println("ERROR: The current agent is not part of the list of str agents!!!")
        end
    end

    setting_UL_bounds = set_UL_bounds()

    ##########SET THE OBJECTIVE ACCORDING TO THE UL AGENT##########
    function return_objective()
        #DECLARATION OF THE CURRENT OBJECTIVE
        if current == "ESS_STR"
            if objective_type == "primal"
                obj = sum(- MPPDC[:dch_st][t] * MPPDC[:λ][t]
                + MPPDC[:ch_st][t] * MPPDC[:λ][t] for t in t_steps)
            elseif objective_type == "dual"
                obj =
                - (- sum(MPPDC[:g_st][t] * dat_str.GEN_STR.MC_g[t] for t in t_steps) #GEN_STR
                - sum(MPPDC[:w_st][t] * dat_str.WIND_STR.MC_g[t] for t in t_steps) #WIND_STR
                - sum(sum(MPPDC[:g][i,t] * dat_nonstr[i].c_g[t] for t in t_steps) for i in non_str_gens) #GENs
                - sum(sum(MPPDC[:w][k,t] * dat_nonstr[k].c_g[t] for t in t_steps) for k in non_str_res) #GENs
                + sum(MPPDC[:d][t] .* p_cap for t in t_steps)
                - sum(MPPDC[:β⁺g_st][t] * dat_str.GEN_STR.P_max[t] for t in t_steps) #GEN_STR
                - sum(MPPDC[:β⁺w_st][t] * dat_str.WIND_STR.P_max[t] for t in t_steps) #WIND_STR
                - sum(sum(MPPDC[:β⁺g][i,t] * dat_nonstr[i].G[t] for t in t_steps) for i in non_str_gens) #GENs
                - sum(sum(MPPDC[:β⁺w][k,t] * dat_nonstr[k].G[t] for t in t_steps) for k in non_str_res) #GENs
                - sum(MPPDC[:β⁺d][t] * load[t] for t in t_steps)
                - sum(MPPDC[:κ⁺soc][t] * dat_str.ESS_STR.E_max[t] for t in t_steps)
                - MPPDC[:μ_0] * dat_str.ESS_STR.E_max[1] * 0.5
                )
            else
                println("ERROR: Correct objective type has to be chosen!!!")
            end

        elseif current == "GEN_STR"
            if objective_type == "primal"
                obj = sum(dat_str.GEN_STR.MC_g[t] * MPPDC[:g_st][t] -
                MPPDC[:g_st][t] * MPPDC[:λ][t] for t in t_steps)
            elseif objective_type == "dual"
                obj =
                - (
                # - sum(MPPDC[:g_st][t] * dat_str.GEN_STR.MC_g[t] for t in t_steps) #GEN_STR
                # - sum(MPPDC[:dch_st][t] * dat_str.ESS_STR.c_dch[t] for t in t_steps) #DCH_STR
                # + sum(MPPDC[:ch_st][t] * dat_str.ESS_STR.c_ch[t] for t in t_steps) #DCH_STR
                - sum(MPPDC[:w_st][t] * dat_str.WIND_STR.MC_g[t] for t in t_steps) #WIND_STR
                - sum(sum(MPPDC[:g][i,t] * dat_nonstr[i].c_g[t] for t in t_steps) for i in non_str_gens) #GENs
                - sum(sum(MPPDC[:w][k,t] * dat_nonstr[k].c_g[t] for t in t_steps) for k in non_str_res) #GENs
                + sum(MPPDC[:d][t] .* p_cap for t in t_steps)
                # - sum(MPPDC[:β⁺g_st][t] * dat_str.GEN_STR.P_max[t] for t in t_steps) #GEN_STR
                - sum(MPPDC[:β⁺dch_st][t] * dat_str.ESS_STR.P_max[t] for t in t_steps) #WIND_STR
                - sum(MPPDC[:β⁺ch_st][t] * dat_str.ESS_STR.P_max[t] for t in t_steps) #WIND_STR
                - sum(MPPDC[:β⁺w_st][t] * dat_str.WIND_STR.P_max[t] for t in t_steps) #WIND_STR
                - sum(sum(MPPDC[:β⁺g][i,t] * dat_nonstr[i].G[t] for t in t_steps) for i in non_str_gens) #GENs
                - sum(sum(MPPDC[:β⁺w][k,t] * dat_nonstr[k].G[t] for t in t_steps) for k in non_str_res) #GENs
                - sum(MPPDC[:β⁺d][t] * load[t] for t in t_steps)
                - sum(MPPDC[:κ⁺soc][t] * dat_str.ESS_STR.E_max[t] for t in t_steps)
                - MPPDC[:μ_0] * dat_str.ESS_STR.E_max[1] * 0.5
                )
            else
                println("ERROR: Correct objective type has to be chosen!!!")
            end

        elseif current == "WIND_STR"
            if objective_type == "primal"
                obj = sum(dat_str.WIND_STR.MC_g[t] * MPPDC[:w_st][t]
                - MPPDC[:w_st][t] * MPPDC[:λ][t] for t in t_steps)
            elseif objective_type == "dual"
                obj =
                - (
                - sum(MPPDC[:g_st][t] * dat_str.GEN_STR.MC_g[t] for t in t_steps) #GEN_STR
                # - sum(MPPDC[:dch_st][t] * dat_str.ESS_STR.c_dch[t] for t in t_steps) #DCH_STR
                # + sum(MPPDC[:ch_st][t] * dat_str.ESS_STR.c_ch[t] for t in t_steps) #DCH_STR
                # - sum(MPPDC[:w_st][t] * dat_str.WIND_STR.MC_g[t] for t in t_steps) #WIND_STR
                - sum(sum(MPPDC[:g][i,t] * dat_nonstr[i].c_g[t] for t in t_steps) for i in non_str_gens) #GENs
                - sum(sum(MPPDC[:w][k,t] * dat_nonstr[k].c_g[t] for t in t_steps) for k in non_str_res) #GENs
                + sum(MPPDC[:d][t] .* p_cap for t in t_steps)
                - sum(MPPDC[:β⁺g_st][t] * dat_str.GEN_STR.P_max[t] for t in t_steps) #GEN_STR
                - sum(MPPDC[:β⁺dch_st][t] * dat_str.ESS_STR.P_max[t] for t in t_steps) #WIND_STR
                - sum(MPPDC[:β⁺ch_st][t] * dat_str.ESS_STR.P_max[t] for t in t_steps) #WIND_STR
                # - sum(MPPDC[:β⁺w_st][t] * dat_str.WIND_STR.P_max[t] for t in t_steps) #WIND_STR
                - sum(sum(MPPDC[:β⁺g][i,t] * dat_nonstr[i].G[t] for t in t_steps) for i in non_str_gens) #GENs
                - sum(sum(MPPDC[:β⁺w][k,t] * dat_nonstr[k].G[t] for t in t_steps) for k in non_str_res) #GENs
                - sum(MPPDC[:β⁺d][t] * load[t] for t in t_steps)
                - sum(MPPDC[:κ⁺soc][t] * dat_str.ESS_STR.E_max[t] for t in t_steps)
                - MPPDC[:μ_0] * dat_str.ESS_STR.E_max[1] * 0.5
                )

            else
                println("ERROR: Correct objective type has to be chosen!!!")
            end

        else
            println("ERROR: The current agent is not part of the list of str agents!!!")
        end

        return obj
    end

    obj = return_objective()

    JuMP.set_objective_function(MPPDC, obj)
    JuMP.set_objective_sense(MPPDC, MOI.MIN_SENSE)

    ##########SETTING THE OPTIMIZER#########
    set_optimizer(MPPDC, optimizer)


    return MPPDC
end
