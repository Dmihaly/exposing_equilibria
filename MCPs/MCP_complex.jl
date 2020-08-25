function build_MCP(t_steps, non_str_gens, non_str_res, optimizer)
    ηch = 0.95 #Needs to be passed in
    ηdch = 0.95 #Needs to be passed in
    soc_init = dat_str.ESS_STR.E_max[1] * 0.5
    non_str_agents = vcat(generators, renewables)
    function init_MCP_complex()
        MCP = Model()
        MCP.ext[:param] = Dict()
        param = MCP.ext[:param]
        MCP.ext[:eq] = Dict()
        eq = MCP.ext[:eq]
        #########ADDING VARIABLES#########
        # #PRICE-TAKING AGENTS QUANTITY BIDS
        # param[:G] = @variable(MCP, G[i in non_str_gens, t in t_steps] == 0) #Generation capacity of non-str generator
        # param[:W] = @variable(MCP, W[k in non_str_res, t in t_steps]  == 0) #Generation capacity of non-str RES
        # param[:D] = @variable(MCP, D[t in t_steps]  == 0) #Max load

        #PRICE BIDS, AS QUANTITY BIDDING ONLY
        #########PRICE BIDS#############
        param[:c_g_st] = @variable(MCP, c_g_st[t in t_steps] == 0) #Generation price bid str
        param[:c_w_st] = @variable(MCP, c_w_st[t in t_steps]  == 0)#Wind price bid str
        # @NLparameter(MCP, c_dch_st[t in t_steps] == 0) #ESS discharge price bid str
        # @NLparameter(MCP, c_ch_st[t in t_steps] == 0) #ESS charge price bid str
        param[:c_g] = @variable(MCP, c_g[i in non_str_gens, t in t_steps] == 0) #Generation price bid nonstr
        param[:c_w] = @variable(MCP, c_w[k in non_str_res, t in t_steps] == 0) #Wind price bid nonstr
        param[:c_d] = @variable(MCP, c_d[t in t_steps] == 0) #WTP of the load

        @variables(MCP, begin
            ############UL DECISION VARS################
            ############PRICE BIDS###########
            #NOTE: ONLY IN CASE OF P-Q BIDDING OF THE ESS
            # c_dch_st[t in t_steps]
            # c_ch_st[t in t_steps]
            #########QUANTITY BIDS############
            0 <= DCH_ST[t in t_steps] #ESS discharge quantity bid
            0 <= CH_ST[t in t_steps] #ESS charge quantity bid
            0 <= G_ST[t in t_steps]
            0 <= W_ST[t in t_steps]
            0 <= E_max[t in t_steps] #max capacity of the ESS
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
            μsoc[t in t_steps] #Dual assiciated with the SoC rules
            μ_0 #Dual assiciated with the SOC init
            0 <= SOC[t in t_steps] #ESS State of Charge var
        end)
        ##########ENERGY BALANCE##########
        eq[:en_balance] = balance_eq = @expression(MCP, [t in t_steps],
        - g_st[t] - w_st[t] -dch_st[t] + ch_st[t]
        - sum(g[i,t] for i in non_str_gens) - sum(w[k,t] for k in non_str_res) + d[t])
        #######STRONG DUALITY EQ##########
        #PRIMAL OBJ - DUAL OBJ
        eq[:SDE_eq] = SDE_eq = @expression(MCP,
                            sum(
                            g_st[t] * dat_str.GEN_STR.MC_g[t]  +
                            w_st[t] * dat_str.WIND_STR.MC_g[t] +
                            # dch_st[t] * c_dch_st[t] +
                            # - ch_st[t] * c_ch_st[t] +
                            sum(w[k,t] * dat_nonstr[k].c_g[t] for k in non_str_res) +
                            sum(g[i,t] * dat_nonstr[i].c_g[t] for i in non_str_gens) -
                            d[t] .* p_cap +
                            sum(β⁺g[i,t] * dat_nonstr[i].G[t] for i in non_str_gens) +
                            sum(β⁺w[k,t] * dat_nonstr[k].G[t]  for k in non_str_res) +
                            β⁺g_st[t] * dat_str.GEN_STR.P_max[t] +
                            β⁺w_st[t] * dat_str.WIND_STR.P_max[t] +
                            β⁺ch_st[t] * dat_str.ESS_STR.P_max[t] +
                            β⁺dch_st[t] * dat_str.ESS_STR.P_max[t] +
                            κ⁺soc[t] * dat_str.ESS_STR.E_max[t] +
                            β⁺d[t] * load[t]
                            for t in t_steps)
                            + μ_0 * soc_init
        )
        #################ADDING CONSTRAINTS#############
        @constraints(MCP, begin
            # ENERGY BALANCE CSTR
            en_balance[t in t_steps], balance_eq[t] == 0
            #SOC RULES
            soc_rule_0[t = 1], SOC[end] - SOC[t] + ch_st[t] * ηch - (dch_st[t]/ηdch) == 0
            soc_rule[t = 2:t_steps[end]], SOC[t-1] - SOC[t] + ch_st[t] * ηch - (dch_st[t]/ηdch) == 0
            soc_0, SOC[end] == soc_init
            #LIMITS OF THE DISPATCHED QUANTITIES OF THE STRATEGIC AGENTS, THEY CANT BE SET FROM OUTSIDE THE MODEL
            [t in t_steps], g_st[t] <= dat_str.GEN_STR.P_max[t]
            [t in t_steps], w_st[t] <= dat_str.WIND_STR.P_max[t]
            [t in t_steps], dch_st[t] <= dat_str.ESS_STR.P_max[t]
            [t in t_steps], ch_st[t] <= dat_str.ESS_STR.P_max[t]
            [t in t_steps], SOC[t] <= dat_str.ESS_STR.E_max[t]
            #FOR THE NON STR AGENTS
            [i in non_str_gens, t in t_steps], g[i,t] <= dat_nonstr[i].G[t]
            [k in non_str_res, t in t_steps], w[k,t] <= dat_nonstr[k].G[t]
            [t in t_steps], d[t] <= load[t]
            #LAGRANGIAN STATIONARITY
            #w.r.t g_{i,t}:
            [i in non_str_gens, t in t_steps], dat_nonstr[i].c_g[t] - λ[t] - β⁻g[i,t] + β⁺g[i,t] == 0
            #w.r.t g_st_{i,t}:
            [t in t_steps], dat_str.GEN_STR.MC_g[t] - λ[t] - β⁻g_st[t] + β⁺g_st[t] == 0
            #w.r.t w_{k,t}:
            [k in non_str_res, t in t_steps], dat_nonstr[k].c_g[t] - λ[t] - β⁻w[k,t] + β⁺w[k,t] == 0
            #w.r.t w_st_{t}:
            [t in t_steps], dat_str.WIND_STR.MC_g[t] - λ[t] - β⁻w_st[t] + β⁺w_st[t] == 0
            #w.r.t dch_st_{t}:
            [t in t_steps], - λ[t] - β⁻dch_st[t] + β⁺dch_st[t] - μsoc[t]/ηdch == 0
            #w.r.t ch_st_{t}:
            [t in t_steps], + λ[t] - β⁻ch_st[t] + β⁺ch_st[t] + μsoc[t] * ηch == 0
            #w.r.t. SoC{l,t}
            [t = t_steps[end]], + μsoc[1] - μsoc[t] - κ⁻soc[t] + κ⁺soc[t] + μ_0 == 0
            [t in 1:t_steps[end-1]], + μsoc[t+1] - μsoc[t] - κ⁻soc[t] + κ⁺soc[t]== 0
            #w.r.t d_{t}:
            [t in t_steps], - p_cap + λ[t] - β⁻d[t] + β⁺d[t] == 0
            #ENFORCING STRONG DUALITY
            SDE, SDE_eq == 0
        end)

        return MCP
    end
    MCP = init_MCP_complex()
    ##########SETTING THE OPTIMIZER#########
    set_optimizer(MCP, optimizer)

    return MCP
end
