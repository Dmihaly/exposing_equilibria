"Diagonalization for EPEC complex"

include("MILPs/MILP_complex.jl") #The MILP model with complex bidding

function run_diagonaluization(EPEC)

    function compute_profits(model)
        lambda = value.(model[:λ][:])
        #NON-STRATEGIC GENERATORS
        gen_1 = value.(model[:g][:GEN1,:])
        gen_2 = value.(model[:g][:GEN2,:])
        gen_3 = value.(model[:g][:GEN3,:])
        gen_5 = value.(model[:g][:GEN5,:])
        gen_6 = value.(model[:g][:GEN6,:])
        gen_7 = value.(model[:g][:GEN7,:])
        gen_8 = value.(model[:g][:GEN8,:])

        #DEMAND
        d = value.(model[:d])
        #NON-STRATEGIC RES
        res_PV = value.(model[:w][:PV,:])
        #DISPATCHED QUANTITIES FROMS STR AGENTS
        dch = value.(model[:dch_st])
        ch = value.(model[:ch_st])
        gen = value.(model[:g_st])
        res = value.(model[:w_st])
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

    ###########READING THE SOLUTIONS OF THE EPEC#############
    mod_dat_str = deepcopy(dat_str)
    mod_dat_str.GEN_STR.G[:] = value.(EPEC[:G_ST])
    mod_dat_str.WIND_STR.G[:] = value.(EPEC[:W_ST])
    mod_dat_str.ESS_STR.DCH[:] = value.(EPEC[:DCH_ST])
    mod_dat_str.ESS_STR.CH[:] = value.(EPEC[:CH_ST])
    #CALCULATING PROFITS OF THE EPEC
    global profit_ESS_EPEC, profit_GEN_EPEC, profit_WIND_EPEC = compute_profits(EPEC)
    global iter_diag = 0 #iteration counter
    global tol_diag = 1.0 #We consider any deviation smaller than this as negligible
    global max_iter_diag = 10

    while iter_diag <= max_iter_diag
        println("...Starting the diagonalization loop...")
        #####INIT ESS MILP######
        current = "ESS_STR"
        objective_type = "dual"
        optimizer = with_optimizer(Gurobi.Optimizer, OutputFlag=0)

        println("SOLUTION METHOD:
                    ....Agent: $current
                    ....Objective_type: $objective_type")

        #INIT MODEL
        if bidding == "complex"
            MILP_ESS = build_MILP_complex(t_steps, non_str_gens,
            non_str_res, current, optimizer, objective_type)
        elseif bidding ==  "quantity"
            MILP_ESS = build_MILP_quantity(t_steps, non_str_gens,
            non_str_res, current, optimizer, objective_type)
        else
            println("ERROR: Pick 'complex' or 'quantity' for bidding.")
        end

        #UPDATING OTHER UL AGENT DECISIONS
        for t in t_steps
            unfix(MILP_ESS[:G_ST][t])
            unfix(MILP_ESS[:W_ST][t])
            fix(MILP_ESS[:G_ST][t], mod_dat_str.GEN_STR.G[t]; force = true)
            fix(MILP_ESS[:W_ST][t], mod_dat_str.WIND_STR.G[t]; force = true)
        end

        #OPTIMIZE
        optimize!(MILP_ESS)
        objective_value(MILP_ESS)

        #READ MILP RESULT
        global profit_ESS_MILP_ESS, profit_GEN_MILP_ESS, profit_RES_MILP_ESS = compute_profits(MILP_ESS)
        println("ESS_PROFIT:$profit_ESS_MILP_ESS")

        #####INIT GEN MILP######
        current = "GEN_STR"
        println("SOLUTION METHOD:
                    ....Agent: $current
                    ....Objective_type: $objective_type")

        #INIT MODEL
        if bidding == "complex"
            MILP_GEN = build_MILP_complex(t_steps, non_str_gens,
            non_str_res, current, optimizer, objective_type)
        elseif bidding ==  "quantity"
            MILP_GEN = build_MILP_quantity(t_steps, non_str_gens,
            non_str_res, current, optimizer, objective_type)
        else
            println("ERROR: Pick 'complex' or 'quantity' for bidding.")
        end

        #UPDATING OTHER UL AGENT DECISIONS
        for t in t_steps
            unfix(MILP_GEN[:DCH_ST][t])
            unfix(MILP_GEN[:CH_ST][t])
            unfix(MILP_GEN[:W_ST][t])
            fix(MILP_GEN[:DCH_ST][t], mod_dat_str.ESS_STR.DCH[t]; force = true)
            fix(MILP_GEN[:CH_ST][t], mod_dat_str.ESS_STR.CH[t]; force = true)
            fix(MILP_GEN[:W_ST][t], mod_dat_str.WIND_STR.G[t]; force = true)
        end

        #OPTIMIZE
        optimize!(MILP_GEN)
        objective_value(MILP_GEN)

        #READ MILP RESULT
        global profit_ESS_MILP_GEN, profit_GEN_MILP_GEN, profit_RES_MILP_GEN = compute_profits(MILP_GEN)
        println("GEN_PROFIT:$profit_GEN_MILP_GEN")

        #####INIT WIND MILP######
        current = "WIND_STR"
        println("SOLUTION METHOD:
                    ....Agent: $current
                    ....Objective_type: $objective_type")

        #INIT MODEL
        if bidding == "complex"
            MILP_WIND = build_MILP_complex(t_steps, non_str_gens,
            non_str_res, current, optimizer, objective_type)
        elseif bidding ==  "quantity"
            MILP_WIND = build_MILP_quantity(t_steps, non_str_gens,
            non_str_res, current, optimizer, objective_type)
        else
            println("ERROR: Pick 'complex' or 'quantity' for bidding.")
        end

        #UPDATING OTHER UL AGENT DECISIONS
        for t in t_steps
            unfix(MILP_WIND[:DCH_ST][t])
            unfix(MILP_WIND[:CH_ST][t])
            unfix(MILP_WIND[:G_ST][t])
            fix(MILP_WIND[:DCH_ST][t], mod_dat_str.ESS_STR.DCH[t]; force = true)
            fix(MILP_WIND[:CH_ST][t], mod_dat_str.ESS_STR.CH[t]; force = true)
            fix(MILP_WIND[:G_ST][t], mod_dat_str.GEN_STR.G[t]; force = true)
        end

        #OPTIMIZE
        optimize!(MILP_WIND)
        objective_value(MILP_WIND)

        #READ MILP RESULT
        global profit_ESS_MILP_WIND, profit_GEN_MILP_WIND, profit_RES_MILP_WIND = compute_profits(MILP_WIND)
        println("WIND_PROFIT:$profit_RES_MILP_WIND")
        ######SAVE RESULTS#########
        #USE DF HERE
        println("FINAL_PROFITS:
                        PROFIT_ESS: $profit_ESS_MILP_ESS,
                        PROFIT_GEN: $profit_GEN_MILP_GEN,
                        PROFIT_RES: $profit_RES_MILP_WIND")
        ####UPDATING DECISIONS#########
        #Update mod_dat_str
        mod_dat_str.GEN_STR.G[:] = value.(MILP_GEN[:G_ST])
        mod_dat_str.WIND_STR.G[:] = value.(MILP_WIND[:W_ST])
        mod_dat_str.ESS_STR.DCH[:] = value.(MILP_ESS[:DCH_ST])
        mod_dat_str.ESS_STR.CH[:] = value.(MILP_ESS[:CH_ST])

        #####CALCULATING ERRORS########
        #W.R.T. thair dec. vars
        #NOTE: THIS MAY BE MORE REPRESENTATIVE
        # global err_bids
        # value.(MILP_GEN[:G_ST]) - value.(EPEC[:G_ST])
        global err = 0
        if iter_diag == 0
            err_ESS = profit_ESS_MILP_ESS - profit_ESS_EPEC
            err_GEN = profit_GEN_MILP_GEN - profit_GEN_EPEC
            err_WIND = profit_RES_MILP_WIND - profit_WIND_EPEC
            if err_ESS > 0
                err += err_ESS
            end
            if err_GEN > 0
                err += err_GEN
            end
            if err_WIND > 0
                err += err_WIND
            end
        else
            err_ESS = profit_ESS_MILP_ESS - profit_ESS_MILP_ESS_0
            err_GEN = profit_GEN_MILP_GEN - profit_GEN_MILP_GEN_0
            err_WIND = profit_RES_MILP_WIND - profit_RES_MILP_WIND_0
            if err_ESS > 0
                err += err_ESS
            end
            if err_GEN > 0
                err += err_GEN
            end
            if err_WIND > 0
                err += err_WIND
            end
        end
        #IF TOLERANCE HAS REACHED -> BREAK
        if err <= tol_diag
            println("Diagonalization has converged to a solution after the $iter_diag-th iteration!")
            if iter_diag == 0
                println("
                @@@HURRAY: EPEC SOLUTION WAS NASH-EQUILIBRIUM@@@
                ")
            end
            println("PROFITS:
                                PROFIT_ESS: $profit_ESS_MILP_ESS,
                                PROFIT_GEN: $profit_GEN_MILP_GEN,
                                PROFIT_RES: $profit_RES_MILP_WIND")
            println("DEVIATIONS:
                            DEV_ESS: $err_ESS,
                            DEV_GEN: $err_GEN,
                            DEV_RES: $err_WIND")
            break

        end



        println("DEVIATIONS:
                        DEV_ESS: $err_ESS,
                        DEV_GEN: $err_GEN,
                        DEV_RES: $err_WIND")

        global profit_ESS_MILP_ESS_0 = profit_ESS_MILP_ESS
        global profit_GEN_MILP_GEN_0 = profit_GEN_MILP_GEN
        global profit_RES_MILP_WIND_0 = profit_RES_MILP_WIND
        iter_diag += 1
    end
    return iter_diag, profit_ESS_MILP_ESS, profit_GEN_MILP_GEN, profit_RES_MILP_WIND
end
