# date: 9/13/24
# author: Yuezhe Li 
# purpose of this code: to remove T cells in spleen & BM

using DifferentialEquations, ModelingToolkit 

include("param.jl");

function TDB_homo7(; name)
    @independent_variables t
    Dt = Differential(t)

    pars = @parameters begin 
        ## B cells 
        kBapop = 0.02                   # Apoptosis rate of B-cells [1/d]
        kBkill = 275.19                 # killing rate of B-cells [1/d]
        fBkill = 0.1                    # kBkill is multiplied by this factor (<1) to capture the less efficient killing of tissue B-cells
        kBprolif = 0.05                 # proliferation rate of B cells [1/d]
        kBtiss3exit = 0.03
        B19no20_B1920_ratio = 0.25      # CD19+CD20- B-cell:CD19+CD20+ B-cell ratio in BM
        kBmat_kBapop_ratio = 0.25       # The ratio of maturation rate of CD19+CD20- B-cells in BM to B-cell apoptosis rate
        ## TDB : B-cell kill
        VmB = 0.95                      # maximum rate of B-cell killing 
        KmTB_kill = 0.75                # the Tact/B ratio that the rate if half of VmB
        fKmTB_kill = 10.0               # KmTB_kill is multiplied by this factor (>1) to compensate for the fact the tissues are not well-stirred environments like PB 
        nkill = 1.03                    # Hill coefficient for the Tact/B portion of the Michaelis-Menten equation
        KdrugB = 1.3                    # the TDB conc at which B-cell killing rate if half of VmB [ng/mL]
        ## T cells 
        kTprolif = 0.7                  # fraction of activated cells that proliferate
        kTaexit = 0.12                  # trafficking rate of activated CD8+ T-cells from PB to tissue [1/d]
        kTact = 9.83                    # conversion rate of resting or post-activated CD8+ T cells (and vice versa) [1/d]
        fTadeact = 0.01                 # fraction of activated CD8+ T-cells that deactivate to resting or "exhausted" CD8+ T-cells
        fTap = 3.77                     # Ratio of the partition coeff for activated CD8+ T-cells to that of resting CD8+ T-cells
        kTaapop = 0.06                  # apoptosis rate of CD69+ CD8+ T-cells (activated CD8+ T-cells) [1/d]
        fTaprolif = 2.11                # proliferation rate of activated CD8+ T-cells [1/d]
        fTrapop = 0.2                   # Ratio of the apoptosis rate of resting CD8+ T-cells to that of CD69+CD8+ T-cells
        kTrexit = 0.05                  # trafficking rate of resting/ post-activated CD8+ T-cells from PB to tissue [1/d]
        fAICD = 1.5                     # Effect of activation-induced cell death (AICD) on CD69+CD8+ T-cells
        fTa0apop = 2.02                 # Ratio of the apoptosis rate of post-activated CD8+ T-cells to that of CD69+CD8+ T-cells
        fTa0deact = 0.                  # Conversion rate of post-activated CD8+ T-cells to resting CD8+ T-cells [1/d]
        kTgen = 1.
        ## T cells : TDB activation 
        VmT = 0.9                       # maximum rate of T-cell activation [unitless]
        fTact = 0.25                    # VmT is multiplied by this factor (<1) to capture the less efficient activation of tissue T-cells
        KmBT_act = 0.72                 # The B/T rario, at which the rate is half of VmT (at high concentrations of TDB)
        ndrugactT = 0.8                 # Hill coefficient for the TDB portion of the Michaelis-Menten equation
        S = 1.4                         # Hill coefficient for the B:T portion of the Michaelis-Menten equation
        KdrugactT = 130.08              # the TDB conc at which T-cell activation was half of VmT [ng/mL]
        ## IL6 
        thalfIL6 = 1/3                  # [hr]
        ## other 
        tinjhalf = 1.
        Trpbo_perml = 2.0E6             # [1/mL]
        Bpbo_perml = 1.0E6              # [1/mL]
        act0on = 1.
        tissue2on = 1.
        tissue3on = 1.
        Trpbref_perml = 2.0e6           # baseline conc of resting CD8+ T-cells in PB [cells/mL]
        Bpbref_perml = 0.5e6            # baseline conc of CD19+CD20+ B-cells in PB [cells/mL]
        Kp2 = 0.07
        Kptumor = 0.05
        # tumor 
        tumor_on = 0.
        Init_Vtumor = 49.6              # [mL] 
        kBtumorprolif = 0.025           # proliferation rate of B-cells in tumor [1/d]
        kIL6prod = 14.0                 # production rate of IL6 [pg/d] (homo)
        Vpb = 5000.                     # physiological volume of peripheral blood [mL]
        Vtissue2 = 400.                 # total physiological volume of LNs [mL]
        Vtissue3 = 500.0                # Physiology volume of BM [mL]
        KTrp2 = 190.
        KTrptumor = 125.
        KBp2 = 190.
        KBp3 = 80.0
        KBptumor = 6.5E3
        ## PK   
        Vm_tdb = 212                    # [μg/d]
        Km_tdb = 4.74                   # [μg/d]
        kabs_TDB = 0.26
        CL_TDB = 637.5              # [mL/day]
        Q_TDB = 1812.75             # [mL/day]
        V1_TDB = 2758.5             # [mL]
        V2_TDB = 13010.25           # [mL]
    end 
     
    vars = @variables begin 
        ## T cells
        restTpb(t) = Trpbo_perml * Vpb
        actTpb(t) = 0. 
        act0Tpb(t) = 0.
        act0Ttiss(t) = 0.
        restTtiss2(t) = KTrp2 * Trpbo_perml * Vtissue2
        actTtiss2(t) = 0.
        act0Ttiss2(t) = 0.
        restTtumor(t) = KTrptumor * Trpbo_perml * Init_Vtumor
        actTtumor(t) = 0.
        act0Ttumor(t) = 0.
        ## B cells 
        Bpb(t) = Bpbo_perml * Vpb
        B1920tiss3(t) = KBp3 * Bpbref_perml * Vtissue3
        B19no20tiss3(t) = KBp3 * Bpbref_perml * Vtissue3 * B19no20_B1920_ratio
        ## Tumor
        Btumor(t) = Init_Vtumor*0.375/Vc
        ## IL6
        IL6pb(t) = 0.
        IL6tiss2(t) = 0.
        IL6tumor(t) = 0.
        ## PK 
        TDBc_ugperml(t) = 0.
        TDBp_ugperml(t) = 0.
        TDBdepot_ug(t) = 0.
        TCEinjection_effect(t) = 0.
    end

    Vtumor = Vc * max(Btumor, 1.) / 0.375

    # PK 
    TDBt2_ugperml = Kp2*max(TDBc_ugperml, 0.)
    TDBtumor_ugperml = tumor_on*Kptumor*max(TDBc_ugperml, 0.)
    # T cells
    restTpb_perml = restTpb/Vpb
    act0Tpb_perml = act0Tpb/Vpb
    actTpb_perml = actTpb/Vpb
    totTpb_perml = restTpb_perml + act0Tpb_perml + actTpb_perml
    # B cells
    Bpb_perml = Bpb/Vpb
    Btumor_perml = Btumor/Vtumor
    # B:T ratio
    BTrratio_pb = Bpb/max(restTpb+act0Tpb,1)
    BTrratio_tumor = Btumor/max(restTtumor+act0Ttumor,1)
    TaBratio_pb = actTpb/max(Bpb,1)
    TaBratio_tumor = actTtumor/max(Btumor,1)
    # TDB-induced T cell activation 
    drugTpbact = VmT*(BTrratio_pb^S/(KmBT_act^S+BTrratio_pb^S))*((max(TDBc_ugperml, 0.)*1000)^ndrugactT/(KdrugactT^ndrugactT+(max(TDBc_ugperml, 0.)*1000)^ndrugactT))
    drugTtumoract = fTact*VmT*(max(BTrratio_tumor, 0)^S/(KmBT_act^S+max(BTrratio_tumor, 0)^S))*((TDBtumor_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBtumor_ugperml*1000)^ndrugactT))
    # TDB-induced B cell killing 
    drugBpbkill = VmB*max(0, TaBratio_pb)^nkill/(max(0, KmTB_kill)^nkill+max(0, TaBratio_pb)^nkill)*((max(TDBc_ugperml, 0.)*1000)/(KdrugB+(max(TDBc_ugperml, 0.)*1000)))
    drugBtumorkill = VmB*max(0, TaBratio_tumor)^nkill/(max(0, KmTB_kill*fKmTB_kill)^nkill+max(0, TaBratio_tumor)^nkill)*((TDBtumor_ugperml*1000)/(KdrugB+(TDBtumor_ugperml*1000)))
    # T cell activation 
    T_act_rest_PB = (kTact*(drugTpbact*restTpb-fTadeact*actTpb)); 
    T_act_rest_Tumor = (kTact*(drugTtumoract*restTtumor-fTadeact*actTtumor)); 
    T_act_postact_PB = (act0on*kTact*(drugTpbact*act0Tpb-fTadeact*actTpb)); 
    T_act_postact_tumor = (act0on*kTact*(drugTtumoract*act0Ttumor-fTadeact*actTtumor));
    # T cell conversion from post-activated to resting state
    T_postact_rest_conv_PB = (act0on*fTa0deact*act0Tpb)
    T_postact_rest_conv_LN = (act0on*fTa0deact*act0Ttiss2); 
    T_postact_rest_conv_Tumor = (act0on*fTa0deact*act0Ttumor);
    # T cell trafficking 
    rest_T_enter_LN = (tissue2on*kTrexit*Vpb*((1+TCEinjection_effect)*restTpb/Vpb*KTrp2-restTtiss2/Vtissue2));
    rest_T_enter_Tumor = (tumor_on*kTrexit*Vpb*((1+TCEinjection_effect)*restTpb/Vpb*KTrptumor-restTtumor/Vtumor));
    T_act_traffick_LN = (tissue2on*kTaexit*Vpb*(actTpb/Vpb*KTrp2*fTap-actTtiss2/Vtissue2));
    T_act_traffick_Tumor = (tumor_on*kTaexit*Vpb*(actTpb/Vpb*KTrptumor*fTap-actTtumor/Vtumor));
    postact_T_traffick_LN = (tissue2on*act0on*kTrexit*Vpb*((1+TCEinjection_effect)*act0Tpb/Vpb*KTrp2-act0Ttiss2/Vtissue2));
    postact_T_traffick_Tumor = (tumor_on*act0on*kTrexit*Vpb*((1+TCEinjection_effect)*act0Tpb/Vpb*KTrptumor-act0Ttumor/Vtumor));
    # B cell trafficking 
    B_leaving_BM = tissue3on*kBtiss3exit*Vpb*max(0,B1920tiss3/Vtissue3-Bpb/Vpb*KBp3);
    # diff equations
    eqs = [
        # T cells
        Dt(restTpb) ~ ( T_postact_rest_conv_PB - rest_T_enter_LN - rest_T_enter_Tumor - T_act_rest_PB + 
                        (kTgen*Vpb*Trpbref_perml*(max(0,1-totTpb_perml/Trpbref_perml)^2)) - fTrapop*kTaapop*max(0,restTpb-Trpbref_perml*Vpb)  ), 
        Dt(actTpb) ~ ( T_act_rest_PB + T_act_postact_PB - T_act_traffick_LN - T_act_traffick_Tumor - (kTaapop*(actTpb+fAICD*actTpb^2/(Vpb*Trpbref_perml))) ), 
        Dt(act0Tpb) ~ ( (act0on*fTaprolif*kTprolif*actTpb) - (fTa0apop*kTaapop*act0Tpb) - T_postact_rest_conv_PB - T_act_postact_PB
                          - postact_T_traffick_LN - postact_T_traffick_Tumor),
        Dt(act0Ttiss2) ~ ((act0on*fTaprolif*kTprolif*actTtiss2) + postact_T_traffick_LN - T_postact_rest_conv_LN - (fTa0apop*kTaapop*act0Ttiss2)), 
        Dt(restTtiss2) ~ (T_postact_rest_conv_LN + rest_T_enter_LN ), 
        Dt(actTtiss2) ~ ( T_act_traffick_LN ), 
        Dt(restTtumor) ~ (T_postact_rest_conv_Tumor - T_act_rest_Tumor + rest_T_enter_Tumor), 
        Dt(actTtumor) ~ (T_act_postact_tumor + T_act_rest_Tumor  + T_act_traffick_Tumor), 
        Dt(act0Ttumor) ~ ( 0 - T_postact_rest_conv_Tumor - T_act_postact_tumor + (act0on*fTaprolif*kTprolif*actTtumor) + postact_T_traffick_Tumor - (fTa0apop*kTaapop*act0Ttumor) ), 
        # PK 
        Dt(TDBc_ugperml) ~ ( kabs_TDB * TDBdepot_ug/V1_TDB - 1/V1_TDB * Q_TDB * (TDBc_ugperml - TDBp_ugperml) - CL_TDB/V1_TDB * TDBc_ugperml - Vm_tdb/V1_TDB / (Km_tdb + TDBc_ugperml) * TDBc_ugperml ), 
        Dt(TDBp_ugperml) ~ ( 1/V2_TDB * Q_TDB * (TDBc_ugperml - TDBp_ugperml) ), 
        Dt(TDBdepot_ug) ~ ( - kabs_TDB * TDBdepot_ug ), 
        # IL6 
        Dt(IL6pb) ~ (
            kIL6prod*actTpb/Vpb*(Bpb_perml/Bpbref_perml)*(max(TDBc_ugperml, 0.)*1000)^ndrugactT/(KdrugactT^ndrugactT+(max(TDBc_ugperml, 0.)*1000)^ndrugactT) - 
            log(2)/(thalfIL6/DayToHour)*IL6pb
            ), 
        Dt(IL6tumor) ~ (
            kIL6prod*actTtumor/Vtumor*(Btumor_perml/(Bpbref_perml*KBptumor))*((TDBtumor_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBtumor_ugperml*1000)^ndrugactT)) - 
            log(2)/(thalfIL6/DayToHour)*IL6tumor
            ), 
        # B cells 
        Dt(Bpb) ~ ( tissue3on*kBapop*Bpb + tissue3on*kBprolif*Bpbref_perml*Vpb*max(0,1-Bpb/(Bpbref_perml*Vpb)) - 
                    kBapop*Bpb - kBkill*drugBpbkill*Bpb + B_leaving_BM ), 
        Dt(B1920tiss3) ~ ( kBprolif*KBp3*Bpbref_perml*Vtissue3*max(0,1-B1920tiss3/(KBp3*Bpbref_perml*Vtissue3)) - B_leaving_BM + 
            kBapop/kBmat_kBapop_ratio*B19no20tiss3 - kBapop*B1920tiss3 ), 
            Dt(B19no20tiss3) ~ ( (Bpbref_perml*KBp3*B19no20_B1920_ratio*Vtissue3*(kBapop+(kBapop/kBmat_kBapop_ratio))*(1+5*max(0,1-Bpb/(Bpbref_perml*Vpb)))) - ((kBapop/kBmat_kBapop_ratio)*B19no20tiss3) - (kBapop*B19no20tiss3) + (kBprolif*KBp3*Bpbref_perml*B19no20_B1920_ratio*Vtissue3*(max(0,1-B19no20tiss3/(KBp3*Bpbref_perml*B19no20_B1920_ratio*Vtissue3)))^1)),
        Dt(Btumor) ~ ( (kBtumorprolif*Btumor) - (fBkill*kBkill*drugBtumorkill*Btumor ) ), 
        # other 
        Dt(TCEinjection_effect) ~ (-(log(2)/tinjhalf*TCEinjection_effect))
    ]; 
    ODESystem(eqs, t, vars, pars; name=name)
end

