# date: 4/10/24
# author: Yuezhe Li 
# purpose of this code: to breakdown Hosseini model in multiple tissue compartment

using Pkg; Pkg.activate("")

using DifferentialEquations, ModelingToolkit, Catalyst 
using DataFrames, DataFramesMeta, CSV
using AlgebraOfGraphics, CairoMakie

const DayToHour = 24.
const HourToMinute = 60.

# define tumor cell volume
const Vc = 4/3 * pi * (4E-6)^3 * 1E6; # [mL]

# PK module 
function PK(; name)
    pars = @parameters begin 
        Vm_tdb = 212                    # [μg/d]
        Km_tdb = 4.74                   # [μg/d]
        kabs_TDB = 0.26
        CL_mosun = 28.9                 # [mL/day]
        Q_mosun = 82.187                # [mL/day]
        V1_mosun = 125.052              # [mL]
        V2_mosun = 589.798              # [mL]
    end
    @variables t; 
    vars = @variables begin 
        TDBc_ugperml(t) = 0.
        TDBp_ugperml(t) = 0.
        TDBdepot_ug(t) = 0.
    end
    Dt = Differential(t); 
    eqs = [
        Dt(TDBc_ugperml) ~ ( kabs_TDB * TDBdepot_ug/V1_mosun - 1/V1_mosun * Q_mosun * (TDBc_ugperml - TDBp_ugperml) - CL_mosun/V1_mosun * TDBc_ugperml - Vm_tdb/V1_mosun / (Km_tdb + TDBc_ugperml) * TDBc_ugperml ), 
        Dt(TDBp_ugperml) ~ ( 1/V2_mosun * Q_mosun * (TDBc_ugperml - TDBp_ugperml) ), 
        Dt(TDBdepot_ug) ~ ( - kabs_TDB * TDBdepot_ug ),
    ];
    ODESystem(eqs, t, vars, pars; name=name);
end

@named pk = PK();

# blood cell module 
function blood_cells(; name)
    pars = @parameters begin 
        ## B cells 
        kBapop = 0.02                   # Apoptosis rate of B-cells [1/d]
        kBkill = 275.19                 # killing rate of B-cells [1/d]
        kBprolif = 0.05                 # proliferation rate of B cells [1/d]
        fBexit = 1.
        kBtiss3exit = 0.03
        B19no20_B1920_ratio = 0.25      # CD19+CD20- B-cell:CD19+CD20+ B-cell ratio in BM
        fBtissue3_v1 = 1.
        kBapop_cll = 1.
        ## TDB : B-cell kill
        VmB = 0.95                      # maximum rate of B-cell killing 
        KmTB_kill = 0.75                # the Tact/B ratio that the rate if half of VmB
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
        fa0 = 0.                        # Ratio of the placebo/TDB injection effect on traficking of activated CD8+ T-cells from PB to tissues to that of resting CD8+ T-cells
        fTa0deact = 0.                  # Conversion rate of post-activated CD8+ T-cells to resting CD8+ T-cells [1/d]
        finj = 1.                       # The placebo injection effect on traficking of resting or "exhausted" CD8+ T-cells from PB to tissues
        kTgen = 1.
        ## T cells : TDB activation 
        VmT = 0.9                       # maximum rate of T-cell activation [unitless]
        KmBT_act = 0.72                 # The B/T rario, at which the rate is half of VmT (at high concentrations of TDB)
        ndrugactT = 0.8                 # Hill coefficient for the TDB portion of the Michaelis-Menten equation
        S = 1.4                         # Hill coefficient for the B:T portion of the Michaelis-Menten equation
        KdrugactT = 130.08              # the TDB conc at which T-cell activation was half of VmT [ng/mL]
        ## IL6 
        thalfIL6 = 20.                  # [minute]
        ## other 
        tinjhalf = 1.
        Trpbo_perml = 2.0E6             # [1/mL]
        Bpbo_perml = 1.0E6              # [1/mL]
        act0on = 1.
        tissue1on = 1.
        tissue2on = 1.
        tissue3on = 1.
        Bcell_tumor_trafficking_on = 1.
        Trpbref_perml = 2.0e6           # baseline conc of resting CD8+ T-cells in PB [cells/mL]
        Bpbref_perml = 1.0e6            # baseline conc of CD19+CD20+ B-cells in PB [cells/mL]
        # tumor 
        tumor_on = 0.
        Vtumor = 1.0                    # [mL] 
        # params that changed between cyno and human
        # kIL6prod = 14.0               # production rate of IL6 [pg/d] (homo)
        kIL6prod = 4.0                  # production rate of IL6 [pg/d] (cyno)
        Vpb = 380.                     # physiological volume of peripheral blood [mL]
        Vtissue = 7.                  # physiological volume of spleen [mL]
        Vtissue2 = 25.                 # total physiological volume of LNs [mL]
        Vtissue3 = 50.0                # Physiology volume of BM [mL]
        # Vpb = 5000.                     # physiological volume of peripheral blood [mL]
        # Vtissue = 210.                  # physiological volume of spleen [mL]
        # Vtissue2 = 400.                 # total physiological volume of LNs [mL]
        # Vtissue3 = 500.0                # Physiology volume of BM [mL]
        KTrp = 500.0                    # Partition coeff: ratio of the baseline conc. of resting CD8+ T-cells in spleen to that in PB
        KTrp2 = 500.
        KTrp3 = 50.0
        KTrptumor = 1.
        # KTrp = 200.0                    # Partition coeff: ratio of the baseline conc. of resting CD8+ T-cells in spleen to that in PB
        # KTrp2 = 190.
        # KTrp3 = 60.0
        # KTrptumor = 125.
        KBp = 900.                      # Partition coeff: ratio of the baseline conc. of B-cells in spleen to that in PB
        KBp2 = 600.
        KBp3 = 60.0
        KBptumor = 1.
        # KBp = 333.                      # Partition coeff: ratio of the baseline conc. of B-cells in spleen to that in PB
        # KBp2 = 190.
        # KBp3 = 80.0
        # KBptumor = 6.5E3
    end
    @variables t; 
    Dt = Differential(t); 
    vars = @variables begin 
        TDBc_ugperml(t) = 0.
        injection_effect(t) = 0.
        IL6pb(t) = 0.
        ## T cells
        restTpb(t) = Trpbo_perml * Vpb
        actTpb(t) = 0. 
        act0Tpb(t) = 0.
        restTtiss(t) = KTrp * Trpbo_perml * Vtissue
        actTtiss(t) = 0.
        act0Ttiss(t) = 0.
        restTtiss2(t) = KTrp2 * Trpbo_perml * Vtissue2
        actTtiss2(t) = 0.
        act0Ttiss2(t) = 0.
        restTtiss3(t) = KTrp3 * Trpbo_perml * Vtissue3
        actTtiss3(t) = 0.
        act0Ttiss3(t) = 0.
        restTtumor(t) = KTrptumor * Trpbo_perml * Vtumor
        actTtumor(t) = 0.
        act0Ttumor(t) = 0.
        ## B cells 
        Bpb(t) = Bpbo_perml * Vpb
        Btiss(t) = KBp * Bpbo_perml * Vtissue
        Btiss2(t) = KBp2 * Bpbo_perml * Vtissue2
        B1920tiss3(t) = KBp3 * Bpbref_perml * Vtissue3
        B19no20tiss3(t) = KBp3 * Bpbref_perml * Vtissue3 * B19no20_B1920_ratio
        Btumor(t) = Vtumor*0.375/Vc
    end
    # T cells
    restTpb_perml = restTpb/Vpb
    act0Tpb_perml = act0Tpb/Vpb
    actTpb_perml = actTpb/Vpb
    totTpb_perml = restTpb_perml + act0Tpb_perml + actTpb_perml
    # B cells
    Bpb_perml = Bpb/Vpb
    # B:T ratio
    BTrratio_pb = Bpb/max(restTpb+act0Tpb,1)
    TaBratio_pb = actTpb/max(Bpb,1)
    # TDB-induced T cell activation 
    drugTpbact = VmT*(BTrratio_pb^S/(KmBT_act^S+BTrratio_pb^S))*((max(TDBc_ugperml, 0.)*1000)^ndrugactT/(KdrugactT^ndrugactT+(max(TDBc_ugperml, 0.)*1000)^ndrugactT))
    # TDB-induced B cell killing 
    drugBpbkill = VmB*max(0, TaBratio_pb)^nkill/(max(0, KmTB_kill)^nkill+max(0, TaBratio_pb)^nkill)*((max(TDBc_ugperml, 0.)*1000)/(KdrugB+(max(TDBc_ugperml, 0.)*1000)))
    # T cell activation 
    T_act_rest_PB = (kTact*(drugTpbact*restTpb-fTadeact*actTpb)); 
    T_act_postact_PB = (act0on*kTact*(drugTpbact*act0Tpb-fTadeact*actTpb)); 
    # T cell conversion from post-activated to resting state
    T_postact_rest_conv_PB = (act0on*fTa0deact*act0Tpb)
    # T cell trafficking 
    rest_T_enter_spleen = (tissue1on*kTrexit*Vpb*((1+finj*injection_effect)*restTpb/Vpb*KTrp-restTtiss/Vtissue));
    rest_T_enter_LN = (tissue2on*kTrexit*Vpb*((1+finj*injection_effect)*restTpb/Vpb*KTrp2-restTtiss2/Vtissue2));
    rest_T_enter_BM = (tissue3on*kTrexit*Vpb*((1+finj*injection_effect)*restTpb/Vpb*KTrp3-restTtiss3/Vtissue3));
    rest_T_enter_Tumor = (tumor_on*kTrexit*Vpb*((1+finj*injection_effect)*restTpb/Vpb*KTrptumor-restTtumor/Vtumor));
    T_act_traffick_spleen = (tissue1on*kTaexit*Vpb*((1+fa0*finj*injection_effect)*actTpb/Vpb*KTrp*fTap-actTtiss/Vtissue)); 
    T_act_traffick_LN = (tissue2on*kTaexit*Vpb*((1+fa0*finj*injection_effect)*actTpb/Vpb*KTrp2*fTap-actTtiss2/Vtissue2));
    T_act_traffick_BM = (tissue3on*kTaexit*Vpb*((1+fa0*finj*injection_effect)*actTpb/Vpb*KTrp3*fTap-actTtiss3/Vtissue3));
    T_act_traffick_Tumor = (tumor_on*kTaexit*Vpb*((1+fa0*finj*injection_effect)*actTpb/Vpb*KTrptumor*fTap-actTtumor/Vtumor));
    postact_T_traffick_spleen = (tissue1on*act0on*kTrexit*Vpb*((1+finj*injection_effect)*act0Tpb/Vpb*KTrp-act0Ttiss/Vtissue)); 
    postact_T_traffick_LN = (tissue2on*act0on*kTrexit*Vpb*((1+finj*injection_effect)*act0Tpb/Vpb*KTrp2-act0Ttiss2/Vtissue2));
    postact_T_traffick_BM = (tissue3on*act0on*kTrexit*Vpb*((1+finj*injection_effect)*act0Tpb/Vpb*KTrp3-act0Ttiss3/Vtissue3)); 
    postact_T_traffick_Tumor = (tumor_on*act0on*kTrexit*Vpb*((1+finj*injection_effect)*act0Tpb/Vpb*KTrptumor-act0Ttumor/Vtumor));
    # B cell trafficking 
    B_enter_spleen = tissue1on*fBexit*kTaexit*Vpb*max(0,Bpb/Vpb*KBp-Btiss/Vtissue);
    B_enter_LN = tissue2on*fBexit*kTaexit*Vpb*max(0,Bpb/Vpb*KBp2-Btiss2/Vtissue2);
    B_leaving_BM = tissue3on*fBtissue3_v1*kBtiss3exit*Vpb*max(0,B1920tiss3/Vtissue3-Bpb/Vpb*KBp3);
    B_enter_Tumor = Bcell_tumor_trafficking_on*tumor_on*fBexit*kTaexit*Vpb*max(0,Bpb/Vpb*KBptumor-Btumor/Vtumor);
    # diff equations
    eqs = [
        # T cells
        Dt(restTpb) ~ ( T_postact_rest_conv_PB - rest_T_enter_spleen - rest_T_enter_LN - rest_T_enter_BM - rest_T_enter_Tumor - T_act_rest_PB + 
                        (kTgen*Vpb*Trpbref_perml*(max(0,1-totTpb_perml/Trpbref_perml)^2)) - fTrapop*kTaapop*max(0,restTpb-Trpbref_perml*Vpb)  ), 
        Dt(actTpb) ~ ( T_act_rest_PB + T_act_postact_PB - T_act_traffick_spleen - T_act_traffick_LN - T_act_traffick_BM - T_act_traffick_Tumor - (kTaapop*(actTpb+fAICD*actTpb^2/(Vpb*Trpbref_perml))) ), 
        Dt(act0Tpb) ~ ( (act0on*fTaprolif*kTprolif*actTpb) - (fTa0apop*kTaapop*act0Tpb) - T_postact_rest_conv_PB - T_act_postact_PB
                         - postact_T_traffick_spleen - postact_T_traffick_LN - postact_T_traffick_BM - postact_T_traffick_Tumor),
        # IL6 
        Dt(IL6pb) ~ (
            kIL6prod*actTpb/Vpb*(Bpb_perml/Bpbref_perml)*(max(TDBc_ugperml, 0.)*1000)^ndrugactT/(KdrugactT^ndrugactT+(max(TDBc_ugperml, 0.)*1000)^ndrugactT) - 
            log(2)/(thalfIL6/HourToMinute/DayToHour)*IL6pb
            ), 
        # B cells  
        Dt(Bpb) ~ ( tissue3on*fBtissue3_v1*kBapop*Bpb*kBapop_cll + tissue3on*fBtissue3_v1*kBprolif*kBapop_cll*Bpbref_perml*Vpb*max(0,1-Bpb/(Bpbref_perml*Vpb)) - 
            kBapop*Bpb*kBapop_cll - kBkill*drugBpbkill*Bpb - B_enter_spleen - B_enter_LN + B_leaving_BM - B_enter_Tumor ), 
        # other 
        Dt(injection_effect) ~ (-(log(2)/tinjhalf*injection_effect))
    ];
    ODESystem(eqs, t, vars, pars; name=name);
end

@named blood_b_t = blood_cells();

# tumor module
function tumor(; name)
    pars1 = @parameters begin 
        # things unique to the tissue 
        fBkill = 0.1                    # kBkill is multiplied by this factor (<1) to capture the less efficient killing of tissue B-cells
        fKmTB_kill = 10.0               # KmTB_kill is multiplied by this factor (>1) to compensate for the fact the tissues are not well-stirred environments like PB 
        fTact = 0.25                    # VmT is multiplied by this factor (<1) to capture the less efficient activation of tissue T-cells
        # tumor 
        tumor_on = 0.
        Vtumor = 1.0                    # [mL] 
        kBtumorprolif = 0.025           # proliferation rate of B-cells in tumor [1/d]
        Bcell_tumor_trafficking_on = 1.
        KTrptumor = 1.
        KBptumor = 1.
        Kptumor = 0.05
    end
    pars_b = @parameters kBkill, fBexit 
    pars_t = @parameters kTprolif, kTaexit, kTact, fTadeact, fTap, kTaapop, fTaprolif, kTrexit, fTa0apop, fa0, fTa0deact, finj
    pars_drug = @parameters VmB, KmTB_kill, nkill, KdrugB, VmT, KmBT_act, S, ndrugactT, KdrugactT
    pars_other = @parameters thalfIL6, Trpbo_perml, Bpbo_perml, act0on, kIL6prod, Vpb, Bpbref_perml
    pars = vcat(pars1, pars_b, pars_t, pars_drug, pars_other);
    
    @variables t; 
    Dt = Differential(t); 

    vars_other = @variables restTpb(t) actTpb(t) act0Tpb(t) Bpb(t) TDBc_ugperml(t) injection_effect(t) 
    vars_tumor = @variables begin 
        restTtumor(t) = KTrptumor * Trpbo_perml * Vtumor
        actTtumor(t) = 0.
        act0Ttumor(t) = 0.
        Btumor(t) = Vtumor*0.375/Vc
        IL6tumor(t) = 0.
    end
    vars = vcat(vars_other, vars_tumor);

    TDBtumor_ugperml = tumor_on*Kptumor*max(TDBc_ugperml, 0.);
    Btumor_perml = Btumor/Vtumor;
    BTrratio_tumor = Btumor/max(restTtumor+act0Ttumor,1)
    TaBratio_tumor = actTtumor/max(Btumor,1)
    drugTtumoract = fTact*VmT*(max(BTrratio_tumor, 0)^S/(KmBT_act^S+max(BTrratio_tumor, 0)^S))*((TDBtumor_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBtumor_ugperml*1000)^ndrugactT))
    drugBtumorkill = VmB*max(0, TaBratio_tumor)^nkill/(max(0, KmTB_kill*fKmTB_kill)^nkill+max(0, TaBratio_tumor)^nkill)*((TDBtumor_ugperml*1000)/(KdrugB+(TDBtumor_ugperml*1000)))
    T_act_rest_Tumor = (kTact*(drugTtumoract*restTtumor-fTadeact*actTtumor)); 
    T_act_postact_tumor = (act0on*kTact*(drugTtumoract*act0Ttumor-fTadeact*actTtumor));
    T_postact_rest_conv_Tumor = (act0on*fTa0deact*act0Ttumor);
    rest_T_enter_Tumor = (tumor_on*kTrexit*Vpb*((1+finj*injection_effect)*restTpb/Vpb*KTrptumor-restTtumor/Vtumor));
    T_act_traffick_Tumor = (tumor_on*kTaexit*Vpb*((1+fa0*finj*injection_effect)*actTpb/Vpb*KTrptumor*fTap-actTtumor/Vtumor));
    postact_T_traffick_Tumor = (tumor_on*act0on*kTrexit*Vpb*((1+finj*injection_effect)*act0Tpb/Vpb*KTrptumor-act0Ttumor/Vtumor));
    B_enter_Tumor = Bcell_tumor_trafficking_on*tumor_on*fBexit*kTaexit*Vpb*max(0,Bpb/Vpb*KBptumor-Btumor/Vtumor);

    eqs = [
        Dt(restTtumor) ~ (T_postact_rest_conv_Tumor - T_act_rest_Tumor + rest_T_enter_Tumor), 
        Dt(actTtumor) ~ (T_act_postact_tumor + T_act_rest_Tumor  + T_act_traffick_Tumor), 
        Dt(act0Ttumor) ~ ( 0 - T_postact_rest_conv_Tumor - T_act_postact_tumor + (act0on*fTaprolif*kTprolif*actTtumor) + postact_T_traffick_Tumor - (fTa0apop*kTaapop*act0Ttumor) ), 
        Dt(IL6tumor) ~ (
            kIL6prod*actTtumor/Vtumor*(Btumor_perml/(Bpbref_perml*KBptumor))*((TDBtumor_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBtumor_ugperml*1000)^ndrugactT)) - 
            log(2)/(thalfIL6/HourToMinute/DayToHour)*IL6tumor
            ), 
        Dt(Btumor) ~ ( (kBtumorprolif*Btumor) - (fBkill*kBkill*drugBtumorkill*Btumor) + B_enter_Tumor ), 
    ];
    ODESystem(eqs, t, vars, pars; name=name);
end

@named tm = tumor();

# bone marrow module
function bone_marrow(; name)
    pars1 = @parameters begin 
        B19no20_B1920_ratio = 0.25      # CD19+CD20- B-cell:CD19+CD20+ B-cell ratio in BM
        kBmat_kBapop_ratio = 0.25       # The ratio of maturation rate of CD19+CD20- B-cells in BM to B-cell apoptosis rate
        Kp3 = 0.07
    end
    pars_b = @parameters kBkill, fBexit, KBp3, kBtiss3exit, fBtissue3_v1, fBkill, fKmTB_kill, kBprolif, kBapop, kBapop_cll, KBp, KBp2
    pars_t = @parameters kTprolif, kTaexit, kTact, fTadeact, fTap, kTaapop, fTaprolif, kTrexit, fTa0apop, fa0, fTa0deact, finj, KTrp3, fTact
    pars_drug = @parameters VmB, KmTB_kill, nkill, KdrugB, VmT, KmBT_act, S, ndrugactT, KdrugactT
    pars_other = @parameters thalfIL6, Trpbo_perml, Bpbo_perml, act0on, kIL6prod, Vpb, Bpbref_perml, Vtissue, Vtissue2, Vtissue3, tissue3on
    pars = vcat(pars1, pars_b, pars_t, pars_drug, pars_other);
    @variables t; 
    vars_other = @variables restTpb(t) actTpb(t) act0Tpb(t) Bpb(t) Btiss(t) Btiss2(t) TDBc_ugperml(t) injection_effect(t)
    vars_bm = @variables begin 
        restTtiss3(t) = KTrp3 * Trpbo_perml * Vtissue3
        actTtiss3(t) = 0.
        act0Ttiss3(t) = 0.
        B1920tiss3(t) = KBp3 * Bpbref_perml * Vtissue3
        B19no20tiss3(t) = KBp3 * Bpbref_perml * Vtissue3 * B19no20_B1920_ratio
        IL6tiss3(t) = 0.
    end
    vars = vcat(vars_other, vars_bm);
    Dt = Differential(t); 
    TDBt3_ugperml = tissue3on*Kp3*max(TDBc_ugperml, 0.);
    B1920tiss3_perml = B1920tiss3/Vtissue3;
    B1920Trratio_tiss3 = B1920tiss3/max(restTtiss3+act0Ttiss3,1);
    TaB1920ratio_tiss3 = actTtiss3/max(B1920tiss3,1);
    drugTtissact3 = fTact*VmT*(B1920Trratio_tiss3^S/(KmBT_act^S+B1920Trratio_tiss3^S))*((TDBt3_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBt3_ugperml*1000)^ndrugactT));
    drugBtisskill3 = VmB*max(0,TaB1920ratio_tiss3)^nkill/(max(0,KmTB_kill*fKmTB_kill)^nkill+max(0,TaB1920ratio_tiss3)^nkill)*((TDBt3_ugperml*1000)/(KdrugB+(TDBt3_ugperml*1000)));
    # T cell dynamics
    T_act_rest_BM = (kTact*(drugTtissact3*restTtiss3-fTadeact*actTtiss3)); 
    T_act_postact_BM = (act0on*kTact*(drugTtissact3*act0Ttiss3-fTadeact*actTtiss3)); 
    T_postact_rest_conv_BM = (act0on*fTa0deact*act0Ttiss3);
    # T cell trafficking 
    rest_T_enter_BM = (tissue3on*kTrexit*Vpb*((1+finj*injection_effect)*restTpb/Vpb*KTrp3-restTtiss3/Vtissue3));
    T_act_traffick_BM = (tissue3on*kTaexit*Vpb*((1+fa0*finj*injection_effect)*actTpb/Vpb*KTrp3*fTap-actTtiss3/Vtissue3));
    postact_T_traffick_BM = (tissue3on*act0on*kTrexit*Vpb*((1+finj*injection_effect)*act0Tpb/Vpb*KTrp3-act0Ttiss3/Vtissue3)); 
    # B cell trafficking 
    B_leaving_BM = tissue3on*fBtissue3_v1*kBtiss3exit*Vpb*max(0,B1920tiss3/Vtissue3-Bpb/Vpb*KBp3);
    eqs = [
        Dt(restTtiss3) ~ (rest_T_enter_BM + T_postact_rest_conv_BM - T_act_rest_BM ), 
        Dt(act0Ttiss3) ~ (postact_T_traffick_BM - (fTa0apop*kTaapop*act0Ttiss3) - T_act_postact_BM - T_postact_rest_conv_BM + (act0on*fTaprolif*kTprolif*actTtiss3)), 
        Dt(actTtiss3) ~ (T_act_traffick_BM + T_act_postact_BM + T_act_rest_BM ), 
        Dt(B1920tiss3) ~ ( kBprolif*KBp3*Bpbref_perml*Vtissue3*max(0,1-B1920tiss3/(KBp3*Bpbref_perml*Vtissue3)) - B_leaving_BM + 
            kBapop/kBmat_kBapop_ratio*B19no20tiss3 - kBapop*B1920tiss3*kBapop_cll - fBkill*kBkill*drugBtisskill3*B1920tiss3 ), 
        Dt(B19no20tiss3) ~ ((Bpbref_perml*KBp3*B19no20_B1920_ratio*Vtissue3*(kBapop+(kBapop/kBmat_kBapop_ratio))*(1+5*(max(0,1-(Bpb+Btiss+Btiss2)/(Bpbref_perml*Vpb+KBp*Bpbo_perml*Vtissue+KBp2*Bpbo_perml*Vtissue2)))^1)) - ((kBapop/kBmat_kBapop_ratio)*B19no20tiss3) - (kBapop*B19no20tiss3*kBapop_cll) + (kBprolif*KBp3*Bpbref_perml*B19no20_B1920_ratio*Vtissue3*(max(0,1-B19no20tiss3/(KBp3*Bpbref_perml*B19no20_B1920_ratio*Vtissue3)))^1)),
        Dt(IL6tiss3) ~ (
            kIL6prod*actTtiss3/Vtissue3*(  B1920tiss3_perml/(Bpbref_perml*KBp3)*(TDBt3_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBt3_ugperml*1000)^ndrugactT)  ) - 
            log(2)/(thalfIL6/HourToMinute/DayToHour)*IL6tiss3), 
    ];
    ODESystem(eqs, t, vars, pars; name=name)
end

@named bm = bone_marrow();

# lymph node module 
function lymph_nodes(; name)
    pars1 = @parameters Kp2 = 0.07 tissue2on KBp2 KTrp2 Vtissue2
    pars_b = @parameters kBkill, fBexit, kBtiss3exit, fBtissue3_v1, fBkill, fKmTB_kill, kBprolif, kBapop, kBapop_cll
    pars_t = @parameters kTprolif, kTaexit, kTact, fTadeact, fTap, kTaapop, fTaprolif, kTrexit, fTa0apop, fa0, fTa0deact, finj, fTact
    pars_drug = @parameters VmB, KmTB_kill, nkill, KdrugB, VmT, KmBT_act, S, ndrugactT, KdrugactT
    pars_other = @parameters thalfIL6, Trpbo_perml, Bpbo_perml, act0on, kIL6prod, Vpb, Bpbref_perml
    pars = vcat(pars1, pars_b, pars_t, pars_drug, pars_other);
    @variables t; 
    vars_other = @variables restTpb(t) actTpb(t) act0Tpb(t) Bpb(t) Btiss(t) Btiss2(t) TDBc_ugperml(t) injection_effect(t)
    vars_ln = @variables begin 
        restTtiss2(t) = KTrp2 * Trpbo_perml * Vtissue2
        actTtiss2(t) = 0.
        act0Ttiss2(t) = 0.
        Btiss2(t) = KBp2 * Bpbo_perml * Vtissue2
        IL6tiss2(t) = 0.
    end
    vars = vcat(vars_other, vars_ln);
    Dt = Differential(t); 
    TDBt2_ugperml = Kp2*max(TDBc_ugperml, 0.);
    Btiss2_perml = Btiss2/Vtissue2;
    BTrratio_tiss2 = Btiss2/max(restTtiss2+act0Ttiss2,1); 
    TaBratio_tiss2 = actTtiss2/max(Btiss2,1);
    drugTtissact2 = fTact*VmT*(BTrratio_tiss2^S/(KmBT_act^S+BTrratio_tiss2^S))*((TDBt2_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBt2_ugperml*1000)^ndrugactT));
    drugBtisskill2 = VmB*max(0, TaBratio_tiss2)^nkill/(max(0, KmTB_kill*fKmTB_kill)^nkill+max(0, TaBratio_tiss2)^nkill)*((TDBt2_ugperml*1000)/(KdrugB+(TDBt2_ugperml*1000)));
    T_act_rest_LN = (kTact*(drugTtissact2*restTtiss2-fTadeact*actTtiss2));
    T_act_postact_LN = (act0on*kTact*(drugTtissact2*act0Ttiss2-fTadeact*actTtiss2)); 
    T_postact_rest_conv_LN = (act0on*fTa0deact*act0Ttiss2); 
    # T/B cell trafficking 
    rest_T_enter_LN = (tissue2on*kTrexit*Vpb*((1+finj*injection_effect)*restTpb/Vpb*KTrp2-restTtiss2/Vtissue2));
    T_act_traffick_LN = (tissue2on*kTaexit*Vpb*((1+fa0*finj*injection_effect)*actTpb/Vpb*KTrp2*fTap-actTtiss2/Vtissue2));
    postact_T_traffick_LN = (tissue2on*act0on*kTrexit*Vpb*((1+finj*injection_effect)*act0Tpb/Vpb*KTrp2-act0Ttiss2/Vtissue2));
    B_enter_LN = tissue2on*fBexit*kTaexit*Vpb*max(0,Bpb/Vpb*KBp2-Btiss2/Vtissue2);
    eqs = [
        Dt(act0Ttiss2) ~ ((act0on*fTaprolif*kTprolif*actTtiss2) + postact_T_traffick_LN - T_postact_rest_conv_LN - T_act_postact_LN - (fTa0apop*kTaapop*act0Ttiss2)), 
        Dt(restTtiss2) ~ (T_postact_rest_conv_LN + rest_T_enter_LN - T_act_rest_LN ), 
        Dt(actTtiss2) ~ ( T_act_postact_LN + T_act_traffick_LN + T_act_rest_LN ), 
        Dt(Btiss2) ~ ( B_enter_LN -(fBkill*kBkill*drugBtisskill2*Btiss2) + (kBprolif*KBp2*Bpbref_perml*Vtissue2*(max(0,1-Btiss2/(KBp2*Bpbref_perml*Vtissue2)))^1) ), 
        Dt(IL6tiss2) ~ (
            kIL6prod*actTtiss2/Vtissue2*(Btiss2_perml/(Bpbref_perml*KBp2))*((TDBt2_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBt2_ugperml*1000)^ndrugactT)) - 
            log(2)/(thalfIL6/HourToMinute/DayToHour)*IL6tiss2
            ), 
    ];  
    ODESystem(eqs, t, vars, pars; name=name)
end

@named ln = lymph_nodes();

# spleen module 
function Spleen(; name)
    pars1 = @parameters Kp = 0.14 tissue1on KBp KTrp Vtissue
    pars_b = @parameters kBkill, fBexit, kBtiss3exit, fBtissue3_v1, fBkill, fKmTB_kill, kBprolif, kBapop, kBapop_cll
    pars_t = @parameters kTprolif, kTaexit, kTact, fTadeact, fTap, kTaapop, fTaprolif, kTrexit, fTa0apop, fa0, fTa0deact, finj, fTact
    pars_drug = @parameters VmB, KmTB_kill, nkill, KdrugB, VmT, KmBT_act, S, ndrugactT, KdrugactT
    pars_other = @parameters thalfIL6, Trpbo_perml, Bpbo_perml, act0on, kIL6prod, Vpb, Bpbref_perml
    pars = vcat(pars1, pars_b, pars_t, pars_drug, pars_other);
    @variables t; 
    vars_other = @variables restTpb(t) actTpb(t) act0Tpb(t) Bpb(t) TDBc_ugperml(t) injection_effect(t)
    vars_sp = @variables begin 
        restTtiss(t) = KTrp * Trpbo_perml * Vtissue
        actTtiss(t) = 0.
        act0Ttiss(t) = 0.
        Btiss(t) = KBp * Bpbo_perml * Vtissue
        IL6tiss(t) = 0.
    end
    vars = vcat(vars_other, vars_sp);
    Dt = Differential(t); 
    TDBt_ugperml = Kp*max(TDBc_ugperml, 0.); 
    Btiss_perml = Btiss/Vtissue;
    BTrratio_tiss = Btiss/max(restTtiss+act0Ttiss,1);
    TaBratio_tiss = actTtiss/max(Btiss,1);
    drugTtissact = fTact*VmT*(BTrratio_tiss^S/(KmBT_act^S+BTrratio_tiss^S))*((TDBt_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBt_ugperml*1000)^ndrugactT));
    drugBtisskill = VmB*max(0, TaBratio_tiss)^nkill/(max(0, KmTB_kill*fKmTB_kill)^nkill+max(0, TaBratio_tiss)^nkill)*((TDBt_ugperml*1000)/(KdrugB+(TDBt_ugperml*1000)));
    T_act_rest_spleen = kTact*(drugTtissact*restTtiss-fTadeact*actTtiss);
    T_act_postact_spleen = (act0on*kTact*(drugTtissact*act0Ttiss-fTadeact*actTtiss));
    T_postact_rest_conv_spleen = (act0on*fTa0deact*act0Ttiss); 
    # T/B cell trafficking 
    rest_T_enter_spleen = (tissue1on*kTrexit*Vpb*((1+finj*injection_effect)*restTpb/Vpb*KTrp-restTtiss/Vtissue));
    T_act_traffick_spleen = (tissue1on*kTaexit*Vpb*((1+fa0*finj*injection_effect)*actTpb/Vpb*KTrp*fTap-actTtiss/Vtissue)); 
    postact_T_traffick_spleen = (tissue1on*act0on*kTrexit*Vpb*((1+finj*injection_effect)*act0Tpb/Vpb*KTrp-act0Ttiss/Vtissue)); 
    B_enter_spleen = tissue1on*fBexit*kTaexit*Vpb*max(0,Bpb/Vpb*KBp-Btiss/Vtissue);
    eqs = [
        Dt(restTtiss) ~ ( 0 - T_act_rest_spleen + rest_T_enter_spleen + T_postact_rest_conv_spleen ), 
        Dt(actTtiss) ~ ( T_act_rest_spleen + T_act_traffick_spleen + T_act_postact_spleen ), 
        Dt(act0Ttiss) ~ ( -(fTa0apop*kTaapop*act0Ttiss) - T_act_postact_spleen - T_postact_rest_conv_spleen + postact_T_traffick_spleen + (act0on*fTaprolif*kTprolif*actTtiss)), 
        Dt(Btiss) ~ ( B_enter_spleen - (fBkill*kBkill*drugBtisskill*Btiss) + (kBprolif*KBp*Bpbref_perml*Vtissue*(max(0,1-Btiss/(KBp*Bpbref_perml*Vtissue)))^1) ), 
        Dt(IL6tiss) ~ (
            kIL6prod*actTtiss/Vtissue*(Btiss_perml/(Bpbref_perml*KBp))*(TDBt_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBt_ugperml*1000)^ndrugactT) - 
            log(2)/(thalfIL6/HourToMinute/DayToHour)*IL6tiss
            ), 
    ];
    ODESystem(eqs, t, vars, pars; name=name)
end

@named sp = Spleen();

# join models together 
@variables t;

pk_b = extend(pk, blood_b_t);
pk_b_sp = extend(pk_b, sp);
pk_b_sp_ln = extend(pk_b_sp, ln);
pk_b_sp_ln_bm = extend(pk_b_sp_ln, bm);
pk_b_sp_ln_bm_tumor = extend(pk_b_sp_ln_bm, tm); 

full = structural_simplify(pk_b_sp_ln_bm_tumor);
# states(full)
# parameters(full)

prob = ODEProblem(full, [], (-1E-12, 80.));

# add dosing 
injection_effect_init = 10.0;
TDBc_ugperkg = 10.  # [ug]
subject_dose_time = 7. *[0 1 2 3];
tspan = (0., 80.);
subject_dose_amt = 3.4 * [1, 0.6, 1e-9, 0.025] * TDBc_ugperkg; # [ug]

cbs = [];
if length(subject_dose_time) > 1
    for i in 1:length(subject_dose_time)
        function affect!(integrator)
            integrator.u[28] += subject_dose_amt[i]/ 125.052;
            integrator.u[27] += injection_effect_init;
        end
        cb = PresetTimeCallback(subject_dose_time[i],affect!);
        global cbs = push!(cbs, cb);
    end
end
cbset = CallbackSet(cbs...);

sol = solve(prob, reltol=1e-18, saveat = 0.2, callback = cbset);
sdf = DataFrame(sol);
@select!(sdf, :timestamp, :restTpb, :actTpb, :Bpb, :IL6pb, :TDBc_ugperml, :injection_effect); 
sdf.Source .= "new";
sdf_long = stack(sdf, Not([:Source,:timestamp]));

# read in sims from the other simulation 
sims_old = CSV.read("../cyno/old-sims.csv", DataFrame, header=true);
@select!(sims_old, :timestamp, :restTpb, :actTpb, :Bpb, :IL6pb, :TDBc_ugperml, :injection_effect); 
sims_old.Source .= "old";
sdf_old = stack(sims_old, Not([:Source,:timestamp]));

# visual comparison 
df_long = vcat(sdf_old, sdf_long);

layer = visual(Lines,linewidth=4.0, alpha = 0.5);
xy = data(df_long) * mapping(:timestamp, :value, color=:Source,layout=:variable, linestyle = :Source);
# with_theme(theme_default())
fg = draw(layer * xy, axis=(; titlesize=16,xlabelsize=16,ylabelsize=16,xticklabelsize=16,yticklabelsize=16,xlabel="Time (days)", ylabel="Value (Mixed Units)"),
            figure=(; size=(1000, 500)), legend = (; labelsize=16,titlesize=16), facet = (; linkyaxes = :none))
                      
save("../figure/validation-new2.png", fg);

#---------------------------------- cyno to human ----------------------------------#

mosu2 = CSV.read("../data/lunsumio-pk-genetech.csv",DataFrame);

p_human = deepcopy(prob.p); 
p_human[findfirst(==("CL_mosun"), string.(parameters(full)))] = 637.5       # [mL/day]
p_human[findfirst(==("Q_mosun"), string.(parameters(full)))] = 1812.75      # [mL/day]
p_human[findfirst(==("V1_mosun"), string.(parameters(full)))] = 2758.5      # [mL]
p_human[findfirst(==("V2_mosun"), string.(parameters(full)))] = 13010.25    # [mL]
p_human[findfirst(==("Vpb"), string.(parameters(full)))] = 5000.            # [ml] blood volume; 
p_human[findfirst(==("Vtissue"), string.(parameters(full)))] = 210.         # [ml] spleen volume;
p_human[findfirst(==("Vtissue2"), string.(parameters(full)))] = 400.        # [ml] lymph nodes volume; 
p_human[findfirst(==("Vtissue3"), string.(parameters(full)))] = 500.        # [ml] bone marrow volume; 
p_human[findfirst(==("KTrp"), string.(parameters(full)))] = 200.
p_human[findfirst(==("KTrp2"), string.(parameters(full)))] = 190.
p_human[findfirst(==("KTrp3"), string.(parameters(full)))] = 60.
p_human[findfirst(==("KBp"), string.(parameters(full)))] = 333.
p_human[findfirst(==("KBp2"), string.(parameters(full)))] = 190.
p_human[findfirst(==("KBp3"), string.(parameters(full)))] = 80. 

dose_amt = [1., 2., 60., 60., 30.] .* 1e3;  # [ug]
dose_time = [0, 7, 14, 21, 42];                 # [day] 

u0_human = deepcopy(prob.u0);
u0_human[findfirst(==("injection_effect(t)"), string.(states(full)))] = injection_effect_init; 
u0_human[findfirst(==("TDBc_ugperml(t)"), string.(states(full)))] = dose_amt[1]/p_human[73]; 

cbs_human = [];
if length(dose_time) > 1
    for i in 2:length(dose_time)
        function affect!(integrator)
            integrator.u[findfirst(==("TDBc_ugperml(t)"), string.(states(full)))] += dose_amt[i]/ p_human[73];
            integrator.u[findfirst(==("injection_effect(t)"), string.(states(full)))] += injection_effect_init;
        end
        cb = PresetTimeCallback(dose_time[i],affect!);
        global cbs_human = push!(cbs_human, cb);
    end
end
cbset_human = CallbackSet(cbs_human...);

prob_human = remake(prob, u0 = u0_human, p = p_human, tspan = (0., 63.)); 

sol_human = solve(prob_human, reltol=1e-18, saveat = 0.2, callback = cbset_human);
sdf_human = DataFrame(sol_human);

using Plots
phigh = Plots.plot(title = "Mosunetuzumab PK (dose = 1,2,60,60,30mg)", titlefont = 8, legend = :outerright, xticks = ([0, 7, 14, 21, 42, 63], string.([0, 7, 14, 21, 42, 63])));
Plots.scatter!(mosu2.day, mosu2.mosun_ug_mL, label = "Lunsumio data from label", alpha = 0.5, color = :blue, markerstrokewidth=0, markersize=5);
Plots.plot!(sdf_human.timestamp, sdf_human.TDBc_ugperml, label = "sims", color = :blue);
Plots.xlabel!("time (day)", xguidefontsize = 8); Plots.ylabel!("Mosunetuzumab plasma conc (ug/mL)", yguidefontsize = 8);
Plots.plot!(yscale = :log10); Plots.ylims!(0.01, 100)

save("../figure/human-pk-validation.png", phigh);


#---------------------------------- alternative methods for param update ----------------------------------#

using SymbolicIndexingInterface: parameter_values, state_values
using SciMLStructures: Tunable, replace, replace!

# update parameters 
p2 = Dict([full.CL_mosun => 637.5, full.Q_mosun => 1812.75, full.V1_mosun => 2758.5, full.V2_mosun => 13010.25, 
    full.Vpb => 5000., full.Vtissue => 210., full.Vtissue2 => 400., full.Vtissue3 => 500., 
    full.KTrp => 200., full.KTrp2 => 190., full.KTrp3 => 60., full.KBp => 333., full.KBp2 => 190., full.KBp3 => 80.]); 
# p2[full.Vpb]

# update initial condition 
u0_2 = Dict([full.injection_effect => injection_effect_init, full.TDBc_ugperml => dose_amt[1]/p2[full.V1_mosun]]);

prob2 = remake(prob, u0 = u0_2, p = p2, tspan = (0., 63.)); 

sdf2 = DataFrame(solve(prob2, reltol=1e-18, saveat = 0.2, callback = cbset_human));

using Plots
Plots.plot(title = "Mosunetuzumab PK (dose = 1,2,60,60,30mg)", titlefont = 8, legend = :outerright, xticks = ([0, 7, 14, 21, 42, 63], string.([0, 7, 14, 21, 42, 63])));
Plots.scatter!(mosu2.day, mosu2.mosun_ug_mL, label = "Lunsumio data from label", alpha = 0.5, color = :blue, markerstrokewidth=0, markersize=5);
Plots.plot!(sdf2.timestamp, sdf2.TDBc_ugperml, label = "sims", color = :blue);
Plots.xlabel!("time (day)", xguidefontsize = 8); Plots.ylabel!("Mosunetuzumab plasma conc (ug/mL)", yguidefontsize = 8);
Plots.plot!(yscale = :log10); Plots.ylims!(0.01, 100)
