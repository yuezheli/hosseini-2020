# author: Yuezhe LI
# date: 10/3/2023
# source of this code: ADC0304F (ref the posters here)

using ComponentArrays
using Parameters: @unpack


function bs_lonca_ode_3!(du, u, p, t)
    # unpack variables
    @unpack actTtiss, Btiss, actTpb, restTtiss, restTpb,
    Bpb, act0Tpb, act0Ttiss, injection_effect, Btiss2, act0Ttiss2, restTtiss2, 
    actTtiss2, B1920tiss3, B19no20tiss3, 
    restTtiss3, act0Ttiss3, actTtiss3, restTtumor, actTtumor, Btumor_cd19neg, Btumor_cd19neg_trans, Btumor, Btumor_trans, act0Ttumor, 
    IL6pb, IL6tiss, IL6tiss2, IL6tiss3, IL6tumor, TDBc_ugperml, TDBp_ugperml, TDBdepot_ug, LONCAc_ugperml, LONCAp_ugperml, TDBKill, LoncaKill = u

    # unpack parameters
    @unpack Kp, Kp2, Kp3, Kptumor, 
    # B cell
    VmB , KmTB_kill, KdrugB, kBapop, kBprolif, kBkill, KBp2, Bpbo_perml, fBexit, KBp, Bpbref_perml, kBtiss3exit, KBp3, B19no20_B1920_ratio, fBtissue3_v1, 
    kBmat_kBapop_ratio, kBapop_cll, fBkill, fKmTB_kill, KBptumor, kBtumorprolif, 
    # T cell
    kTprolif, KdrugactT, VmT, KmBT_act, kTaexit, kTact, fTadeact, fTap, kTaapop, KTrp, fTaprolif, kTgen, fTrapop, Trpbref_perml, KTrp2, fa0, fAICD, fTa0deact, fTa0apop, finj, 
    kTrexit, KTrp3, fTact, KTrptumor, S, nkill, ndrugactT, 
    # IL6 
    kIL6prod, thalfIL6, IL6_tiss_contribution, 
    # parameters undocumented 
    act0on, tissue1on, tissue2on, tissue3on, tumor_on, Bcell_tumor_trafficking_on, tinjhalf, Vtumor, 
    # parameters not used in the ODE 
    Trpbo_perml, 
    # params updated for human 
    Vpb, Vtissue, Vtissue2, Vtissue3, 
    # params for TDB PK
    Vm_tdb, Km_tdb, kabs_TDB, fbio_TDB, BW, CL_TDB, V1_TDB, V2_TDB, Q_TDB, infusion_TDB, 
    # params related to Lonca 
    CL_lonca, V1_lonca, V2_lonca, Q_lonca, infusion_lonca, lonca_kill, lonca_EC50, n_kill_lonca, k_trans, lonca_kill_neg, k_trans_neg, lonca_EC50_neg = p

    Btumor_cd19neg = max(0., Btumor_cd19neg);
    Btumor_cd19neg_trans = max(0., Btumor_cd19neg_trans);
    Btumor = max(1., Btumor);
    Btumor_trans = max(0., Btumor_trans);
    Btumor_total = Btumor_cd19neg + Btumor
    # note the tumor volume computation was updated based on Lonca PBPK model; the 0.375 is a tumor cellularity parameter
    # note in Hosseini model, tumor volume is a fixed number, which didn't make much sense
    Vtumor = Vc * (Btumor + Btumor_trans + Btumor_cd19neg + Btumor_cd19neg_trans) / 0.375
    
    ## Repeated Assignments
    TDBc_ugperml = max(TDBc_ugperml, 0.)
    TDBt_ugperml = Kp*TDBc_ugperml
    TDBt2_ugperml = Kp2*TDBc_ugperml
    TDBt3_ugperml = tissue3on*Kp3*TDBc_ugperml
    TDBtumor_ugperml = tumor_on*Kptumor*TDBc_ugperml
    restTpb_perml = restTpb/Vpb
    act0Tpb_perml = act0Tpb/Vpb
    actTpb_perml = actTpb/Vpb
    Bpb_perml = Bpb/Vpb
    Btiss_perml = Btiss/Vtissue
    Btiss2_perml = Btiss2/Vtissue2
    B1920tiss3_perml = B1920tiss3/Vtissue3
    Btumor_perml = Btumor_total/Vtumor
    BTrratio_pb = Bpb/max(restTpb+act0Tpb,1)
    BTrratio_tiss = Btiss/max(restTtiss+act0Ttiss,1)
    BTrratio_tiss2 = Btiss2/max(restTtiss2+act0Ttiss2,1)
    B1920Trratio_tiss3 = B1920tiss3/max(restTtiss3+act0Ttiss3,1)
    BTrratio_tumor = Btumor_total/max(restTtumor+act0Ttumor,1)
    TaBratio_pb = actTpb/max(Bpb,1)
    TaBratio_tiss = actTtiss/max(Btiss,1)
    TaBratio_tiss2 = actTtiss2/max(Btiss2,1)
    TaB1920ratio_tiss3 = actTtiss3/max(B1920tiss3,1)
    TaBratio_tumor = actTtumor/Btumor_total
    totTpb_perml = restTpb_perml+act0Tpb_perml + actTpb_perml


    # Lonca-associated term 
    LONCAc_ugperml = max(LONCAc_ugperml, 0.)
    lonca1_ugperml = tissue1on*Kp*LONCAc_ugperml
    lonca2_ugperml = tissue2on*Kp2*LONCAc_ugperml
    lonca3_ugperml = tissue3on*Kp3*LONCAc_ugperml
    loncatumor_ugperml = tumor_on*Kptumor*LONCAc_ugperml

    loncaBbp = lonca_kill * LONCAc_ugperml^n_kill_lonca/(LONCAc_ugperml^n_kill_lonca+lonca_EC50^n_kill_lonca)
    loncaBtissueact1 = lonca_kill * lonca1_ugperml^n_kill_lonca/(lonca1_ugperml^n_kill_lonca+lonca_EC50^n_kill_lonca)
    loncaBtissueact2 = lonca_kill * lonca2_ugperml^n_kill_lonca/(lonca2_ugperml^n_kill_lonca+lonca_EC50^n_kill_lonca)
    loncaBtissueact3 = lonca_kill * lonca3_ugperml^n_kill_lonca/(lonca3_ugperml^n_kill_lonca+lonca_EC50^n_kill_lonca)
    loncaBtumoract = lonca_kill * loncatumor_ugperml^n_kill_lonca/(loncatumor_ugperml^n_kill_lonca+lonca_EC50^n_kill_lonca)
    loncaBtumoractneg = lonca_kill_neg * loncatumor_ugperml^n_kill_lonca/(loncatumor_ugperml^n_kill_lonca+lonca_EC50_neg^n_kill_lonca)

    drugTpbact = VmT*(BTrratio_pb^S/(KmBT_act^S+BTrratio_pb^S))*((TDBc_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBc_ugperml*1000)^ndrugactT))
    drugTtissact = fTact*VmT*(BTrratio_tiss^S/(KmBT_act^S+BTrratio_tiss^S))*((TDBt_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBt_ugperml*1000)^ndrugactT))
    drugTtissact2 = fTact*VmT*(BTrratio_tiss2^S/(KmBT_act^S+BTrratio_tiss2^S))*((TDBt2_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBt2_ugperml*1000)^ndrugactT))
    drugTtissact3 = fTact*VmT*(B1920Trratio_tiss3^S/(KmBT_act^S+B1920Trratio_tiss3^S))*((TDBt3_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBt3_ugperml*1000)^ndrugactT))
    drugTtumoract = fTact*VmT*(max(BTrratio_tumor, 0)^S/(KmBT_act^S+max(BTrratio_tumor, 0)^S))*((TDBtumor_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBtumor_ugperml*1000)^ndrugactT))
    drugBpbkill = VmB*max(0, TaBratio_pb)^nkill/(max(0, KmTB_kill)^nkill+max(0, TaBratio_pb)^nkill)*((TDBc_ugperml*1000)/(KdrugB+(TDBc_ugperml*1000)))
    drugBtisskill = VmB*max(0, TaBratio_tiss)^nkill/(max(0, KmTB_kill*fKmTB_kill)^nkill+max(0, TaBratio_tiss)^nkill)*((TDBt_ugperml*1000)/(KdrugB+(TDBt_ugperml*1000)))
    drugBtisskill2 = VmB*max(0, TaBratio_tiss2)^nkill/(max(0, KmTB_kill*fKmTB_kill)^nkill+max(0, TaBratio_tiss2)^nkill)*((TDBt2_ugperml*1000)/(KdrugB+(TDBt2_ugperml*1000)))
    drugBtisskill3 = VmB*max(0,TaB1920ratio_tiss3)^nkill/(max(0,KmTB_kill*fKmTB_kill)^nkill+max(0,TaB1920ratio_tiss3)^nkill)*((TDBt3_ugperml*1000)/(KdrugB+(TDBt3_ugperml*1000)))
    drugBtumorkill = VmB*max(0, TaBratio_tumor)^nkill/(max(0, KmTB_kill*fKmTB_kill)^nkill+max(0, TaBratio_tumor)^nkill)*((TDBtumor_ugperml*1000)/(KdrugB+(TDBtumor_ugperml*1000)))

    # T cells 
    du.actTtiss = ((kTact*(drugTtissact*restTtiss-fTadeact*actTtiss)) + (tissue1on*kTaexit*Vpb*((1+fa0*finj*injection_effect)*actTpb/Vpb*KTrp*fTap-actTtiss/Vtissue)) + (act0on*kTact*(drugTtissact*act0Ttiss-fTadeact*actTtiss)) )
    du.actTpb = (
        (kTact*(drugTpbact*restTpb-fTadeact*actTpb)) - 
        (tissue1on*kTaexit*Vpb*((1+fa0*finj*injection_effect)*actTpb/Vpb*KTrp*fTap-actTtiss/Vtissue)) - 
        (kTaapop*(actTpb+fAICD*actTpb^2/(Vpb*Trpbref_perml))) + (act0on*kTact*(drugTpbact*act0Tpb-fTadeact*actTpb)) - 
        (tissue2on*kTaexit*Vpb*((1+fa0*finj*injection_effect)*actTpb/Vpb*KTrp2*fTap-actTtiss2/Vtissue2)) - 
        (tissue3on*kTaexit*Vpb*((1+fa0*finj*injection_effect)*actTpb/Vpb*KTrp3*fTap-actTtiss3/Vtissue3)) - 
        (tumor_on*kTaexit*Vpb*((1+fa0*finj*injection_effect)*actTpb/Vpb*KTrptumor*fTap-actTtumor/Vtumor))
        )
    du.restTtiss = (-(kTact*(drugTtissact*restTtiss-fTadeact*actTtiss)) + (tissue1on*kTrexit*Vpb*((1+finj*injection_effect)*restTpb/Vpb*KTrp-restTtiss/Vtissue)) + (act0on*fTa0deact*act0Ttiss) )
    du.restTpb = (
        -(tissue1on*kTrexit*Vpb*((1+finj*injection_effect)*restTpb/Vpb*KTrp-restTtiss/Vtissue)) + 
        (kTgen*Vpb*Trpbref_perml*(max(0,1-totTpb_perml/Trpbref_perml)^2)) - 
        (kTact*(drugTpbact*restTpb-fTadeact*actTpb)) - fTrapop*kTaapop*max(0,restTpb-Trpbref_perml*Vpb) + (act0on*fTa0deact*act0Tpb) - 
        (tissue2on*kTrexit*Vpb*((1+finj*injection_effect)*restTpb/Vpb*KTrp2-restTtiss2/Vtissue2)) - 
        (tissue3on*kTrexit*Vpb*((1+finj*injection_effect)*restTpb/Vpb*KTrp3-restTtiss3/Vtissue3)) - 
        (tumor_on*kTrexit*Vpb*((1+finj*injection_effect)*restTpb/Vpb*KTrptumor-restTtumor/Vtumor))
        )
    du.act0Ttiss2 = ((act0on*fTaprolif*kTprolif*actTtiss2) + (tissue2on*act0on*kTrexit*Vpb*((1+finj*injection_effect)*act0Tpb/Vpb*KTrp2-act0Ttiss2/Vtissue2)) - (act0on*fTa0deact*act0Ttiss2) - (act0on*kTact*(drugTtissact2*act0Ttiss2-fTadeact*actTtiss2)) - (fTa0apop*kTaapop*act0Ttiss2))
    du.restTtiss2 = ((act0on*fTa0deact*act0Ttiss2) + (tissue2on*kTrexit*Vpb*((1+finj*injection_effect)*restTpb/Vpb*KTrp2-restTtiss2/Vtissue2)) - (kTact*(drugTtissact2*restTtiss2-fTadeact*actTtiss2)) )
    du.actTtiss2 = ( (act0on*kTact*(drugTtissact2*act0Ttiss2-fTadeact*actTtiss2)) + (tissue2on*kTaexit*Vpb*((1+fa0*finj*injection_effect)*actTpb/Vpb*KTrp2*fTap-actTtiss2/Vtissue2)) + (kTact*(drugTtissact2*restTtiss2-fTadeact*actTtiss2)) )
    du.act0Tpb = (
        -(act0on*fTa0deact*act0Tpb) - (act0on*kTact*(drugTpbact*act0Tpb-fTadeact*actTpb)) - (fTa0apop*kTaapop*act0Tpb) - 
        (tissue1on*act0on*kTrexit*Vpb*((1+finj*injection_effect)*act0Tpb/Vpb*KTrp-act0Ttiss/Vtissue)) + (act0on*fTaprolif*kTprolif*actTpb) - 
        (tissue2on*act0on*kTrexit*Vpb*((1+finj*injection_effect)*act0Tpb/Vpb*KTrp2-act0Ttiss2/Vtissue2)) - 
        (tissue3on*act0on*kTrexit*Vpb*((1+finj*injection_effect)*act0Tpb/Vpb*KTrp3-act0Ttiss3/Vtissue3)) - 
        (tumor_on*act0on*kTrexit*Vpb*((1+finj*injection_effect)*act0Tpb/Vpb*KTrptumor-act0Ttumor/Vtumor))
        )
    du.act0Ttiss = (-(fTa0apop*kTaapop*act0Ttiss) - (act0on*kTact*(drugTtissact*act0Ttiss-fTadeact*actTtiss)) - (act0on*fTa0deact*act0Ttiss) + (tissue1on*act0on*kTrexit*Vpb*((1+finj*injection_effect)*act0Tpb/Vpb*KTrp-act0Ttiss/Vtissue)) + (act0on*fTaprolif*kTprolif*actTtiss))
    du.restTtiss3 = ((tissue3on*kTrexit*Vpb*((1+finj*injection_effect)*restTpb/Vpb*KTrp3-restTtiss3/Vtissue3)) + (act0on*fTa0deact*act0Ttiss3) - (kTact*(drugTtissact3*restTtiss3-fTadeact*actTtiss3)) )
    du.act0Ttiss3 = ((tissue3on*act0on*kTrexit*Vpb*((1+finj*injection_effect)*act0Tpb/Vpb*KTrp3-act0Ttiss3/Vtissue3)) - (fTa0apop*kTaapop*act0Ttiss3) - (act0on*kTact*(drugTtissact3*act0Ttiss3-fTadeact*actTtiss3)) - (act0on*fTa0deact*act0Ttiss3) + (act0on*fTaprolif*kTprolif*actTtiss3))
    du.actTtiss3 = ((tissue3on*kTaexit*Vpb*((1+fa0*finj*injection_effect)*actTpb/Vpb*KTrp3*fTap-actTtiss3/Vtissue3)) + (act0on*kTact*(drugTtissact3*act0Ttiss3-fTadeact*actTtiss3)) + (kTact*(drugTtissact3*restTtiss3-fTadeact*actTtiss3)) )
    du.restTtumor = (
        (act0on*fTa0deact*act0Ttumor) - (kTact*(drugTtumoract*restTtumor-fTadeact*actTtumor)) + 
        (tumor_on*kTrexit*Vpb*((1+finj*injection_effect)*restTpb/Vpb*KTrptumor-restTtumor/Vtumor)) 
        )
    du.actTtumor = (
        (act0on*kTact*(drugTtumoract*act0Ttumor-fTadeact*actTtumor)) + 
        (kTact*(drugTtumoract*restTtumor-fTadeact*actTtumor))  + 
        (tumor_on*kTaexit*Vpb*((1+fa0*finj*injection_effect)*actTpb/Vpb*KTrptumor*fTap-actTtumor/Vtumor))
        )
    du.act0Ttumor = (
        -(act0on*fTa0deact*act0Ttumor) - (act0on*kTact*(drugTtumoract*act0Ttumor-fTadeact*actTtumor)) + (act0on*fTaprolif*kTprolif*actTtumor) + 
        (tumor_on*act0on*kTrexit*Vpb*((1+finj*injection_effect)*act0Tpb/Vpb*KTrptumor-act0Ttumor/Vtumor)) - (fTa0apop*kTaapop*act0Ttumor)
        )

    # B cell
    du.Bpb = (
        -tissue1on*fBexit*kTaexit*Vpb*max(0,Bpb/Vpb*KBp-Btiss/Vtissue) + tissue3on*fBtissue3_v1*kBapop*Bpb*kBapop_cll - kBapop*Bpb*kBapop_cll - 
        kBkill*drugBpbkill*Bpb - 
        loncaBbp*Bpb - 
        tissue2on*fBexit*kTaexit*Vpb*max(0,Bpb/Vpb*KBp2-Btiss2/Vtissue2) + 
        tissue3on*fBtissue3_v1*kBtiss3exit*Vpb*max(0, B1920tiss3/Vtissue3-Bpb/Vpb*KBp3) + 
        tissue3on*fBtissue3_v1*kBprolif*kBapop_cll*Bpbref_perml*Vpb*max(0,1-Bpb/(Bpbref_perml*Vpb)) - 
        Bcell_tumor_trafficking_on*tumor_on*fBexit*kTaexit*Vpb*max(0,Bpb/Vpb*KBptumor-(Btumor + Btumor_trans + Btumor_cd19neg)/Vtumor)
        )
    du.Btiss = (
        (kBprolif*KBp*Bpbref_perml*Vtissue*(max(0,1-Btiss/(KBp*Bpbref_perml*Vtissue)))^1) + 
        (tissue1on*fBexit*kTaexit*Vpb*max(0,Bpb/Vpb*KBp-Btiss/Vtissue)) - 
        (fBkill*kBkill*drugBtisskill*Btiss) - 
        loncaBtissueact1*Btiss
        )
    du.Btiss2 = (
        -(fBkill*kBkill*drugBtisskill2*Btiss2 + loncaBtissueact2*Btiss2) + 
        (tissue2on*fBexit*kTaexit*Vpb*max(0,Bpb/Vpb*KBp2-Btiss2/Vtissue2)) + 
        (kBprolif*KBp2*Bpbref_perml*Vtissue2*(max(0,1-Btiss2/(KBp2*Bpbref_perml*Vtissue2)))^1)
        )
    du.B1920tiss3 = (
        kBapop/kBmat_kBapop_ratio*B19no20tiss3 - tissue3on*fBtissue3_v1*kBtiss3exit*Vpb*max(0,B1920tiss3/Vtissue3-Bpb/Vpb*KBp3) - kBapop*B1920tiss3*kBapop_cll - 
        loncaBtissueact3*B1920tiss3 - 
        fBkill*kBkill*drugBtisskill3*B1920tiss3 + 
        kBprolif*KBp3*Bpbref_perml*Vtissue3*max(0,1-B1920tiss3/(KBp3*Bpbref_perml*Vtissue3))
        )
    du.B19no20tiss3 = (
        Bpbref_perml*KBp3*B19no20_B1920_ratio*Vtissue3*(kBapop+kBapop/kBmat_kBapop_ratio)*(1+5*(max(0,1-(Bpb+Btiss+Btiss2)/(Bpbref_perml*Vpb+KBp*Bpbo_perml*Vtissue+KBp2*Bpbo_perml*Vtissue2)))^1) - 
        kBapop/kBmat_kBapop_ratio*B19no20tiss3 - kBapop*B19no20tiss3*kBapop_cll - fBkill*kBkill*B19no20tiss3 - loncaBtissueact3*B19no20tiss3 + 
        kBprolif*KBp3*Bpbref_perml*B19no20_B1920_ratio*Vtissue3*max(0,1-B19no20tiss3/(KBp3*Bpbref_perml*B19no20_B1920_ratio*Vtissue3))
        )
    du.Btumor = (
        (kBtumorprolif*Btumor) - 
        (fBkill*kBkill*drugBtumorkill*Btumor + loncaBtumoract*Btumor) + 
        (Bcell_tumor_trafficking_on*tumor_on*fBexit*kTaexit*Vpb*max(0,Bpb/Vpb*KBptumor-(Btumor + Btumor_trans + Btumor_cd19neg)/Vtumor))
    )
    du.Btumor_trans = loncaBtumoract*Btumor - k_trans*Btumor_trans
    du.Btumor_cd19neg = ( kBtumorprolif*Btumor_cd19neg - fBkill*kBkill*drugBtumorkill*Btumor_cd19neg - loncaBtumoractneg*Btumor_cd19neg )
    du.Btumor_cd19neg_trans = loncaBtumoractneg*Btumor_cd19neg - k_trans_neg*Btumor_cd19neg_trans
    
    # T-cell dependent bispecific (TDB)
    TDBdepot_ug = max(TDBdepot_ug, 0.);
    du.TDBc_ugperml = infusion_TDB + kabs_TDB * TDBdepot_ug/V1_TDB - 1/V1_TDB * Q_TDB * (TDBc_ugperml - TDBp_ugperml) - CL_TDB/V1_TDB * TDBc_ugperml - Vm_tdb/V1_TDB / (Km_tdb + TDBc_ugperml) * TDBc_ugperml
    du.TDBp_ugperml = 1/V2_TDB * Q_TDB * (TDBc_ugperml - TDBp_ugperml)
    du.TDBdepot_ug = - kabs_TDB * TDBdepot_ug
    du.TDBdepot_ug = TDBdepot_ug >= 0.0 ? du.TDBdepot_ug : 0.0
    du.TDBc_ugperml = TDBc_ugperml >= 0.0 ? du.TDBc_ugperml : 0.0


    # lonca 
    du.LONCAc_ugperml = infusion_lonca - 1/V1_lonca * Q_lonca * (LONCAc_ugperml - LONCAp_ugperml) - CL_lonca/V1_lonca * LONCAc_ugperml
    du.LONCAp_ugperml = 1/V2_lonca * Q_lonca * (LONCAc_ugperml - LONCAp_ugperml);
    du.LONCAc_ugperml = LONCAc_ugperml >= 0.0 ? du.LONCAc_ugperml : 0.
    du.LONCAp_ugperml = LONCAp_ugperml >= 0.0 ? du.LONCAp_ugperml : 0.

    # IL6 
    du.IL6pb = (
        kIL6prod*actTpb/Vpb*(Bpb_perml/Bpbref_perml)*(TDBc_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBc_ugperml*1000)^ndrugactT) - 
        log(2)/(thalfIL6/HourToMinute/DayToHour)*IL6pb
        )
    du.IL6tiss = (
        kIL6prod*actTtiss/Vtissue*(Btiss_perml/(Bpbref_perml*KBp))*(TDBt_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBt_ugperml*1000)^ndrugactT) - 
        log(2)/(thalfIL6/HourToMinute/DayToHour)*IL6tiss
        )
    du.IL6tiss2 = (
        kIL6prod*actTtiss2/Vtissue2*(Btiss2_perml/(Bpbref_perml*KBp2))*((TDBt2_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBt2_ugperml*1000)^ndrugactT)) - 
        log(2)/(thalfIL6/HourToMinute/DayToHour)*IL6tiss2
        )
    du.IL6tiss3 = (
        kIL6prod*actTtiss3/Vtissue3*(  B1920tiss3_perml/(Bpbref_perml*KBp3)*(TDBt3_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBt3_ugperml*1000)^ndrugactT)  ) - 
        log(2)/(thalfIL6/HourToMinute/DayToHour)*IL6tiss3
        )
    du.IL6tumor = (
        kIL6prod*actTtumor/Vtumor*(Btumor_perml/(Bpbref_perml*KBptumor))*((TDBtumor_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBtumor_ugperml*1000)^ndrugactT)) - 
        log(2)/(thalfIL6/HourToMinute/DayToHour)*IL6tumor
        )
    
    #  
    du.injection_effect = (-(log(2)/tinjhalf*injection_effect))

    # track cell loss due to Lonca or TDB 
    du.TDBKill = (fBkill*kBkill*drugBtumorkill*Btumor + fBkill*kBkill*drugBtumorkill*Btumor_cd19neg)
    du.LoncaKill = (loncaBtumoract*Btumor + loncaBtumoractneg*Btumor_cd19neg)

    # println(du)
end