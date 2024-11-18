using ComponentArrays
using Parameters: @unpack

include("Hosseini_params.jl");

function Hosseini_init(p)
    @unpack VmB , KmTB_kill, KdrugB, kBapop, kBprolif, kTprolif, kBkill, KdrugactT, VmT, KmBT_act, kTaexit, kTact, fTadeact, fTap, Cl_tdb, Cld_tdb, 
    Vc_tdb, Vp_tdb, Ktgen, kTaapop, KTrp, Trpbo_perml, fTaprolif, fBprolif, Bpbo_perml, fBexit, kTgen, KBp, Kp, Vpb, Vtissue, fTrapop, Trpbref_perml, 
    nkill, KTrp2, Vm_tdb, Km_tdb, BW, act0on, kIL6prod, fa0, fAICD, kBAFFprod, Bpbref_perml, thBAFF, thalfIL6, fTa0deact, fTa0apop, BAFFo, fBAFFo, 
    tinjhalf, finj, KBp2, Vtissue2, Kp2, tissue2on, Vc_rtx, Vp_rtx, Cl_rtx, Cld_rtx, Vm_rtx, Km_rtx, Kmkill_rtx, fTgenbl, depleteTpb, depleteBpb, 
    kTrexit, PKflag, VPid, end_time, fdrug, Vtissue3, KTrp3, tissue3on, Kp3, kBtiss3exit, KBp3, B19no20_B1920_ratio, fvalidation, fapop_v24 , fBtissue3_v1, 
    kBmat_kBapop_ratio, ndrugactT, Cl_blin, Vz_blin, KdrugactT_blin, ndrugactT_blin, KmTB_kill_blin, KdrugB_blin, nkill_blin, KmBT_act_blin, kBapop_cll, 
    kBgen_cll, fBkill, fTact, fKmTB_kill, Kptumor, Vtumor, KBptumor, KTrptumor, tumor_on, IL6_tiss_contribution, kBtumorprolif, S, tissue1on, kabs_TDB, 
    fbio_TDB, Bcell_tumor_trafficking_on = p

    actTtiss = 0.0
    Btiss = KBp*Bpbo_perml*Vtissue
    TDBc_ugperkg = 1000.0   # [μg/kg]
    actTpb = 0.0
    restTtiss = KTrp*Trpbo_perml*Vtissue
    TDBp_ugperkg = 0.0   # [μg/kg]
    restTpb = (1-depleteTpb)*Trpbo_perml*Vpb
    Bpb = (1-depleteBpb)*Bpbo_perml*Vpb
    BAFF = BAFFo
    act0Tpb = 0.0
    act0Ttiss = 0.0
    injection_effect = 0.0
    Btiss2 = KBp2*Bpbo_perml*Vtissue2
    act0Ttiss2 = 0.0
    restTtiss2 = KTrp2*Trpbo_perml*Vtissue2
    actTtiss2 = 0.0
    RTXc_ugperkg = 0.0
    RTXp_ugperkg = 0.0
    drug_effect = 0.0
    B1920tiss3 = Bpbref_perml*KBp3*Vtissue3
    B19no20tiss3 = Bpbref_perml*KBp3*B19no20_B1920_ratio*Vtissue3
    restTtiss3 = KTrp3*Trpbo_perml*Vtissue3
    act0Ttiss3 = 0.0
    actTtiss3 = 0.0
    Blinc_ug = 0.0
    restTtumor = KTrptumor*Trpbo_perml*Vtumor
    actTtumor = 0.0
    Btumor = KBptumor*Bpbo_perml*Vtumor
    act0Ttumor = 0.0
    TDBsc_ugperkg = 0.0
    TDBc_ugperml_AUC = 0.0
    IL6pb = 0.0
    IL6tiss = 0.0
    IL6tiss2 = 0.0
    IL6tiss3 = 0.0
    IL6tumor = 0.0

    u0 = ComponentArray(actTtiss=actTtiss, Btiss=Btiss, TDBc_ugperkg=TDBc_ugperkg, actTpb=actTpb, restTtiss=restTtiss, TDBp_ugperkg=TDBp_ugperkg, restTpb=restTpb,
    Bpb=Bpb, BAFF=BAFF, act0Tpb=act0Tpb, act0Ttiss=act0Ttiss, injection_effect=injection_effect, Btiss2=Btiss2, act0Ttiss2=act0Ttiss2, restTtiss2=restTtiss2, 
    actTtiss2=actTtiss2, RTXc_ugperkg=RTXc_ugperkg, RTXp_ugperkg=RTXp_ugperkg, drug_effect=drug_effect, B1920tiss3=B1920tiss3, B19no20tiss3=B19no20tiss3, 
    restTtiss3=restTtiss3, act0Ttiss3=act0Ttiss3, actTtiss3=actTtiss3, Blinc_ug=Blinc_ug, restTtumor=restTtumor, actTtumor=actTtumor, Btumor=Btumor, act0Ttumor=act0Ttumor, 
    TDBsc_ugperkg=TDBsc_ugperkg, TDBc_ugperml_AUC=TDBc_ugperml_AUC, IL6pb=IL6pb, IL6tiss=IL6tiss, IL6tiss2=IL6tiss2, IL6tiss3=IL6tiss3, IL6tumor=IL6tumor);

    return u0;
end
