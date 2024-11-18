using ComponentArrays
using Parameters: @unpack

include("PK.jl")

function Hosseini_ode!(du, u, p, t)
    # unpack variables
    @unpack actTtiss, Btiss, TDBc_ugperkg, actTpb, restTtiss, TDBp_ugperkg, restTpb,
    Bpb, BAFF, act0Tpb, act0Ttiss, injection_effect, Btiss2, act0Ttiss2, restTtiss2, 
    actTtiss2, RTXc_ugperkg, RTXp_ugperkg, drug_effect, B1920tiss3, B19no20tiss3, 
    restTtiss3, act0Ttiss3, actTtiss3, Blinc_ug, restTtumor, actTtumor, Btumor, act0Ttumor, 
    TDBsc_ugperkg, TDBc_ugperml_AUC, IL6pb, IL6tiss, IL6tiss2, IL6tiss3, IL6tumor = u

    # unpack parameters
    @unpack VmB , KmTB_kill, KdrugB, kBapop, kBprolif, kTprolif, kBkill, KdrugactT, VmT, KmBT_act, kTaexit, kTact, fTadeact, fTap, Cl_tdb, Cld_tdb, 
    Vc_tdb, Vp_tdb, Ktgen, kTaapop, KTrp, Trpbo_perml, fTaprolif, fBprolif, Bpbo_perml, fBexit, kTgen, KBp, Kp, Vpb, Vtissue, fTrapop, Trpbref_perml, 
    nkill, KTrp2, Vm_tdb, Km_tdb, BW, act0on, kIL6prod, fa0, fAICD, kBAFFprod, Bpbref_perml, thBAFF, thalfIL6, fTa0deact, fTa0apop, BAFFo, fBAFFo, 
    tinjhalf, finj, KBp2, Vtissue2, Kp2, tissue2on, Vc_rtx, Vp_rtx, Cl_rtx, Cld_rtx, Vm_rtx, Km_rtx, Kmkill_rtx, fTgenbl, depleteTpb, depleteBpb, 
    kTrexit, PKflag, VPid, end_time, fdrug, Vtissue3, KTrp3, tissue3on, Kp3, kBtiss3exit, KBp3, B19no20_B1920_ratio, fvalidation, fapop_v24 , fBtissue3_v1, 
    kBmat_kBapop_ratio, ndrugactT, Cl_blin, Vz_blin, KdrugactT_blin, ndrugactT_blin, KmTB_kill_blin, KdrugB_blin, nkill_blin, KmBT_act_blin, kBapop_cll, 
    kBgen_cll, fBkill, fTact, fKmTB_kill, Kptumor, Vtumor, KBptumor, KTrptumor, tumor_on, IL6_tiss_contribution, kBtumorprolif, S, tissue1on, kabs_TDB, 
    fbio_TDB, Bcell_tumor_trafficking_on = p

    ## Repeated Assignments
    TDBc_ugperml = PK_v26(TDBc_ugperkg, Vc_tdb, PKflag, VPid, time, end_time, fvalidation)
    println(TDBc_ugperml)
    TDBt_ugperml = Kp*TDBc_ugperml
    TDBt2_ugperml = Kp2*TDBc_ugperml
    TDBt3_ugperml = tissue3on*Kp3*TDBc_ugperml
    TDBtumor_ugperml = tumor_on*Kptumor*TDBc_ugperml
    RTXc_ugperml = (RTXc_ugperkg/Vc_rtx>1e-5)*RTXc_ugperkg/Vc_rtx
    RTXt_ugperml = RTXc_ugperml*Kp
    RTXt2_ugperml = RTXc_ugperml*Kp2
    RTXt3_ugperml = tissue3on*RTXc_ugperml*Kp3
    RTXtumor_ugperml = tumor_on*RTXc_ugperml*Kptumor
    Blinc_ngperml = Blinc_ug/Vz_blin*1000
    Blint_ngperml = Kp*Blinc_ngperml
    Blint2_ngperml = Kp2*Blinc_ngperml
    Blint3_ngperml = tissue3on*Kp3*Blinc_ngperml
    Blintumor_ngperml = tumor_on*Kptumor*Blinc_ngperml
    restTpb_perml = restTpb/Vpb
    restTtiss_perml = restTtiss/Vtissue
    restTtiss2_perml = restTtiss2/Vtissue2
    restTtiss3_perml = restTtiss3/Vtissue3
    restTtumor_perml = restTtumor/Vtumor
    act0Tpb_perml = act0Tpb/Vpb
    act0Ttiss_perml = act0Ttiss/Vtissue
    act0Ttiss2_perml = act0Ttiss2/Vtissue2
    act0Ttiss3_perml = act0Ttiss3/Vtissue3
    act0Ttumor_perml = act0Ttumor/Vtumor
    actTpb_perml = actTpb/Vpb
    actTtiss_perml = actTtiss/Vtissue
    actTtiss2_perml = actTtiss2/Vtissue2
    actTtiss3_perml = actTtiss3/Vtissue3
    actTtumor_perml = actTtumor/Vtumor
    B19tiss3 = B19no20tiss3+B1920tiss3
    Bpb_perml = Bpb/Vpb
    Btiss_perml = Btiss/Vtissue
    Btiss2_perml = Btiss2/Vtissue2
    B19tiss3_perml = B19tiss3/Vtissue3
    B1920tiss3_perml = B1920tiss3/Vtissue3
    B19no20tiss3_perml = B19no20tiss3/Vtissue3
    Btumor_perml = Btumor/Vtumor
    BTrratio_pb = Bpb/max(restTpb+act0Tpb,1)
    BTrratio_tiss = Btiss/max(restTtiss+act0Ttiss,1)
    BTrratio_tiss2 = Btiss2/max(restTtiss2+act0Ttiss2,1)
    B19Trratio_tiss3 = B19tiss3/max(restTtiss3+act0Ttiss3,1)
    B1920Trratio_tiss3 = B1920tiss3/max(restTtiss3+act0Ttiss3,1)
    BTrratio_tumor = Btumor/max(restTtumor+act0Ttumor,1)
    TaBratio_pb = actTpb/max(Bpb,1)
    TaBratio_tiss = actTtiss/max(Btiss,1)
    TaBratio_tiss2 = actTtiss2/max(Btiss2,1)
    TaB19ratio_tiss3 = actTtiss3/max(B19tiss3,1)
    TaB1920ratio_tiss3 = actTtiss3/max(B1920tiss3,1)
    TaBratio_tumor = actTtumor/max(Btumor,1)
    totTpb_perml = restTpb_perml+act0Tpb_perml + actTpb_perml
    totTtiss_perml = restTtiss_perml + act0Ttiss_perml + actTtiss_perml
    totTtiss2_perml = restTtiss2_perml + act0Ttiss2_perml + actTtiss2_perml
    totTtiss3_perml = restTtiss3_perml + act0Ttiss3_perml + actTtiss3_perml
    totTtumor_perml = restTtumor_perml + act0Ttumor_perml + actTtumor_perml
    Tafraction_pb = actTpb_perml/max(totTpb_perml,1)
    Tafraction_tiss = actTtiss_perml/max(totTtiss_perml,1)
    Tafraction_tiss2 = actTtiss2_perml/max(totTtiss2_perml,1)
    Tafraction_tiss3 = actTtiss3_perml/max(totTtiss3_perml,1)
    Tafraction_tumor = actTtumor_perml/max(totTtumor_perml,1)

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
    BlinTpbact = VmT*(BTrratio_pb^S/(KmBT_act_blin^S+BTrratio_pb^S))*(Blinc_ngperml^ndrugactT_blin/(Blinc_ngperml^ndrugactT_blin+KdrugactT_blin^ndrugactT_blin))
    BlinTtissact = fTact*VmT*(BTrratio_tiss^S/(KmBT_act_blin^S+BTrratio_tiss^S))*(Blint_ngperml^ndrugactT_blin/(Blint_ngperml^ndrugactT_blin+KdrugactT_blin^ndrugactT_blin))
    BlinTtissact2 = fTact*VmT*(BTrratio_tiss2^S/(KmBT_act_blin^S+BTrratio_tiss2^S))*(Blint2_ngperml^ndrugactT_blin/(Blint2_ngperml^ndrugactT_blin+KdrugactT_blin^ndrugactT_blin))
    BlinTtissact3 = fTact*VmT*(B19Trratio_tiss3^S/(KmBT_act_blin^S+B19Trratio_tiss3^S))*(Blint3_ngperml^ndrugactT_blin/(Blint3_ngperml^ndrugactT_blin+KdrugactT_blin^ndrugactT_blin))
    BlinTtumoract = fTact*VmT*(max(BTrratio_tumor,0)^S/(KmBT_act_blin^S+max(BTrratio_tumor, 0)^S))*(Blintumor_ngperml^ndrugactT_blin/(Blintumor_ngperml^ndrugactT_blin+KdrugactT_blin^ndrugactT_blin))
    BlinBpbkill = VmB*max(0, TaBratio_pb)^nkill_blin/(max(0, KmTB_kill_blin)^nkill_blin+max(0, TaBratio_pb)^nkill_blin)*(Blinc_ngperml/(KdrugB_blin+Blinc_ngperml))
    BlinBtisskill = VmB*max(0, TaBratio_tiss)^nkill_blin/(max(0, KmTB_kill_blin*fKmTB_kill)^nkill_blin+max(0, TaBratio_tiss)^nkill_blin)*(Blint_ngperml/(KdrugB_blin+Blint_ngperml))
    BlinBtisskill2 = VmB*max(0, TaBratio_tiss2)^nkill_blin/(max(0, KmTB_kill_blin*fKmTB_kill)^nkill_blin+max(0, TaBratio_tiss2)^nkill_blin)*(Blint2_ngperml/(KdrugB_blin+Blint2_ngperml))
    BlinBtisskill3 = VmB*max(0,TaB19ratio_tiss3)^nkill_blin/(max(0,KmTB_kill_blin*fKmTB_kill)^nkill_blin+max(0,TaB19ratio_tiss3)^nkill_blin)*(Blint3_ngperml/(KdrugB_blin+Blint3_ngperml))
    BlinBtumorkill = VmB*max(0, TaBratio_tumor)^nkill_blin/(max(0, KmTB_kill_blin*fKmTB_kill)^nkill_blin+max(0, TaBratio_tumor)^nkill_blin)*(Blintumor_ngperml/(KdrugB_blin+Blintumor_ngperml))
    Bpb_norm = Bpb_perml/Bpbref_perml
    totTtiss = restTtiss + act0Ttiss + actTtiss
    totTtiss2 = restTtiss2 + act0Ttiss2 + actTtiss2
    totTtiss3 = restTtiss3 + act0Ttiss3 + actTtiss3
    totTtumor = restTtumor + act0Ttumor + actTtumor
    Baffconsumption = log(2)/thBAFF*BAFF*(Btiss+Bpb + Btiss2)/(Bpbref_perml*(Vpb+KBp*Vtissue+KBp2*Vtissue2))
    BTtotRatio_tiss = Btiss_perml/totTtiss_perml
    BTtotRatio_tiss2 = Btiss2_perml/totTtiss2_perml
    IL6combo = IL6pb+IL6_tiss_contribution*(IL6tiss*Vtissue+IL6tiss2*Vtissue2+IL6tiss3*Vtissue3+IL6tumor*Vtumor)/Vpb

    # T cells
    du.actTtiss = ((kTact*((drugTtissact+BlinTtissact)*restTtiss-fTadeact*actTtiss)) + (tissue1on*kTaexit*Vpb*((1+fa0*(finj*injection_effect+fdrug*drug_effect))*actTpb/Vpb*KTrp*fTap-actTtiss/Vtissue)) + (act0on*kTact*((drugTtissact+BlinTtissact)*act0Ttiss-fTadeact*actTtiss)) - (fapop_v24*kTaapop*(actTtiss+fAICD*actTtiss^2/(KTrp*Vtissue*Trpbref_perml))))
    du.actTpb = ((kTact*((drugTpbact+BlinTpbact)*restTpb-fTadeact*actTpb)) - (tissue1on*kTaexit*Vpb*((1+fa0*(finj*injection_effect+fdrug*drug_effect))*actTpb/Vpb*KTrp*fTap-actTtiss/Vtissue)) - (kTaapop*(actTpb+fAICD*actTpb^2/(Vpb*Trpbref_perml))) + (act0on*kTact*((drugTpbact+BlinTpbact)*act0Tpb-fTadeact*actTpb)) - (tissue2on*kTaexit*Vpb*((1+fa0*(finj*injection_effect+fdrug*drug_effect))*actTpb/Vpb*KTrp2*fTap-actTtiss2/Vtissue2)) - (tissue3on*kTaexit*Vpb*((1+fa0*(finj*injection_effect+fdrug*drug_effect))*actTpb/Vpb*KTrp3*fTap-actTtiss3/Vtissue3)) - (tumor_on*kTaexit*Vpb*((1+fa0*(finj*injection_effect+fdrug*drug_effect))*actTpb/Vpb*KTrptumor*fTap-actTtumor/Vtumor)))
    du.restTtiss = (-(kTact*((drugTtissact+BlinTtissact)*restTtiss-fTadeact*actTtiss)) + (0) + (tissue1on*kTrexit*Vpb*((1+finj*injection_effect+fdrug*drug_effect)*restTpb/Vpb*KTrp-restTtiss/Vtissue)) + (act0on*fTa0deact*act0Ttiss) - (fapop_v24*fTrapop*kTaapop*(Trpbref_perml*Vtissue*KTrp)*max(0,restTtiss/(Trpbref_perml*Vtissue*KTrp)-1)))
    du.restTpb = (-(tissue1on*kTrexit*Vpb*((1+finj*injection_effect+fdrug*drug_effect)*restTpb/Vpb*KTrp-restTtiss/Vtissue)) + (kTgen*Vpb*Trpbref_perml*(fTgenbl+max(0,1-totTpb_perml/Trpbref_perml)^2)) - (kTact*((drugTpbact+BlinTpbact)*restTpb-fTadeact*actTpb)) - (kTgen*fTgenbl*restTpb+fTrapop*kTaapop*max(0,restTpb-Trpbref_perml*Vpb)) + (act0on*fTa0deact*act0Tpb) - (tissue2on*kTrexit*Vpb*((1+finj*injection_effect+fdrug*drug_effect)*restTpb/Vpb*KTrp2-restTtiss2/Vtissue2)) - (tissue3on*kTrexit*Vpb*((1+finj*injection_effect+fdrug*drug_effect)*restTpb/Vpb*KTrp3-restTtiss3/Vtissue3)) - (tumor_on*kTrexit*Vpb*((1+finj*injection_effect+fdrug*drug_effect)*restTpb/Vpb*KTrptumor-restTtumor/Vtumor)))
    du.act0Ttiss2 = ((act0on*fTaprolif*kTprolif*actTtiss2) + (tissue2on*act0on*kTrexit*Vpb*((1+finj*injection_effect+fdrug*drug_effect)*act0Tpb/Vpb*KTrp2-act0Ttiss2/Vtissue2)) - (act0on*fTa0deact*act0Ttiss2) - (act0on*kTact*((drugTtissact2+BlinTtissact2)*act0Ttiss2-fTadeact*actTtiss2)) - (fTa0apop*kTaapop*act0Ttiss2))
    du.restTtiss2 = ((act0on*fTa0deact*act0Ttiss2) + (tissue2on*kTrexit*Vpb*((1+finj*injection_effect+fdrug*drug_effect)*restTpb/Vpb*KTrp2-restTtiss2/Vtissue2)) - (kTact*((drugTtissact2+BlinTtissact2)*restTtiss2-fTadeact*actTtiss2)) - (fapop_v24*fTrapop*kTaapop*(Trpbref_perml*Vtissue2*KTrp2)*max(0,restTtiss2/(Trpbref_perml*Vtissue2*KTrp2)-1)))
    du.actTtiss2 = ((act0on*kTact*((drugTtissact2+BlinTtissact2)*act0Ttiss2-fTadeact*actTtiss2)) + (tissue2on*kTaexit*Vpb*((1+fa0*(finj*injection_effect+fdrug*drug_effect))*actTpb/Vpb*KTrp2*fTap-actTtiss2/Vtissue2)) + (kTact*((drugTtissact2+BlinTtissact2)*restTtiss2-fTadeact*actTtiss2)) - (fapop_v24*kTaapop*(actTtiss2+fAICD*actTtiss2^2/(Vtissue2*KTrp2*Trpbref_perml))))
    du.act0Tpb = (-(act0on*fTa0deact*act0Tpb) - (act0on*kTact*((drugTpbact+BlinTpbact)*act0Tpb-fTadeact*actTpb)) - (fTa0apop*kTaapop*act0Tpb) - (tissue1on*act0on*kTrexit*Vpb*((1+finj*injection_effect+fdrug*drug_effect)*act0Tpb/Vpb*KTrp-act0Ttiss/Vtissue)) + (act0on*fTaprolif*kTprolif*actTpb) - (tissue2on*act0on*kTrexit*Vpb*((1+finj*injection_effect+fdrug*drug_effect)*act0Tpb/Vpb*KTrp2-act0Ttiss2/Vtissue2)) - (tissue3on*act0on*kTrexit*Vpb*((1+finj*injection_effect+fdrug*drug_effect)*act0Tpb/Vpb*KTrp3-act0Ttiss3/Vtissue3)) - (tumor_on*act0on*kTrexit*Vpb*((1+finj*injection_effect+fdrug*drug_effect)*act0Tpb/Vpb*KTrptumor-act0Ttumor/Vtumor)))
    du.act0Ttiss = (-(fTa0apop*kTaapop*act0Ttiss) - (act0on*kTact*((drugTtissact+BlinTtissact)*act0Ttiss-fTadeact*actTtiss)) - (act0on*fTa0deact*act0Ttiss) + (tissue1on*act0on*kTrexit*Vpb*((1+finj*injection_effect+fdrug*drug_effect)*act0Tpb/Vpb*KTrp-act0Ttiss/Vtissue)) + (act0on*fTaprolif*kTprolif*actTtiss))
    du.restTtiss3 = ((tissue3on*kTrexit*Vpb*((1+finj*injection_effect+fdrug*drug_effect)*restTpb/Vpb*KTrp3-restTtiss3/Vtissue3)) + (act0on*fTa0deact*act0Ttiss3) - (kTact*((drugTtissact3+BlinTtissact3)*restTtiss3-fTadeact*actTtiss3)) - (fapop_v24*fTrapop*kTaapop*(Trpbref_perml*Vtissue3*KTrp3)*max(0,restTtiss3/(Trpbref_perml*Vtissue3*KTrp3)-1)))
    du.act0Ttiss3 = ((tissue3on*act0on*kTrexit*Vpb*((1+finj*injection_effect+fdrug*drug_effect)*act0Tpb/Vpb*KTrp3-act0Ttiss3/Vtissue3)) - (fTa0apop*kTaapop*act0Ttiss3) - (act0on*kTact*((drugTtissact3+BlinTtissact3)*act0Ttiss3-fTadeact*actTtiss3)) - (act0on*fTa0deact*act0Ttiss3) + (act0on*fTaprolif*kTprolif*actTtiss3))
    du.actTtiss3 = ((tissue3on*kTaexit*Vpb*((1+fa0*(finj*injection_effect+fdrug*drug_effect))*actTpb/Vpb*KTrp3*fTap-actTtiss3/Vtissue3)) + (act0on*kTact*((drugTtissact3+BlinTtissact3)*act0Ttiss3-fTadeact*actTtiss3)) + (kTact*((drugTtissact3+BlinTtissact3)*restTtiss3-fTadeact*actTtiss3)) - (fapop_v24*kTaapop*(actTtiss3+fAICD*actTtiss3^2/(Vtissue3*KTrp3*Trpbref_perml))))
    du.restTtumor = ((act0on*fTa0deact*act0Ttumor) - (kTact*((drugTtumoract+BlinTtumoract)*restTtumor-fTadeact*actTtumor)) + (tumor_on*kTrexit*Vpb*((1+finj*injection_effect+fdrug*drug_effect)*restTpb/Vpb*KTrptumor-restTtumor/Vtumor)) - (fapop_v24*fTrapop*kTaapop*(Trpbref_perml*Vtumor*KTrptumor)*max(0,restTtumor/(Trpbref_perml*Vtumor*KTrptumor)-1)))
    du.actTtumor = ((act0on*kTact*((drugTtumoract+BlinTtumoract)*act0Ttumor-fTadeact*actTtumor)) + (kTact*((drugTtumoract+BlinTtumoract)*restTtumor-fTadeact*actTtumor)) - (fapop_v24*kTaapop*(actTtumor+fAICD*actTtumor^2/(Vtumor*KTrptumor*Trpbref_perml))) + (tumor_on*kTaexit*Vpb*((1+fa0*(finj*injection_effect+fdrug*drug_effect))*actTpb/Vpb*KTrptumor*fTap-actTtumor/Vtumor)))
    du.act0Ttumor = (-(act0on*fTa0deact*act0Ttumor) - (act0on*kTact*((drugTtumoract+BlinTtumoract)*act0Ttumor-fTadeact*actTtumor)) + (act0on*fTaprolif*kTprolif*actTtumor) + (tumor_on*act0on*kTrexit*Vpb*((1+finj*injection_effect+fdrug*drug_effect)*act0Tpb/Vpb*KTrptumor-act0Ttumor/Vtumor)) - (fTa0apop*kTaapop*act0Ttumor))

    # B cell
    du.Bpb = (-(tissue1on*fBexit*kTaexit*Vpb*max(0,Bpb/Vpb*KBp-Btiss/Vtissue)) + (tissue3on*fBtissue3_v1*kBapop*Bpb*kBapop_cll) - (kBapop*Bpb*kBapop_cll) - (kBkill*(drugBpbkill+BlinBpbkill+RTXc_ugperml/(Kmkill_rtx/1000+RTXc_ugperml))*Bpb) - (tissue2on*fBexit*kTaexit*Vpb*max(0,Bpb/Vpb*KBp2-Btiss2/Vtissue2)) + (tissue3on*fBtissue3_v1*kBtiss3exit*Vpb*max(0,B1920tiss3/Vtissue3-Bpb/Vpb*KBp3)) + (tissue3on*fBtissue3_v1*kBprolif*kBapop_cll*Bpbref_perml*Vpb*(max(0,1-Bpb/(Bpbref_perml*Vpb)))^1) - (Bcell_tumor_trafficking_on*tumor_on*fBexit*kTaexit*Vpb*max(0,Bpb/Vpb*KBptumor-Btumor/Vtumor)))
    du.Btiss = ((kBprolif*KBp*Bpbref_perml*Vtissue*(max(0,1-Btiss/(KBp*Bpbref_perml*Vtissue)))^1) - (0*fBprolif*kBprolif*Btiss*kBapop_cll) + (tissue1on*fBexit*kTaexit*Vpb*max(0,Bpb/Vpb*KBp-Btiss/Vtissue)) - (fBkill*kBkill*(drugBtisskill+BlinBtisskill+RTXt_ugperml/(Kmkill_rtx/1000+RTXt_ugperml))*Btiss))
    du.Btiss2 = (-(fBkill*kBkill*(drugBtisskill2+BlinBtisskill2+RTXt2_ugperml/(Kmkill_rtx/1000+RTXt2_ugperml))*Btiss2) + (tissue2on*fBexit*kTaexit*Vpb*max(0,Bpb/Vpb*KBp2-Btiss2/Vtissue2)) - (0*fBprolif*kBprolif*Btiss2*kBapop_cll) + (kBprolif*KBp2*Bpbref_perml*Vtissue2*(max(0,1-Btiss2/(KBp2*Bpbref_perml*Vtissue2)))^1))
    du.B1920tiss3 = (((kBapop/kBmat_kBapop_ratio)*B19no20tiss3) - (tissue3on*fBtissue3_v1*kBtiss3exit*Vpb*max(0,B1920tiss3/Vtissue3-Bpb/Vpb*KBp3)) - (kBapop*B1920tiss3*kBapop_cll) - (fBkill*kBkill*(drugBtisskill3+BlinBtisskill3+RTXt3_ugperml/(Kmkill_rtx/1000+RTXt3_ugperml))*B1920tiss3) + (kBprolif*KBp3*Bpbref_perml*Vtissue3*(max(0,1-B1920tiss3/(KBp3*Bpbref_perml*Vtissue3)))^1))
    du.B19no20tiss3 = ((Bpbref_perml*KBp3*B19no20_B1920_ratio*Vtissue3*(kBapop+(kBapop/kBmat_kBapop_ratio))*(1+5*(max(0,1-(Bpb+Btiss+Btiss2)/(Bpbref_perml*Vpb+KBp*Bpbo_perml*Vtissue+KBp2*Bpbo_perml*Vtissue2)))^1)) - ((kBapop/kBmat_kBapop_ratio)*B19no20tiss3) - (kBapop*B19no20tiss3*kBapop_cll) - (fBkill*kBkill*BlinBtisskill3*B19no20tiss3) + (kBprolif*KBp3*Bpbref_perml*B19no20_B1920_ratio*Vtissue3*(max(0,1-B19no20tiss3/(KBp3*Bpbref_perml*B19no20_B1920_ratio*Vtissue3)))^1))
    du.Btumor = ((kBtumorprolif*Btumor+0*kBprolif*KBptumor*Bpbref_perml*Vtumor*(max(0,1-Btumor/(KBptumor*Bpbref_perml*Vtumor)))^1) - (0*fBprolif*kBprolif*Btumor) - (fBkill*kBkill*(drugBtumorkill+BlinBtumorkill+RTXtumor_ugperml/(Kmkill_rtx/1000+RTXtumor_ugperml))*Btumor) + (Bcell_tumor_trafficking_on*tumor_on*fBexit*kTaexit*Vpb*max(0,Bpb/Vpb*KBptumor-Btumor/Vtumor)))

    # T-cell dependent bispecific (TDB)
    du.TDBc_ugperkg = (-((Cl_tdb+Vm_tdb/(Km_tdb+TDBc_ugperkg/Vc_tdb)/BW)/Vc_tdb*TDBc_ugperkg) - ((Cld_tdb*(TDBc_ugperkg/Vc_tdb-TDBp_ugperkg/Vp_tdb))) + (kabs_TDB*fbio_TDB*TDBsc_ugperkg))
    du.TDBp_ugperkg = (((Cld_tdb*(TDBc_ugperkg/Vc_tdb-TDBp_ugperkg/Vp_tdb))))
    du.TDBsc_ugperkg = (-(kabs_TDB*fbio_TDB*TDBsc_ugperkg) - (kabs_TDB*(1-fbio_TDB)*TDBsc_ugperkg))
    du.TDBc_ugperml_AUC = ((TDBc_ugperml))

    # RTX
    du.RTXc_ugperkg = (-(Cld_rtx*(RTXc_ugperkg/Vc_rtx-RTXp_ugperkg/Vp_rtx)) - ((Cl_rtx+Vm_rtx/(Km_rtx+RTXc_ugperkg/Vc_rtx)/BW)/Vc_rtx*RTXc_ugperkg))
    du.RTXp_ugperkg = ((Cld_rtx*(RTXc_ugperkg/Vc_rtx-RTXp_ugperkg/Vp_rtx)))
    du.drug_effect = (-(log(2)/tinjhalf*drug_effect))

    # blinatumomab
    du.Blinc_ug = (-(Cl_blin*Blinc_ug/Vz_blin))

    # IL6 
    du.IL6pb = ((kIL6prod*actTpb/Vpb*(Bpb_perml/Bpbref_perml)*((TDBc_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBc_ugperml*1000)^ndrugactT)+Blinc_ngperml^ndrugactT_blin/(Blinc_ngperml^ndrugactT_blin+KdrugactT_blin^ndrugactT_blin))) - (log(2)/(thalfIL6/60/24)*IL6pb))
    du.IL6tiss = ((kIL6prod*actTtiss/Vtissue*(Btiss_perml/(Bpbref_perml*KBp))*((TDBt_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBt_ugperml*1000)^ndrugactT)+Blint_ngperml^ndrugactT_blin/(Blint_ngperml^ndrugactT_blin+KdrugactT_blin^ndrugactT_blin))) - (log(2)/(thalfIL6/60/24)*IL6tiss))
    du.IL6tiss2 = ((kIL6prod*actTtiss2/Vtissue2*(Btiss2_perml/(Bpbref_perml*KBp2))*((TDBt2_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBt2_ugperml*1000)^ndrugactT)+Blint2_ngperml^ndrugactT_blin/(Blint2_ngperml^ndrugactT_blin+KdrugactT_blin^ndrugactT_blin))) - (log(2)/(thalfIL6/60/24)*IL6tiss2))
    du.IL6tiss3 = ((kIL6prod*actTtiss3/Vtissue3*((B1920tiss3_perml/(Bpbref_perml*KBp3))*((TDBt3_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBt3_ugperml*1000)^ndrugactT))+(B19tiss3_perml/(Bpbref_perml*KBp3*(1+B19no20_B1920_ratio)))*(Blint3_ngperml^ndrugactT_blin/(Blint3_ngperml^ndrugactT_blin+KdrugactT_blin^ndrugactT_blin)))) - (log(2)/(thalfIL6/60/24)*IL6tiss3))
    du.IL6tumor = ((kIL6prod*actTtumor/Vtumor*(Btumor_perml/(Bpbref_perml*KBptumor))*((TDBtumor_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBtumor_ugperml*1000)^ndrugactT)+Blintumor_ngperml^ndrugactT_blin/(Blintumor_ngperml^ndrugactT_blin+KdrugactT_blin^ndrugactT_blin))) - (log(2)/(thalfIL6/60/24)*IL6tumor))
    
    # 
    du.BAFF = (-(log(2)/(thBAFF/24/60)*(fBAFFo*BAFF/BAFFo+(1-fBAFFo)*(Btiss+Bpb+Btiss2)/(Bpbref_perml*(Vpb+KBp*Vtissue+KBp2*Vtissue2)))) + (log(2)/(thBAFF/24/60)*(fBAFFo+(1-fBAFFo)*Bpbo_perml/Bpbref_perml)))
    du.injection_effect = (-(log(2)/tinjhalf*injection_effect))

end