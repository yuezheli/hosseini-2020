

using ComponentArrays
using Parameters: @unpack

function blin_ode!(du, u, p, t)
    # unpack variables
    @unpack Blinc_ug, Btumor, 
    actTtiss, Btiss, actTpb, restTtiss, restTpb,
    Bpb, act0Tpb, act0Ttiss, injection_effect, Btiss2, act0Ttiss2, restTtiss2, 
    actTtiss2, B1920tiss3, B19no20tiss3, 
    restTtiss3, act0Ttiss3, actTtiss3, restTtumor, actTtumor, act0Ttumor, 
    IL6pb, IL6tiss, IL6tiss2, IL6tiss3, IL6tumor = u

    # unpack parameters
    @unpack Kp, Kp2, Kp3, Kptumor, 
    # B cell
    VmB , kBapop, kBprolif, kBkill, KBp2, Bpbo_perml, fBexit, KBp, Bpbref_perml, kBtiss3exit, KBp3, B19no20_B1920_ratio, fBtissue3_v1, 
    kBmat_kBapop_ratio, kBapop_cll, fBkill, fKmTB_kill, KBptumor, kBtumorprolif, 
    # T cell
    kTprolif, VmT, kTaexit, kTact, fTadeact, fTap, kTaapop, KTrp, fTaprolif, kTgen, fTrapop, Trpbref_perml, KTrp2, fa0, fAICD, fTa0deact, fTa0apop, finj, 
    kTrexit, KTrp3, fTact, KTrptumor, S, 
    # IL6 
    kIL6prod, thalfIL6, IL6_tiss_contribution, 
    # parameters undocumented 
    act0on, tissue1on, tissue2on, tissue3on, tumor_on, tinjhalf, Vtumor, 
    # parameters not used in the ODE 
    Trpbo_perml, 
    # params updated for human 
    Vpb, Vtissue, Vtissue2, Vtissue3, 
    # params for blincyto PK
    Cl_blin, Vz_blin, KdrugactT_blin, KmTB_kill_blin, KdrugB_blin, nkill_blin, KmBT_act_blin, blin_infusion = p
    
    
    Blinc_ngperml = max(Blinc_ug/Vz_blin*1000, 1E-18)
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

    BlinTpbact = VmT*(BTrratio_pb^S/(KmBT_act_blin^S+BTrratio_pb^S))*(Blinc_ngperml/(Blinc_ngperml+KdrugactT_blin))
    BlinTtissact = fTact*VmT*(BTrratio_tiss^S/(KmBT_act_blin^S+BTrratio_tiss^S))*(Blint_ngperml/(Blint_ngperml+KdrugactT_blin))
    BlinTtissact2 = fTact*VmT*(BTrratio_tiss2^S/(KmBT_act_blin^S+BTrratio_tiss2^S))*(Blint2_ngperml/(Blint2_ngperml+KdrugactT_blin))
    BlinTtissact3 = fTact*VmT*(B19Trratio_tiss3^S/(KmBT_act_blin^S+B19Trratio_tiss3^S))*(Blint3_ngperml/(Blint3_ngperml+KdrugactT_blin))
    BlinTtumoract = fTact*VmT*(max(BTrratio_tumor,0)^S/(KmBT_act_blin^S+max(BTrratio_tumor, 0)^S))*(Blintumor_ngperml/(Blintumor_ngperml+KdrugactT_blin))
    BlinBpbkill = VmB*max(0, TaBratio_pb)^nkill_blin/(max(0, KmTB_kill_blin)^nkill_blin+max(0, TaBratio_pb)^nkill_blin)*(Blinc_ngperml/(KdrugB_blin+Blinc_ngperml))
    BlinBtisskill = VmB*max(0, TaBratio_tiss)^nkill_blin/(max(0, KmTB_kill_blin*fKmTB_kill)^nkill_blin+max(0, TaBratio_tiss)^nkill_blin)*(Blint_ngperml/(KdrugB_blin+Blint_ngperml))
    BlinBtisskill2 = VmB*max(0, TaBratio_tiss2)^nkill_blin/(max(0, KmTB_kill_blin*fKmTB_kill)^nkill_blin+max(0, TaBratio_tiss2)^nkill_blin)*(Blint2_ngperml/(KdrugB_blin+Blint2_ngperml))
    BlinBtisskill3 = VmB*max(0,TaB19ratio_tiss3)^nkill_blin/(max(0,KmTB_kill_blin*fKmTB_kill)^nkill_blin+max(0,TaB19ratio_tiss3)^nkill_blin)*(Blint3_ngperml/(KdrugB_blin+Blint3_ngperml))
    BlinBtumorkill = VmB*max(0, TaBratio_tumor)^nkill_blin/(max(0, KmTB_kill_blin*fKmTB_kill)^nkill_blin+max(0, TaBratio_tumor)^nkill_blin)*(Blintumor_ngperml/(KdrugB_blin+Blintumor_ngperml))

    # T cells
    du.actTtiss = ((kTact*(BlinTtissact*restTtiss-fTadeact*actTtiss)) + (tissue1on*kTaexit*Vpb*((1+fa0*finj*injection_effect)*actTpb/Vpb*KTrp*fTap-actTtiss/Vtissue)) + (act0on*kTact*(BlinTtissact*act0Ttiss-fTadeact*actTtiss)) )
    du.actTpb = (
        (kTact*(BlinTpbact*restTpb-fTadeact*actTpb)) - 
        (tissue1on*kTaexit*Vpb*((1+fa0*finj*injection_effect)*actTpb/Vpb*KTrp*fTap-actTtiss/Vtissue)) - 
        (kTaapop*(actTpb+fAICD*actTpb*actTpb/(Vpb*Trpbref_perml))) + (act0on*kTact*(BlinTpbact*act0Tpb-fTadeact*actTpb)) - 
        (tissue2on*kTaexit*Vpb*((1+fa0*finj*injection_effect)*actTpb/Vpb*KTrp2*fTap-actTtiss2/Vtissue2)) - 
        (tissue3on*kTaexit*Vpb*((1+fa0*finj*injection_effect)*actTpb/Vpb*KTrp3*fTap-actTtiss3/Vtissue3)) - 
        (tumor_on*kTaexit*Vpb*((1+fa0*finj*injection_effect)*actTpb/Vpb*KTrptumor*fTap-actTtumor/Vtumor))
        )
    du.restTtiss = (-(kTact*(BlinTtissact*restTtiss-fTadeact*actTtiss)) + (tissue1on*kTrexit*Vpb*((1+finj*injection_effect)*restTpb/Vpb*KTrp-restTtiss/Vtissue)) + (act0on*fTa0deact*act0Ttiss) )
    du.restTpb = (
        -(tissue1on*kTrexit*Vpb*((1+finj*injection_effect)*restTpb/Vpb*KTrp-restTtiss/Vtissue)) + 
        (kTgen*Vpb*Trpbref_perml*(max(0,1-totTpb_perml/Trpbref_perml)^2)) - 
        (kTact*(BlinTpbact*restTpb-fTadeact*actTpb)) - fTrapop*kTaapop*max(0,restTpb-Trpbref_perml*Vpb) + (act0on*fTa0deact*act0Tpb) - 
        (tissue2on*kTrexit*Vpb*((1+finj*injection_effect)*restTpb/Vpb*KTrp2-restTtiss2/Vtissue2)) - 
        (tissue3on*kTrexit*Vpb*((1+finj*injection_effect)*restTpb/Vpb*KTrp3-restTtiss3/Vtissue3)) - 
        (tumor_on*kTrexit*Vpb*((1+finj*injection_effect)*restTpb/Vpb*KTrptumor-restTtumor/Vtumor))
        )
    du.act0Ttiss2 = ((act0on*fTaprolif*kTprolif*actTtiss2) + (tissue2on*act0on*kTrexit*Vpb*((1+finj*injection_effect)*act0Tpb/Vpb*KTrp2-act0Ttiss2/Vtissue2)) - (act0on*fTa0deact*act0Ttiss2) - (act0on*kTact*(BlinTtissact2*act0Ttiss2-fTadeact*actTtiss2)) - (fTa0apop*kTaapop*act0Ttiss2))
    du.restTtiss2 = ((act0on*fTa0deact*act0Ttiss2) + (tissue2on*kTrexit*Vpb*((1+finj*injection_effect)*restTpb/Vpb*KTrp2-restTtiss2/Vtissue2)) - (kTact*(BlinTtissact2*restTtiss2-fTadeact*actTtiss2)) )
    du.actTtiss2 = ( (act0on*kTact*(BlinTtissact2*act0Ttiss2-fTadeact*actTtiss2)) + (tissue2on*kTaexit*Vpb*((1+fa0*finj*injection_effect)*actTpb/Vpb*KTrp2*fTap-actTtiss2/Vtissue2)) + (kTact*(BlinTtissact2*restTtiss2-fTadeact*actTtiss2)) )
    du.act0Tpb = (
        -(act0on*fTa0deact*act0Tpb) - (act0on*kTact*(BlinTpbact*act0Tpb-fTadeact*actTpb)) - (fTa0apop*kTaapop*act0Tpb) - 
        (tissue1on*act0on*kTrexit*Vpb*((1+finj*injection_effect)*act0Tpb/Vpb*KTrp-act0Ttiss/Vtissue)) + (act0on*fTaprolif*kTprolif*actTpb) - 
        (tissue2on*act0on*kTrexit*Vpb*((1+finj*injection_effect)*act0Tpb/Vpb*KTrp2-act0Ttiss2/Vtissue2)) - 
        (tissue3on*act0on*kTrexit*Vpb*((1+finj*injection_effect)*act0Tpb/Vpb*KTrp3-act0Ttiss3/Vtissue3)) - 
        (tumor_on*act0on*kTrexit*Vpb*((1+finj*injection_effect)*act0Tpb/Vpb*KTrptumor-act0Ttumor/Vtumor))
        )
    du.act0Ttiss = (-(fTa0apop*kTaapop*act0Ttiss) - (act0on*kTact*(BlinTtissact*act0Ttiss-fTadeact*actTtiss)) - (act0on*fTa0deact*act0Ttiss) + (tissue1on*act0on*kTrexit*Vpb*((1+finj*injection_effect)*act0Tpb/Vpb*KTrp-act0Ttiss/Vtissue)) + (act0on*fTaprolif*kTprolif*actTtiss))
    du.restTtiss3 = ((tissue3on*kTrexit*Vpb*((1+finj*injection_effect)*restTpb/Vpb*KTrp3-restTtiss3/Vtissue3)) + (act0on*fTa0deact*act0Ttiss3) - (kTact*(BlinTtissact3*restTtiss3-fTadeact*actTtiss3)) )
    du.act0Ttiss3 = ((tissue3on*act0on*kTrexit*Vpb*((1+finj*injection_effect)*act0Tpb/Vpb*KTrp3-act0Ttiss3/Vtissue3)) - (fTa0apop*kTaapop*act0Ttiss3) - (act0on*kTact*(BlinTtissact3*act0Ttiss3-fTadeact*actTtiss3)) - (act0on*fTa0deact*act0Ttiss3) + (act0on*fTaprolif*kTprolif*actTtiss3))
    du.actTtiss3 = ((tissue3on*kTaexit*Vpb*((1+fa0*finj*injection_effect)*actTpb/Vpb*KTrp3*fTap-actTtiss3/Vtissue3)) + (act0on*kTact*(BlinTtissact3*act0Ttiss3-fTadeact*actTtiss3)) + (kTact*(BlinTtissact3*restTtiss3-fTadeact*actTtiss3)) )
    du.restTtumor = (
        (act0on*fTa0deact*act0Ttumor) - (kTact*(BlinTtumoract*restTtumor-fTadeact*actTtumor)) + 
        (tumor_on*kTrexit*Vpb*((1+finj*injection_effect)*restTpb/Vpb*KTrptumor-restTtumor/Vtumor)) 
        )
    du.actTtumor = (
        (act0on*kTact*(BlinTtumoract*act0Ttumor-fTadeact*actTtumor)) + 
        (kTact*(BlinTtumoract*restTtumor-fTadeact*actTtumor))  + 
        (tumor_on*kTaexit*Vpb*((1+fa0*finj*injection_effect)*actTpb/Vpb*KTrptumor*fTap-actTtumor/Vtumor))
        )
    du.act0Ttumor = (
        -(act0on*fTa0deact*act0Ttumor) - (act0on*kTact*(BlinTtumoract*act0Ttumor-fTadeact*actTtumor)) + (act0on*fTaprolif*kTprolif*actTtumor) + 
        (tumor_on*act0on*kTrexit*Vpb*((1+finj*injection_effect)*act0Tpb/Vpb*KTrptumor-act0Ttumor/Vtumor)) - (fTa0apop*kTaapop*act0Ttumor)
        )
    # B cell
    du.Bpb = (
        -tissue1on*fBexit*kTaexit*Vpb*max(0,Bpb/Vpb*KBp-Btiss/Vtissue) + tissue3on*fBtissue3_v1*kBapop*Bpb*kBapop_cll - kBapop*Bpb*kBapop_cll - 
        kBkill*BlinBpbkill*Bpb - 
        tissue2on*fBexit*kTaexit*Vpb*max(0,Bpb/Vpb*KBp2-Btiss2/Vtissue2) + 
        tissue3on*fBtissue3_v1*kBtiss3exit*Vpb*max(0, B1920tiss3/Vtissue3-Bpb/Vpb*KBp3) + 
        tissue3on*fBtissue3_v1*kBprolif*kBapop_cll*Bpbref_perml*Vpb*max(0,1-Bpb/(Bpbref_perml*Vpb)) 
        )
    du.Btiss = (
        (kBprolif*KBp*Bpbref_perml*Vtissue*(max(0,1-Btiss/(KBp*Bpbref_perml*Vtissue)))) + 
        (tissue1on*fBexit*kTaexit*Vpb*max(0,Bpb/Vpb*KBp-Btiss/Vtissue)) - 
        (fBkill*kBkill*BlinBtisskill*Btiss) 
        )
    du.Btiss2 = (
        -fBkill*kBkill*BlinBtisskill2*Btiss2 + 
        (tissue2on*fBexit*kTaexit*Vpb*max(0,Bpb/Vpb*KBp2-Btiss2/Vtissue2)) + 
        (kBprolif*KBp2*Bpbref_perml*Vtissue2*(max(0,1-Btiss2/(KBp2*Bpbref_perml*Vtissue2))))
        )
    du.B1920tiss3 = (
        kBapop/kBmat_kBapop_ratio*B19no20tiss3 - tissue3on*fBtissue3_v1*kBtiss3exit*Vpb*max(0,B1920tiss3/Vtissue3-Bpb/Vpb*KBp3) - kBapop*B1920tiss3*kBapop_cll - 
        fBkill*kBkill*BlinBtisskill3*B1920tiss3 + 
        kBprolif*KBp3*Bpbref_perml*Vtissue3*max(0,1-B1920tiss3/(KBp3*Bpbref_perml*Vtissue3))
        )
    du.B19no20tiss3 = (
        Bpbref_perml*KBp3*B19no20_B1920_ratio*Vtissue3*(kBapop+kBapop/kBmat_kBapop_ratio)*(1+5*(max(0,1-(Bpb+Btiss+Btiss2)/(Bpbref_perml*Vpb+KBp*Bpbo_perml*Vtissue+KBp2*Bpbo_perml*Vtissue2)))^1) - 
        fBkill*kBkill*BlinBtisskill3*B19no20tiss3 - 
        kBapop/kBmat_kBapop_ratio*B19no20tiss3 - kBapop*B19no20tiss3*kBapop_cll - fBkill*kBkill*B19no20tiss3 + 
        kBprolif*KBp3*Bpbref_perml*B19no20_B1920_ratio*Vtissue3*max(0,1-B19no20tiss3/(KBp3*Bpbref_perml*B19no20_B1920_ratio*Vtissue3))
        )
    du.Btumor = (
        (kBtumorprolif*Btumor) - 
        (fBkill*kBkill*BlinBtumorkill*Btumor) 
    )
    
    # blinatumomab
    du.Blinc_ug = (-Cl_blin/Vz_blin*Blinc_ug + blin_infusion)

    # IL6 
    du.IL6pb = kIL6prod*actTpb/Vpb*(Bpb_perml/Bpbref_perml)*Blinc_ngperml/(Blinc_ngperml+KdrugactT_blin) - log(2)/thalfIL6*IL6pb
    du.IL6tiss = kIL6prod*actTtiss/Vtissue*(Btiss_perml/(Bpbref_perml*KBp))*Blint_ngperml/(Blint_ngperml+KdrugactT_blin) - log(2)/thalfIL6*IL6tiss
    du.IL6tiss2 = kIL6prod*actTtiss2/Vtissue2*(Btiss2_perml/(Bpbref_perml*KBp2))*Blint2_ngperml/(Blint2_ngperml+KdrugactT_blin) - log(2)/thalfIL6*IL6tiss2
    du.IL6tiss3 = kIL6prod*actTtiss3/Vtissue3*(B1920tiss3_perml/(Bpbref_perml*KBp3))* B19tiss3_perml/(Bpbref_perml*KBp3*(1+B19no20_B1920_ratio))*Blint3_ngperml/(Blint3_ngperml+KdrugactT_blin) - log(2)/thalfIL6*IL6tiss3
    du.IL6tumor = kIL6prod*actTtumor/Vtumor*(Btumor_perml/(Bpbref_perml*KBptumor))*Blintumor_ngperml/(Blintumor_ngperml+KdrugactT_blin) - log(2)/thalfIL6*IL6tumor
    
    # 
    du.injection_effect = (-(log(2)/tinjhalf*injection_effect))

end