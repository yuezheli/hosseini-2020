# date: 11/1/2024
# author: Yuezhe Li 
# purpose of this code: to create default initial condition for blinatumomab simulation 

include("param_blin.jl");

# pre-equilibriate init condition
u0_0 = ComponentArray(
    restTpb = p_homo_blin.Trpbo_perml * p_homo_blin.Vpb,
    restTtiss = p_homo_blin.KTrp * p_homo_blin.Trpbo_perml * p_homo_blin.Vtissue, 
    restTtiss2 = p_homo_blin.KTrp2 * p_homo_blin.Trpbo_perml * p_homo_blin.Vtissue2, 
    restTtiss3 = p_homo_blin.KTrp3 * p_homo_blin.Trpbo_perml * p_homo_blin.Vtissue3,
    restTtumor = p_homo_blin.Kptumor * p_homo_blin.Trpbo_perml * p_homo_blin.Vtumor,
    Bpb = p_homo_blin.Bpbo_perml * p_homo_blin.Vpb, 
    Btiss = p_homo_blin.KBp * p_homo_blin.Bpbo_perml * p_homo_blin.Vtissue, 
    Btiss2 = p_homo_blin.KBp2 * p_homo_blin.Bpbo_perml * p_homo_blin.Vtissue2, 
    B1920tiss3 = p_homo_blin.Bpbref_perml * p_homo_blin.KBp3 * p_homo_blin.Vtissue3, 
    B19no20tiss3 = p_homo_blin.Bpbref_perml * p_homo_blin.KBp3 * p_homo_blin.B19no20_B1920_ratio * p_homo_blin.Vtissue3, 
    Btumor = (p_homo_blin.Vtumor*0.375)/Vc, 
    actTpb=0., actTtiss=0., actTtiss2=0., actTtiss3=0., actTtumor=0., act0Tpb=0., act0Ttiss=0., act0Ttiss2=0., act0Ttiss3=0., act0Ttumor=0., 
    IL6pb=0., IL6tiss=0., IL6tiss2=0., IL6tiss3=0., IL6tumor=0., injection_effect = 0., 
    Blinc_ug = 0.);
