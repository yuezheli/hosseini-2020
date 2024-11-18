# author: Yuezhe LI
# date: 10/3/2023
# purpose of this code: to initialize the model 
# this model was based on the cleaned up model from bs-lonca-2.jl

include("param.jl");

# pre-equilibriate init condition, except tumor 
u0_0 = ComponentArray(
    restTpb = p_homo_3.Trpbo_perml * p_homo_3.Vpb,
    restTtiss = p_homo_3.KTrp * p_homo_3.Trpbo_perml * p_homo_3.Vtissue, 
    restTtiss2 = p_homo_3.KTrp2 * p_homo_3.Trpbo_perml * p_homo_3.Vtissue2, 
    restTtiss3 = p_homo_3.KTrp3 * p_homo_3.Trpbo_perml * p_homo_3.Vtissue3,
    restTtumor = 5.7E10,  
    Bpb = p_homo_3.Bpbo_perml * p_homo_3.Vpb, 
    Btiss = p_homo_3.KBp * p_homo_3.Bpbo_perml * p_homo_3.Vtissue, 
    Btiss2 = p_homo_3.KBp2 * p_homo_3.Bpbo_perml * p_homo_3.Vtissue2, 
    B1920tiss3 = p_homo_3.Bpbref_perml * p_homo_3.KBp3 * p_homo_3.Vtissue3, 
    B19no20tiss3 = p_homo_3.Bpbref_perml * p_homo_3.KBp3 * p_homo_3.B19no20_B1920_ratio * p_homo_3.Vtissue3, 
    Btumor_cd19neg = 0., 
    Btumor_cd19neg_trans = 0., 
    Btumor = (p_homo_3.Vtumor*0.375)/Vc, 
    Btumor_trans = 0., 
    actTpb=0., actTtiss=0., actTtiss2=0., actTtiss3=0., actTtumor=0., act0Tpb=0., act0Ttiss=0., act0Ttiss2=0., act0Ttiss3=0., act0Ttumor=0., 
    IL6pb=0., IL6tiss=0., IL6tiss2=0., IL6tiss3=0., IL6tumor=0., injection_effect = 0., 
    TDBc_ugperml= 0., TDBp_ugperml=0. , TDBdepot_ug = 0., LONCAc_ugperml = 0., LONCAp_ugperml = 0., 
    TDBKill = 0., LoncaKill = 0. );
