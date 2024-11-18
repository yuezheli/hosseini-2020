# author: Yuezhe Li
# date: 6/5/23
# purpose: this table was obtained from Hosseini  et al., 2020 supplementary MATLAB code 
# https://www.nature.com/articles/s41540-020-00145-7

using ComponentArrays
using Parameters: @unpack

# MATLAB params
# note not all params matches either cyno or human parameters, it's merely a list of params extracted from MATLAB code 
# this is kept as it is to reproduce MATLAB result; DO NOT MODIFY IN THIS FILE
Hosseini_params = ComponentArray(
    # TDB : T-cell act
    VmT = 0.9,              # maximum rate of T-cell activation [unitless]
    fTact = 1.0,
    KmBT_act = 0.1,
    ndrugactT = 1.0,
    S = 1.0, 
    KdrugactT = 50,         # the TDB conc at which T-cell activation was half of VmT [ng/mL]
    # TDB : B-cell kill
    VmB = 0.95,             # maximum rate of B-cell killing 
    KmTB_kill = 0.75,       # the Tact/B ratio that the rate if half of VmB
    fKmTB_kill = 1.0,
    nkill = 2.0,
    KdrugB = 10,            # the TDB conc at which B-cell killing rate if half of VmB [ng/mL]
    # TDB PK 
    Cl_tdb = 3.0,           # elimination clearance of TDB [mL/d/kg]
    Cld_tdb = 25.0,         # distributon clearance of TDB [mL/d/kg]
    Vc_tdb = 35.0,          # central volume of distribution of TDB [mL/kg]
    Vp_tdb = 110.0,         # peripheral tisue volume of distribution of TDB [mL/kg]
    Vm_tdb = 212 ,          # [μg/d]
    Km_tdb = 4.74 ,         # [μg/d]
    # 
    # Blin : T-cell act
    KmBT_act_blin = 0.716,
    ndrugactT_blin = 1.0,
    KdrugactT_blin = 0.1 ,  # [ng/mL]
    # Blin : B-cell kill
    KmTB_kill_blin = 0.75,
    KdrugB_blin = 0.015  ,  # [ng/mL]
    nkill_blin = 1.027,
    # Blin PK
    Cl_blin = 22300.0,      # systematic clearance of blinatumomab (anti-CD3/CD19) [mL/m2/day]
    Vz_blin = 1610,         # distribution volume of blinatumomab [mL/m^2]
    #
    # B-cells
    kBapop = 0.03,          # [1/d]
    kBkill = 1.0,           # killing rate of B-cells [1/d]
    fBkill = 1.0, 
    kBprolif = 0.7,         # [1/d]
    # 
    # T-cells
    kTprolif = 0.7 ,        
    kTaexit = 1.0,          # trafficking rate of activated CD8+ T-cells from PB to tissue [1/d]
    kTact = 1.0,            # conversion rate of resting or post-activated CD8+ T cells (and vice versa) [1/d]
    fTadeact = 0.0,         # fraction of activated CD8+ T-cells that deactivate to resting or "exhausted" CD8+ T-cells
    fTap = 1.0,              
    kTaapop = 0.5,          # apoptosis rate of CD69+ CD8+ T-cells (activated CD8+ T-cells) [1/d]
    fTaprolif = 2.0,        # proliferation rate of activated CD8+ T-cells [1/d]
    fTrapop = 0.1,          # 
    kTrexit = 1.0,          # trafficking rate of resting/ post-activated CD8+ T-cells from PB to tissue [1/d]
    fdrug = 1.0,
    fa0 = 0.5,
    fAICD = 1.0,
    fTa0deact = 1.0,
    fTa0apop = 1.0,
    finj = 0.1,
    #
    # PB Physiology
    Trpbref_perml = 2.0e6,  # baseline conc of resting CD8+ T-cells in PB [cells/mL]
    Bpbref_perml = 1.0e6,   # baseline conc of CD19+CD20+ B-cells in PB [cells/mL]
    Vpb = 40.0,             # physiological volume of peripheral blood [mL]
    # 
    # Spleen Physiology 
    # Trtiss_perml
    # Btiss_perml
    Vtissue = 40.0,         # physiological volume of spleen [mL]
    KTrp = 1.0,
    KBp = 1.0,
    Kp = 0.14,
    # 
    # LNs Physiology
    # Trtiss2_perml 
    # Btiss2_perml
    Vtissue2 = 1. ,         # total physiological volume of LNs [mL]
    KTrp2 = 0.02, 
    KBp2 = 1.0,
    Kp2 = 0.07,
    # 
    # BM Physiology
    # Trtiss3_perml
    # B1920tiss3_perml
    # B19no20tiss3_perml
    # B19tiss3_perml
    Vtissue3 = 1.0,         # Physiology volume of BM [mL]
    KTrp3 = 500.0,
    KBp3 = 600.0,
    Kp3 = 0.07,
    kBmat_kBapop_ratio = 0.25, 
    B19no20_B1920_ratio = 0.25, 
    # 
    # Tumor Physiology
    tumor_on = 0.0,         # this turns on/ off the tumor compartment
    # Trtumor_perml
    # Btumor_perml
    Vtumor = 1.0, 
    KTrptumor = 1.0,
    KBptumor = 1.0,
    Kptumor = 1.0,
    kBtumorprolif = 0.025,  # proliferation rate of B-cells in tumor [1/d]
    # 
    # Cytokine 
    kIL6prod = 1.0,         # production rate of IL6 [pg/d]
    IL6_tiss_contribution = 0.02,
    thalfIL6 = 20.0,        # [minute]
    #
    # Init-related values
    Trpbo_perml = 2.0e6,    # [1/mL]
    Bpbo_perml = 1.0e6,     # [1/mL]
    BAFFo = 1.0,
    # 
    # RTX 
    Vc_rtx = 28.0,
    Vp_rtx = 9.0,
    Cl_rtx = 8.0,
    Cld_rtx = 23.0,
    Vm_rtx = 0.0,
    Km_rtx = 1.0,
    Kmkill_rtx = 25.0,
    # 
    # not documented in Supp 2 
    act0on = 1.0,
    tissue1on = 1.0,
    tissue2on = 0.0,
    tissue3on = 0.0,
    Bcell_tumor_trafficking_on = 1.0,
    depleteTpb = 0.0,
    depleteBpb = 0.0,
    BW = 3.4 ,  # [kg]
    Ktgen = 1.0,   # [1/d]
    fBprolif = 0.0,
    fBexit = 1.0,
    kTgen = 1.0 ,  # [1/d]
    kBAFFprod = 1.0,
    thBAFF = 30.0,   # [minute]
    fBAFFo = 0.1,
    tinjhalf = 1.0,   # [d]
    fTgenbl = 0.0,
    kBtiss3exit = 0.03,
    fapop_v24 = 0.0,
    fBtissue3_v1 = 1.0,
    kBapop_cll = 1.0,
    kBgen_cll = 1.0,
    kabs_TDB = 1.4,
    fbio_TDB = 0.6,
    PKflag = 1.0,
    VPid = 1.0,
    end_time = 1.0,
    fvalidation = 0.0,
);

# cyno params, if mentioned in Hosseini 2020 Supp Table 2
p_cyno = ComponentArray(
    # TDB : T-cell act
    VmT = 0.9,                      # maximum rate of T-cell activation [unitless]
    fTact = 0.25,                   # VmT is multiplied by this factor (<1) to capture the less efficient activation of tissue T-cells
    KmBT_act = 0.72,                # The B/T rario, at which the rate is half of VmT (at high concentrations of TDB)
    ndrugactT = 0.8,                # Hill coefficient for the TDB portion of the Michaelis-Menten equation
    S = 1.4,                        # Hill coefficient for the B:T portion of the Michaelis-Menten equation
    KdrugactT = 130.08,             # the TDB conc at which T-cell activation was half of VmT [ng/mL]
    # TDB : B-cell kill
    VmB = 0.95,                     # maximum rate of B-cell killing 
    KmTB_kill = 0.75,               # the Tact/B ratio that the rate if half of VmB
    fKmTB_kill = 10.0,              # KmTB_kill and KmTB_kill_blin are multiplied by this factor (>1) to compensate for the fact the tissues are not well-stirred environments like PB
    nkill = 1.03,                   # Hill coefficient for the Tact/B portion of the Michaelis-Menten equation
    KdrugB = 1.3,                   # the TDB conc at which B-cell killing rate if half of VmB [ng/mL]
    # TDB PK 
    Cl_tdb = 8.5,                   # elimination clearance of TDB [mL/d/kg]
    Cld_tdb = 24.17,                # distributon clearance of TDB [mL/d/kg]
    Vc_tdb = 36.78,                 # central volume of distribution of TDB [mL/kg]
    Vp_tdb = 173.47,                # peripheral tisue volume of distribution of TDB [mL/kg]
    # 
    # Blin : T-cell act
    KmBT_act_blin = 0.716,          # The B/T rario, at which the rate is half of VmT (at high concentrations of Blin)
    ndrugactT_blin = 1.0,           # Hill coefficient for the Blin portion of the Michaelis-Menten equation
    KdrugactT_blin = 0.1 ,          # The Blin concentration, at which T-cell activation rate is half of VmT (at high B/T ratios) [ng/mL]
    # Blin : B-cell kill
    KmTB_kill_blin = 0.75,          # The Tact/B rario, at which the rate is half of VmB (at high concentrations of Blin)
    nkill_blin = 1.027,             # Hill coefficient for the Tact/B portion of the Michaelis-Menten equation
    KdrugB_blin = 0.15  ,          # The Blin concentration, at which B-cell killing rate is half of VmB (at high Tact/B ratios) [ng/mL]
    # Blin PK
    Cl_blin = 22300.0,              # systematic clearance of blinatumomab (anti-CD3/CD19); fixed based on literature [mL/m2/day]
    Vz_blin = 1610.,                # distribution volume of blinatumomab [mL/m^2]; fixed based on literature 
    #
    # B-cells
    kBapop = 0.02,                  # Apoptosis rate of B-cells [1/d]
    kBkill = 275.19,                # killing rate of B-cells [1/d]
    fBkill = 0.1,                   # kBkill is multiplied by this factor (<1) to capture the less efficient killing of tissue B-cells
    kBprolif = 0.05,                # proliferation rate of B cells [1/d]
    # 
    # T-cells
    kTprolif = 0.7 ,                # fraction of activated cells that proliferate
    kTaexit = 0.12,                 # trafficking rate of activated CD8+ T-cells from PB to tissue [1/d]
    kTact = 9.83,                   # conversion rate of resting or post-activated CD8+ T cells (and vice versa) [1/d]
    fTadeact = 0.01,                # fraction of activated CD8+ T-cells that deactivate to resting or "exhausted" CD8+ T-cells
    fTap = 3.77,                    # Ratio of the partition coeff for activated CD8+ T-cells to that of resting CD8+ T-cells
    kTaapop = 0.06,                 # apoptosis rate of CD69+ CD8+ T-cells (activated CD8+ T-cells) [1/d]
    fTaprolif = 2.11,               # proliferation rate of activated CD8+ T-cells [1/d]
    fTrapop = 0.2,                  # Ratio of the apoptosis rate of resting CD8+ T-cells to that of CD69+CD8+ T-cells
    kTrexit = 0.05,                 # trafficking rate of resting/ post-activated CD8+ T-cells from PB to tissue [1/d]
    fdrug = 6.38,                   # The TDB injection effect on traficking of resting or "exhausted" CD8+ T-cells from PB to tissues
    fa0 = 0.,                       # Ratio of the placebo/TDB injection effect on traficking of activated CD8+ T-cells from PB to tissues to that of resting CD8+ T-cells
    fAICD = 1.5,                    # Effect of activation-induced cell death (AICD) on CD69+CD8+ T-cells
    fTa0deact = 0,                  # Conversion rate of post-activated CD8+ T-cells to resting CD8+ T-cells [1/d]
    fTa0apop = 2.02,                # Ratio of the apoptosis rate of post-activated CD8+ T-cells to that of CD69+CD8+ T-cells
    finj = 1.,                      # The placebo injection effect on traficking of resting or "exhausted" CD8+ T-cells from PB to tissues
    #
    # PB Physiology
    Trpbref_perml = 2.0e6,          # baseline conc of resting CD8+ T-cells in PB [cells/mL]
    Bpbref_perml = 1.0e6,           # baseline conc of CD19+CD20+ B-cells in PB [cells/mL]
    Vpb = 380.,                     # physiological volume of peripheral blood [mL]
    # 
    # Spleen Physiology 
    # Trtiss_perml = 1e9            # Baseline conc. of resting CD8+ T-cells in spleen [cells/ml]
    # Btiss_perml = 9e8             # Baseline conc. of CD19+CD20+ B-cells in spleen [cells/ml]
    Vtissue = 7.0,                  # physiological volume of spleen [mL]
    KTrp = 500.0,                   # Partition coeff: ratio of the baseline conc. of resting CD8+ T-cells in spleen to that in PB          
    KBp = 900. ,                    # Partition coeff: ratio of the baseline conc. of B-cells in spleen to that in PB
    Kp = 0.14,                      # Partition coeff: ratio of the drug conc. in spleen to that in PB at any given time
    # 
    # LNs Physiology
    # Trtiss2_perml = 1e9           # Baseline conc. of resting CD8+ T-cells in LNs [cells/ml]
    # Btiss2_perml = 6e8            # Baseline conc. of CD19+CD20+ B-cells in LNs [cells/ml]
    Vtissue2 = 25. ,                # total physiological volume of LNs [mL]
    KTrp2 = 500.,                   # Partition coeff: ratio of the baseline conc. of B-cells in LNs to that in PB
    KBp2 = 600.,                    # Partition coeff: ratio of the baseline conc. of resting CD8+ T-cells in LNs to that in PB
    Kp2 = 0.07,                     # Partition coeff: ratio of the drug conc. in LNs to that in PB at any given time
    # 
    # BM Physiology
    # Trtiss3_perml = 1e8           # Baseline conc. of resting CD8+ T-cells in BM [cells/ml]
    # B1920tiss3_perml = 6e7        # Baseline conc. of CD19+CD20+ B-cells in BM [cells/ml]
    # B19no20tiss3_perml = 1.5e7    # Baseline conc. of CD19+CD20- B-cells in BM [cells/ml]
    # B19tiss3_perml = 7.5e7        # Baseline conc. of CD19+ B cells in BM [cells/ml]
    Vtissue3 = 50.0,                # Physiology volume of BM [mL]
    KTrp3 = 50.0,                   # Partition coeff: ratio of the baseline conc. of resting CD8+ T-cells in BM to that in PB
    KBp3 = 60.0,                    # Partition coeff: ratio of the baseline conc. of CD19+CD20+ B-cells in BM to that in PB
    Kp3 = 0.07,                     # Partition coeff: ratio of the drug conc. in BM to that in PB at any given time
    kBmat_kBapop_ratio = 0.25,      # The ratio of maturation rate of CD19+CD20- B-cells in BM to B-cell apoptosis rate
    B19no20_B1920_ratio = 0.25,     # CD19+CD20- B-cell:CD19+CD20+ B-cell ratio in BM
    # 
    # Tumor Physiology (THIS PART IS TURNED OFF IN MONKEY)
    tumor_on = 0.0,                 # this turns on/ off the tumor compartment 
    # Trtumor_perml                 # Baseline conc. of resting CD8+ T-cells in tumor [cells/ml]
    # Btumor_perml                  # Baseline conc. of CD19+CD20+ B-cells in tumor [cells/ml]
    Vtumor = 1.0,                   # Physiological volume of tumor [ml]
    KTrptumor = 1.0,                # Partition coeff: ratio of the baseline conc. of resting CD8+ T-cells in tumor to that in PB
    KBptumor = 1.0,                 # Partition coeff: ratio of the baseline conc. of B-cells in tumor to that in PB
    Kptumor = 1.0,                  # Partition coeff: ratio of the drug conc. in tumor to that in PB at any given time
    kBtumorprolif = 0.025,          # proliferation rate of B-cells in tumor [1/d]
    # 
    # Cytokine 
    kIL6prod = 4.0,                 # production rate of IL6 [pg/d]
    IL6_tiss_contribution = 0.,   # Contribution of tissue-produced IL6 to the systematic IL6 conc
    thalfIL6 = 20.0,                # [minute]
    #
    # Init-related values
    Trpbo_perml = 2.0e6,    # [1/mL]
    Bpbo_perml = 1.0e6,     # [1/mL]
    BAFFo = 1.0,
    # 
    # RTX 
    Vc_rtx = 28.0,
    Vp_rtx = 9.0,
    Cl_rtx = 8.0,
    Cld_rtx = 23.0,
    Vm_rtx = 0.0,
    Km_rtx = 1.0,
    Kmkill_rtx = 25.0,
    # 
    # not documented in Supp 2 
    Vm_tdb = 212 ,          # [μg/d]
    Km_tdb = 4.74 ,         # [μg/d]
    act0on = 1.0,
    tissue1on = 1.0,
    tissue2on = 0.0,
    tissue3on = 0.0,
    Bcell_tumor_trafficking_on = 1.0,
    depleteTpb = 0.0,
    depleteBpb = 0.0,
    BW = 3.4 ,  # [kg]
    Ktgen = 1.0,   # [1/d]
    fBprolif = 0.0,
    fBexit = 1.0,
    kTgen = 1.0 ,  # [1/d]
    kBAFFprod = 1.0,
    thBAFF = 30.0,   # [minute]
    fBAFFo = 0.1,
    tinjhalf = 1.0,   # [d]
    fTgenbl = 0.0,
    kBtiss3exit = 0.03,
    fapop_v24 = 0.0,
    fBtissue3_v1 = 1.0,
    kBapop_cll = 1.0,
    kBgen_cll = 1.0,
    kabs_TDB = 1.4,
    fbio_TDB = 0.6,
    PKflag = 1.0,
    VPid = 1.0,
    end_time = 1.0,
    fvalidation = 0.0,
);


# human params; only update those that are different from cyno 
# based on Supp Table 2, Hosseini et al., 2020
# https://www.nature.com/articles/s41540-020-00145-7

p_homo =  deepcopy(p_cyno);
p_homo.Cl_tdb = 5.4
p_homo.Vpb = 5000.              # [ml] PB volume; Jones 2019 used 3.126L
p_homo.Vtissue = 210.           # [ml] spleen volume; Jones 2019 used 221mL
p_homo.KTrp = 200.
p_homo.KBp = 333.
p_homo.Vtissue2 = 400.          # [ml] lymph nodes volume; Jones 2019 used 274mL
p_homo.KTrp2 = 190.
p_homo.KBp2 = 190.
p_homo.Vtissue3 = 500.          # [ml] bone marrow volume; Jones 2019 used 1500mL
p_homo.KTrp3 = 60.
p_homo.KBp3 = 80.

p_homo.BW = 75. # [kg] this value was presumed 
