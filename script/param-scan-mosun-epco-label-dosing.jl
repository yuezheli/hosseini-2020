# date: 8/9/24
# author: Yuezhe Li 
# purpose of this code: pop simulation for mosun, glofit, and epco monotherapy under their label

using Pkg; Pkg.activate("");

using DifferentialEquations
using DifferentialEquations.EnsembleAnalysis
using DataFrames, CSV
using QuasiMonteCarlo
using ProgressMeter

include("model/param.jl");
include("model/ParamUpdate.jl");
include("model/init.jl");
include("model/bs-lonca-3.jl");
include("script/helper.jl"); 
include("model/dosing-function-TDB-sc.jl");
include("model/dosing-function-TDB-iv.jl");


# simulation setup (use the TDB that is desired)
numSim = 10
# TDB = "mosunetuzumab"
# TDB = "glofitamab"
TDB = "epcoritamab"

# define function for tumor volume post-processing
function ReturnOutput(sol, i)
    sdf0 = DataFrame(sol);
    rename!(sdf0, Symbol.(names(sdf0)[2:end]) .=> collect(keys(sol.u[end])) );
    return sdf0, false
end

# placeholders
tumor_vol = 7.5; #  # [mL]
cd19neg_fraction = 0.5
u0_eql = preequlibrium3(tumor_vol, cd19neg_fraction, p_homo_3, 10.); 
tspan = (0., 12*7.);
prob0 = ODEProblem(bs_lonca_ode_3!, u0_eql, tspan, p_homo_3);
global df =  DataFrame(ID =[], tumor_prof=[], k_trans = [], lonca_kill_neg = [], cd19neg_frac = [], init_tv = [], 
            timestamp=[], restTtumor=[], actTtumor=[], act0Ttumor=[], normalized_tv=[], Btumor_tot = [], Bpb=[], restTpb=[], actTpb=[], act0Tpb=[], 
            TDBc_ugperml=[], LONCAc_ugperml=[]);
# param scan set up the same as population change 

# tumor proliferation rate (kBtumorprolif) range obtained from Susilo et al., 2023 # https://ascpt.onlinelibrary.wiley.com/doi/10.1111/cts.13501
# the rate for dying B cells to be removed from the system (k_trans) was calibrated based on comparison to prior PBPK model simulation outcome
# maximum rate of lonca-induced CD19-/low cell-death (lonca_kill_neg) was calibrated based on comparison to prior PBPK model simulation outcome
# CD19-/low fraction was assumed to be between 0 and 0.1
# tumor volume was assumed to be between 0.01mL and 50mL (sample on exp scale)
lb_hlc = [log10(0.00001), 0.01, 0., 0., log10(0.1)]; 
ub_hlc = [log10(0.15), 0.05, 0.07, 0.1, log10(100.)];
s = QuasiMonteCarlo.sample(numSim, lb_hlc, ub_hlc, LatinHypercubeSample());
param_ranges =  [(s[1,i], s[2,i], s[3,i], s[4,i], s[5,i]) for i in axes(s,2)];

if TDB == "mosunetuzumab"
    p_mosun = TDB_param_update(TDB);
    # IV dosing of mosun
    cbset = tdb_dosing_iv([1., 2., 60., 60., 30., 30.] * 1E3, 7. *[0 1 2 3 6 9], 10., p_mosun.V1_mosun);
    function prob_func(prob,i, repeat)
        p_mosun_backup = deepcopy(p_mosun);
        p_mosun_backup.kBtumorprolif = 10^param_ranges[i][1];
        p_mosun_backup.k_trans = param_ranges[i][2];#
        p_mosun_backup.lonca_kill_neg = param_ranges[i][3];
        u0_new = preequlibrium3(10^param_ranges[i][5], param_ranges[i][4], p_mosun_backup);
        # use IV dosing of mosun  
        u0_new.TDBc_ugperml = 1E3/p_mosun.V1_mosun
        remake(prob,p=p_mosun_backup, u0 = u0_new);
    end
    ensemble_prob = EnsembleProblem(prob0, prob_func=prob_func, output_func=ReturnOutput);
    csvfilename = "data/sims/mono-mosun-pop-size-" * string.(numSim) * ".csv"
elseif TDB == "glofitamab"
    p_glo = TDB_param_update(TDB);
    # Glofit dosing (label and LOTIS-7 are the same)
    cbset  =  tdb_dosing_iv([2.5, 10., 30., 30., 30.] * 1E3, 7. *[1 2 3 6 9], 10., p_glo.V1_mosun); 
    function prob_func(prob,i, repeat)
        p_glo_backup = deepcopy(p_glo);
        p_glo_backup.kBtumorprolif = 10^param_ranges[i][1];
        p_glo_backup.k_trans = param_ranges[i][2];
        p_glo_backup.lonca_kill_neg = param_ranges[i][3];
        u0_new = preequlibrium3(10^param_ranges[i][5], param_ranges[i][4], p_glo_backup);
        u0_new.Bpb = 0.
        u0_new.Btiss = 0.
        u0_new.Btiss2 = 0.
        u0_new.B1920tiss3 = 0.
        remake(prob,p=p_glo_backup, u0 = u0_new);
    end
    ensemble_prob = EnsembleProblem(prob0, prob_func=prob_func, output_func=ReturnOutput);
    csvfilename = "data/sims/mono-glofit-pop-size-" * string.(numSim) * ".csv"
elseif TDB == "epcoritamab"
    p_epco = TDB_param_update(TDB);
    # Epco dosing (label)
    cbset = tdb_dosing_sc([0.16, 0.8, 48., 48., 48., 48., 48., 48., 48., 48., 48., 48.] * 1E3, 7. *[0 1 2 3 4 5 6 7 8 9 10 11], 10., p_epco);
    function prob_func(prob,i, repeat)
        p_epco_backup = deepcopy(p_epco);
        p_epco_backup.kBtumorprolif = 10^param_ranges[i][1];
        p_epco_backup.k_trans = param_ranges[i][2];
        p_epco_backup.lonca_kill_neg = param_ranges[i][3];
        u0_new = preequlibrium3(10^param_ranges[i][5], param_ranges[i][4], p_epco_backup); 
        u0_new.TDBdepot_ug = 0.16 * 1E3 * p_epco_backup.fbio_TDB;
        remake(prob,p=p_epco_backup, u0 = u0_new);
    end
    ensemble_prob = EnsembleProblem(prob0, prob_func=prob_func, output_func=ReturnOutput); 
    csvfilename = "data/sims/mono-epco-pop-size-" * string.(numSim) * ".csv"   
else
    println("Need TDB choice")
end

# run simulation 
try
    sim = solve(ensemble_prob, alg=AutoTsit5(Rosenbrock23(autodiff = false)), EnsembleThreads(), trajectories = numSim, reltol = 1E-18, callback=cbset, saveat = 1.); 
    @showprogress for k = 1:length(sim)
        ntv = (sim[k].Btumor .+ sim[k].Btumor_trans .+ sim[k].Btumor_cd19neg .+ sim[k].Btumor_cd19neg_trans) / (sim[k].Btumor[1] + sim[k].Btumor_cd19neg[1]);
        tmpdf = DataFrame(
            timestamp = sim[k].timestamp,
            restTtumor = sim[k].restTtumor,
            actTtumor = sim[k].actTtumor,
            act0Ttumor = sim[k].act0Ttumor, 
            normalized_tv = ntv, 
            Btumor_tot = max.(sim[k].Btumor .+ sim[k].Btumor_trans .+ sim[k].Btumor_cd19neg .+ sim[k].Btumor_cd19neg_trans, 1.0), 
            Bpb = sim[k].Bpb, 
            restTpb = sim[k].restTpb, 
            actTpb = sim[k].actTpb, 
            act0Tpb = sim[k].act0Tpb, 
            TDBc_ugperml = sim[k].TDBc_ugperml, 
            LONCAc_ugperml = sim[k].LONCAc_ugperml, 
        );
        tmpdf.ID .= k; 
        tmpdf.tumor_prof .= 10^param_ranges[k][1];
        tmpdf.k_trans .= param_ranges[k][2];
        tmpdf.lonca_kill_neg .= param_ranges[k][3];
        tmpdf.cd19neg_frac .= param_ranges[k][4];
        tmpdf.init_tv .= 10^param_ranges[k][5];
        df = vcat(df, tmpdf);
    end
    CSV.write(csvfilename, df)
catch 
    println("check the code")
end
