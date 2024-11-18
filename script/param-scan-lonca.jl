# author: Yuezhe Li 
# date: 8/7/24
# purpose of this code: to test param scan in glofit simulation 

using DifferentialEquations
using DifferentialEquations.EnsembleAnalysis
using DataFrames, CSV
using QuasiMonteCarlo
using ProgressMeter

RunSim = true
# number of virtual patients used in the simulation; can be changed 
numSim = Int(2E4)

include("../model/param.jl");
include("../model/init.jl");
include("../model/bs-lonca-3.jl");
include("helper.jl"); 
include("../model/dosing-function.jl");

p_mosun = deepcopy(p_homo_3);

# define function for tumor volume post-processing
function ReturnOutput(sol, i)
    sdf0 = DataFrame(sol);
    rename!(sdf0, Symbol.(names(sdf0)[2:end]) .=> collect(keys(sol.u[end])) );
    return sdf0, false
end

tspan = (0., 105.);
subject_dose_time_lonca = 7. *[0 3 6 9 12];

#dosescheme = "150ugkg"
#subject_dose_amt_lonca = p_mosun.BW * [150., 150., 75., 75., 75.]; # [ug]

#dosescheme = "120ugkg"
#subject_dose_amt_lonca = p_mosun.BW * [120., 120., 75., 75., 75.]; # [ug]

dosescheme = "90ugkg"
subject_dose_amt_lonca = p_mosun.BW * [90., 90., 90., 90., 90.]; # [ug]

tumor_vol = 1; #  # [mL]
cd19neg_fraction = 0.5

u0_eql = deepcopy(u0_0); 
u0_eql.LONCAc_ugperml = subject_dose_amt_lonca[1] / p_mosun.V1_lonca; 
u0_eql.Btumor = tumor_vol * 0.375/ Vc * (1 - cd19neg_fraction);
u0_eql.Btumor_cd19neg = tumor_vol * 0.375/ Vc * cd19neg_fraction; 

# Lonca dosing 
global cbs = [];
if length(subject_dose_time_lonca) > 1
    for i in 2:length(subject_dose_time_lonca)
        function affect_lonca!(integrator)
            integrator.u.LONCAc_ugperml += subject_dose_amt_lonca[i]/p_mosun.V1_lonca;
        end
        cb = PresetTimeCallback(subject_dose_time_lonca[i],affect_lonca!);
        global cbs = push!(cbs, cb);
    end
end
cbset = CallbackSet(cbs...);

prob0 = ODEProblem(bs_lonca_ode_3!, u0_eql, tspan, p_mosun);

# tumor proliferation rate (kBtumorprolif) range obtained from Susilo et al., 2023 # https://ascpt.onlinelibrary.wiley.com/doi/10.1111/cts.13501
# the rate for dying B cells to be removed from the system (k_trans) was calibrated based on comparison to prior PBPK model simulation outcome
# maximum rate of lonca-induced CD19-/low cell-death (lonca_kill_neg) was calibrated based on comparison to prior PBPK model simulation outcome
# CD19-/low fraction was assumed to be between 0 and 0.1
# tumor volume was assumed to be between 0.01mL and 50mL (sample on exp scale)
lb_hlc = [log10(0.00001), 0.01, 0., 0., log10(0.1)]; 
ub_hlc = [log10(0.15), 0.05, 0.07, 0.1, log10(100.)];
s = QuasiMonteCarlo.sample(numSim, lb_hlc, ub_hlc, LatinHypercubeSample());
param_ranges =  [(s[1,i], s[2,i], s[3,i], s[4,i], s[5,i]) for i in axes(s,2)];

function prob_func(prob,i, repeat)
        p_mosun_backup = deepcopy(p_mosun);
        p_mosun_backup.kBtumorprolif = 10^param_ranges[i][1];
        p_mosun_backup.k_trans = param_ranges[i][2];#
        p_mosun_backup.lonca_kill_neg = param_ranges[i][3];
        u0_new = preequlibrium3(10^param_ranges[i][5], param_ranges[i][4], p_mosun); 
        u0_new.LONCAc_ugperml = subject_dose_amt_lonca[1]/p_mosun.V1_lonca
        remake(prob,p=p_mosun_backup, u0 = u0_new);
end
ensemble_prob = EnsembleProblem(prob0, prob_func=prob_func, output_func=ReturnOutput);

global df =  DataFrame(ID =[], tumor_prof=[], k_trans = [], lonca_kill_neg = [], cd19neg_frac = [], init_tv = [], 
            timestamp=[], normalized_tv=[], Btumor_tot = [], Bpb=[], TDBc_ugperml=[], LONCAc_ugperml=[]);

if RunSim 
    @time sim = solve(ensemble_prob, alg=AutoTsit5(Rosenbrock23()), EnsembleThreads(), trajectories = numSim, reltol=1e-18, callback=cbset, saveat = 1.); 
    @showprogress for k = 1:length(sim)
        tmpdf = DataFrame(
            timestamp = sim[k].timestamp, 
            normalized_tv = (sim[k].Btumor .+ sim[k].Btumor_trans .+ sim[k].Btumor_cd19neg .+ sim[k].Btumor_cd19neg_trans) / (sim[k].Btumor[1] + sim[k].Btumor_cd19neg[1]), 
            Btumor_tot = max.(sim[k].Btumor .+ sim[k].Btumor_trans .+ sim[k].Btumor_cd19neg .+ sim[k].Btumor_cd19neg_trans, 1.0), 
            Bpb = sim[k].Bpb, 
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
    csvfilename = "../data/sims/ADCT-402-Lonca-mono-" * dosescheme * "-pop-size-" * string.(numSim) * ".csv"
    CSV.write(csvfilename, df);
end

