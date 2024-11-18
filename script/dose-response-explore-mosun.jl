# date: 9/4/24
# author: Yuezhe Li 
# purpose of this code: to explore dose-response in mosunetuzumab simulation 

using Pkg; Pkg.activate("");

using DifferentialEquations
using DataFrames, CSV, DataFramesMeta
using ProgressMeter
using Plots
using PDFmerger: append_pdf!

include("model/param.jl");
include("model/init.jl");
include("model/bs-lonca-3.jl");
include("model/dosing-function-TDB-iv.jl");
include("script/helper.jl"); 
include("model/ParamUpdate.jl");

p_mosun = TDB_param_update("mosunetuzumab"); 

# TDB iv dosing
function IndSims(; subject_dose_amt = [1., 2., 60., 60., 30.] * 1E3,  tumor_vol = 49.5, subject_dose_time = 7. *[0 1 2 3 6], tspan = (0., 63.), p_homo = p_homo_3)
    u0_eql = preequlibrium3(tumor_vol, 0., p_homo); 
    u0_eql.TDBc_ugperml = subject_dose_amt[1] / p_homo.V1_TDB
    prob0 = ODEProblem(bs_lonca_ode_3!, u0_eql, tspan, p_homo); 
    iv_tdb_dosing = tdb_dosing_iv(subject_dose_amt, subject_dose_time, 10., p_homo);
    sol0 = solve(prob0, alg=AutoTsit5(Rosenbrock23(autodiff = false)), reltol=1e-12, saveat = 1., callback = iv_tdb_dosing);
    sdf0 = DataFrame(sol0); rename!(sdf0, Symbol.(names(sdf0)[2:end]) .=> collect(keys(sol0.u[end])) );
    @transform!(sdf0, :Btotal = :Btumor .+ :Btumor_trans .+ :Btumor_cd19neg .+ :Btumor_cd19neg_trans);
    @transform!(sdf0, :Vtumor = :Btotal * Vc / 0.375)
    @transform!(sdf0, :ntv = :Vtumor / first(:Vtumor) )
    @transform!(sdf0, :TaBratio_tumor = :actTtumor ./ (:Btumor_cd19neg .+ :Btumor) )
    @transform!(sdf0, :BTrratio_tumor = (:Btumor_cd19neg .+ :Btumor) ./ (:restTtumor .+ :act0Ttumor) )
    @transform!(sdf0, :TaBratio_pb = :actTpb ./ :Bpb )
    @transform!(sdf0, :BTrratio_pb = :Bpb./ (:restTpb .+ :act0Tpb) )
    return sdf0
end 

function TDBSims(max_dose; tumor_vol = 49.5, kBprolif = 0.025, p_homo = p_homo_3)
    subject_dose_amt = [1., 2., max_dose, max_dose, max_dose] * 1E3
    p_tdb_backup = deepcopy(p_homo);
    p_tdb_backup.kBprolif = kBprolif; 
    df = IndSims(subject_dose_amt = subject_dose_amt, tumor_vol = tumor_vol, p_homo = p_tdb_backup)
    df.Dose .= string(max_dose) * "mg"
    return df
end

function TDBSimsBatch(; tumor_vol = 49.5, kBprolif = 0.025, p_homo = p_homo_3)
    df60 = TDBSims(60., tumor_vol = tumor_vol, kBprolif = kBprolif, p_homo = p_homo); 
    df30 = TDBSims(30., tumor_vol = tumor_vol, kBprolif = kBprolif, p_homo = p_homo); 
    df1 = TDBSims(1., tumor_vol = tumor_vol, kBprolif = kBprolif, p_homo = p_homo); 
    dfpoint1 = TDBSims(0.1, tumor_vol = tumor_vol, kBprolif = kBprolif, p_homo = p_homo); 

    # pk check
    p_pk_tdb = plot(legend = :outerright, xticks = ([0, 21, 42, 63], string.([0, 21, 42, 63])));
    plot!(df60.timestamp, df60.TDBc_ugperml, label = "60mg");
    plot!(df30.timestamp, df30.TDBc_ugperml, label = "30mg");
    plot!(df1.timestamp, df1.TDBc_ugperml, label = "1mg");
    plot!(dfpoint1.timestamp, dfpoint1.TDBc_ugperml, label = "0.1mg");
    xlabel!("Time (day)", xguidefontsize = 8); 
    ylabel!("TDB plasma conc (ug/mL)", yguidefontsize = 8);
    plot!(yscale = :log10); 
    ylims!(0.01, 100);

    # pd
    pd_title = "init tv = " * string(tumor_vol) * "mL"
    p_pd_tv = plot(legend = :outerright, xticks = ([0, 21, 42, 63], string.([0, 21, 42, 63])), title = pd_title);
    plot!(df60.timestamp, df60.ntv, label = "60mg");
    plot!(df30.timestamp, df30.ntv, label = "30mg");
    plot!(df1.timestamp, df1.ntv, label = "1mg");
    plot!(dfpoint1.timestamp, dfpoint1.ntv, label = "0.1mg");
    xlabel!("Time (day)", xguidefontsize = 8); 
    ylabel!("Normalized tumor volume", yguidefontsize = 8);

    # activated T cells  
    p_actT_pb = plot(legend = :outerright, xticks = ([0, 21, 42, 63], string.([0, 21, 42, 63])));
    plot!(df60.timestamp, df60.actTpb, label = "60mg");
    plot!(df30.timestamp, df30.actTpb, label = "30mg");
    plot!(df1.timestamp, df1.actTpb, label = "1mg");
    plot!(dfpoint1.timestamp, dfpoint1.actTpb, label = "0.1mg");
    xlabel!("Time (day)", xguidefontsize = 8); 
    ylabel!("Activated T cell count in plasma", yguidefontsize = 8);

    p_actT_tm = plot(legend = :outerright, xticks = ([0, 21, 42, 63], string.([0, 21, 42, 63])));
    plot!(df60.timestamp, df60.actTtumor, label = "60mg");
    plot!(df30.timestamp, df30.actTtumor, label = "30mg");
    plot!(df1.timestamp, df1.actTtumor, label = "1mg");
    plot!(dfpoint1.timestamp, dfpoint1.actTtumor, label = "0.1mg");
    xlabel!("Time (day)", xguidefontsize = 8); 
    ylabel!("Activated T cell count in tumor", yguidefontsize = 8);

    # save figure
    filename = "deliv/figure/dose-response-mosun/" * "initTV=" * string(tumor_vol) * "mL" * "-kBprolif=" * string(kBprolif) * ".pdf"
    savefig(p_pk_tdb, filename);
    append_pdf!(filename, savefig(p_pd_tv, "temp.pdf"), cleanup=true)
    append_pdf!(filename, savefig(p_actT_pb, "temp.pdf"), cleanup=true)
    append_pdf!(filename, savefig(p_actT_tm, "temp.pdf"), cleanup=true)

    # save results 
    df = vcat(df60, df30, df1, dfpoint1);
    csvfilename =  "data/dose-response-mosun/" * "initTV=" * string(tumor_vol) * "mL" * "-kBprolif=" * string(kBprolif) * ".csv"
    CSV.write(csvfilename, df);
end

TDBSimsBatch(tumor_vol = 1., kBprolif = 0.1, p_homo = p_mosun)
TDBSimsBatch(tumor_vol = 10., kBprolif = 0.1, p_homo = p_mosun)
TDBSimsBatch(tumor_vol = 1., kBprolif = 0.025, p_homo = p_mosun)
