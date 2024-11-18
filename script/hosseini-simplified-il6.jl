# date: 9/20/24
# author: Yuezhe Li 
# purpose of this code: to test peak IL6 concentration after dosing 

using Pkg; Pkg.activate("")

using DifferentialEquations, ModelingToolkit 
using ModelingToolkit: getdefault
using Plots
using DataFrames

# observed data; Fig 6 of Hosseini et al., 2020; https://pubmed.ncbi.nlm.nih.gov/32859946/
obs = DataFrame(
    dose_mg = [0.05, 0.2, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 
    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 
    1.2, 1.2, 1.2, 1.2, 1.2, 
    1.6, 1.6, 1.6, 1.6, 1.6, 
    2., 2., 2., 
    2.8, 2.8, 2.8, 2.8, 2.8, 2.8, 2.8], 
    IL6_pgmL = [4.73, 12.84, 258.61, 116.78, 42.86, 33.65, 23.01, 13.69, 9.04, 1307.66, 702.17, 447.81, 432.46, 316.80, 248.77, 170.16, 195.26, 94.53, 112.33, 74.17, 64.59, 41.23, 32.36, 29.16, 22.90, 16.77, 16.78, 11.47, 9.32, 6.16, 5.36, 3.79, 
    2344.45, 1258.01, 194.44, 175.36, 119.86, 90.87, 36.99, 33.34, 21.28, 11.03, 3.29, 
    246.57, 193.56, 66.33, 20.46, 16.63, 
    64214.28, 1644.08, 87.09, 35.43, 17.14, 
    1241.60, 18.94, 4.75, 
    423.22, 332.33, 113.81, 44.76, 21.65, 19.51, 12.44]
)

# original hosseini model 
include("model/tdb_cyno_2.jl");
@mtkbuild tdb = TDB_homo();
# simplified model (T cells not in spleen & BM, B cells not in spleen)
include("model/tdb_cyno_6.jl");
@mtkbuild tdb6 = TDB_homo6();

function il6_after_first_dose(dose_mg; sys = tdb, injection_effect = 10., tspan = (0., 7.), turnontumor = true, Trpbo_perml = 2E6, Bpbo_perml = 1E6)
    if turnontumor
        prob = ODEProblem(sys, [], tspan, [sys.tumor_on => 1., sys.Trpbo_perml => Trpbo_perml, sys.Bpbo_perml => Bpbo_perml]);
    else
        prob = ODEProblem(sys, [], tspan, [sys.Trpbo_perml => Trpbo_perml, sys.Bpbo_perml => Bpbo_perml]);
    end
    # add dosing
    prob_new = remake(prob, u0 = Dict([sys.TDBc_ugperml => dose_mg*1E3/getdefault(sys.V1_TDB), sys.TCEinjection_effect => injection_effect]))
    @time sol_new = solve(prob_new, reltol=1e-18, saveat = 0.1); 
    il6 = DataFrame(t = sol_new.t, IL6 = sol_new[:IL6pb]);
    return il6
end

IL6_dose = [0.05, 0.2, 0.4, 0.8, 1., 1.2, 1.6, 2., 2.8]; # [mg]
IL6_peak_orig = [];
IL6_peak_simp = [];

for il6dose in IL6_dose
    sim1 = il6_after_first_dose(il6dose; sys = tdb, Bpbo_perml = 0.5E6); 
    append!(IL6_peak_orig, maximum(sim1.IL6)); 
    sim2 = il6_after_first_dose(il6dose; sys = tdb6, Bpbo_perml = 0.5E6); 
    append!(IL6_peak_simp, maximum(sim2.IL6)); 
end

pIL6peak = plot(ylims = (1., 1E5), yaxis = :log10, xlabel = "Dose level (mg)", ylabel = "First IL6 peak (pg/mL)", legend = :outerright, dpi = 300);
plot!(IL6_dose, IL6_peak_orig, label = "Hosseini 2020", alpha = 0.7, lw = 3);
plot!(IL6_dose, IL6_peak_simp, label = "Simplified", alpha = 0.7, lw = 3);
plot!(obs.dose_mg, obs.IL6_pgmL, label = "Observed", seriestype = :scatter, alpha = 0.5)
plot!(xticks = (IL6_dose, string.(IL6_dose)), xrotation = 45);
display(pIL6peak)

savefig(pIL6peak, "deliv/figure/hosseini-simplification/IL6.png");

## sensitivity analysis, alt init blood T cell count, original Hosseini model 
IL6_peak_Tscan_1e6 = [];
IL6_peak_Tscan_point1e6 = [];

for il6dose in IL6_dose
    sim_1 = il6_after_first_dose(il6dose; sys = tdb, Bpbo_perml = 0.5E6, Trpbo_perml = 1E6); 
    append!(IL6_peak_Tscan_1e6, maximum(sim_1.IL6)); 
    sim_point1 = il6_after_first_dose(il6dose; sys = tdb, Bpbo_perml = 0.5E6, Trpbo_perml = 0.1E6); 
    append!(IL6_peak_Tscan_point1e6, maximum(sim_point1.IL6)); 
end

pIL6peak_Tscan = plot(ylims = (1., 1E5), yaxis = :log10, xlabel = "Dose level (mg)", ylabel = "First IL6 peak (pg/mL)", legend = :outerright, dpi = 300);
plot!(IL6_dose, IL6_peak_orig, label = "pb T cell = 2E6/mL", alpha = 0.7, lw = 3);
plot!(IL6_dose, IL6_peak_Tscan_1e6, label = "pb T cell = 1E6/mL", alpha = 0.7, lw = 3);
plot!(IL6_dose, IL6_peak_Tscan_point1e6, label = "pb T cell = 0.1E6/mL", alpha = 0.7, lw = 3);
plot!(obs.dose_mg, obs.IL6_pgmL, label = "Observed", seriestype = :scatter, alpha = 0.5);
plot!(xticks = (IL6_dose, string.(IL6_dose)), xrotation = 45);
display(pIL6peak_Tscan);

## sensitivity analysis, alt init blood B cell count, original Hosseini model 
IL6_peak_Bscan_1e5 = [];
IL6_peak_Bscan_1e4 = [];

for il6dose in IL6_dose
    sim_1 = il6_after_first_dose(il6dose; sys = tdb, Bpbo_perml = 1E5); 
    append!(IL6_peak_Bscan_1e5, maximum(sim_1.IL6)); 
    sim_2 = il6_after_first_dose(il6dose; sys = tdb, Bpbo_perml = 1E4); 
    append!(IL6_peak_Bscan_1e4, maximum(sim_2.IL6)); 
end

pIL6peak_Bscan = plot(ylims = (1., 1E5), yaxis = :log10, xlabel = "Dose level (mg)", ylabel = "First IL6 peak (pg/mL)", legend = :outerright, dpi = 300);
plot!(IL6_dose, IL6_peak_orig, label = "pb B cell = 5E5/mL", alpha = 0.7, lw = 3);
plot!(IL6_dose, IL6_peak_Bscan_1e5, label = "pb B cell = 1E5/mL", alpha = 0.7, lw = 3);
plot!(IL6_dose, IL6_peak_Bscan_1e4, label = "pb B cell = 1E4/mL", alpha = 0.7, lw = 3);
plot!(obs.dose_mg, obs.IL6_pgmL, label = "Observed", seriestype = :scatter, alpha = 0.5);
plot!(xticks = (IL6_dose, string.(IL6_dose)), xrotation = 45);
display(pIL6peak_Bscan);

savefig(pIL6peak_Tscan, "deliv/figure/hosseini-simplification/IL6-sens-pb-T.png");
savefig(pIL6peak_Bscan, "deliv/figure/hosseini-simplification/IL6-sens-pb-B.png");

## sensitivity analysis, alt init blood T cell count, simplified Hosseini model 
IL6_2_peak_Tscan_1e6 = [];
IL6_2_peak_Tscan_point1e6 = [];

for il6dose in IL6_dose
    sim_1 = il6_after_first_dose(il6dose; sys = tdb6, Bpbo_perml = 0.5E6, Trpbo_perml = 1E6); 
    append!(IL6_2_peak_Tscan_1e6, maximum(sim_1.IL6)); 
    sim_point1 = il6_after_first_dose(il6dose; sys = tdb6, Bpbo_perml = 0.5E6, Trpbo_perml = 0.1E6); 
    append!(IL6_2_peak_Tscan_point1e6, maximum(sim_point1.IL6)); 
end

pIL6_2_peak_Tscan = plot(ylims = (1., 1E5), yaxis = :log10, xlabel = "Dose level (mg)", ylabel = "First IL6 peak (pg/mL)", legend = :outerright, dpi = 300);
plot!(IL6_dose, IL6_peak_simp, label = "pb T cell = 2E6/mL", alpha = 0.7, lw = 3);
plot!(IL6_dose, IL6_2_peak_Tscan_1e6, label = "pb T cell = 1E6/mL", alpha = 0.7, lw = 3);
plot!(IL6_dose, IL6_2_peak_Tscan_point1e6, label = "pb T cell = 0.1E6/mL", alpha = 0.7, lw = 3);
plot!(obs.dose_mg, obs.IL6_pgmL, label = "Observed", seriestype = :scatter, alpha = 0.5);
plot!(xticks = (IL6_dose, string.(IL6_dose)), xrotation = 45);
display(pIL6_2_peak_Tscan);

## sensitivity analysis, alt init blood B cell count, simplified Hosseini model 
IL6_2_peak_Bscan_1e5 = [];
IL6_2_peak_Bscan_1e4 = [];

for il6dose in IL6_dose
    sim_1 = il6_after_first_dose(il6dose; sys = tdb6, Bpbo_perml = 1E5); 
    append!(IL6_2_peak_Bscan_1e5, maximum(sim_1.IL6)); 
    sim_2 = il6_after_first_dose(il6dose; sys = tdb6, Bpbo_perml = 1E4); 
    append!(IL6_2_peak_Bscan_1e4, maximum(sim_2.IL6)); 
end

pIL6_2_peak_Bscan = plot(ylims = (1., 1E5), yaxis = :log10, xlabel = "Dose level (mg)", ylabel = "First IL6 peak (pg/mL)", legend = :outerright, dpi = 300);
plot!(IL6_dose, IL6_peak_simp, label = "pb B cell = 5E5/mL", alpha = 0.7, lw = 3);
plot!(IL6_dose, IL6_2_peak_Bscan_1e5, label = "pb B cell = 1E5/mL", alpha = 0.7, lw = 3);
plot!(IL6_dose, IL6_2_peak_Bscan_1e4, label = "pb B cell = 1E4/mL", alpha = 0.7, lw = 3);
plot!(obs.dose_mg, obs.IL6_pgmL, label = "Observed", seriestype = :scatter, alpha = 0.5);
plot!(xticks = (IL6_dose, string.(IL6_dose)), xrotation = 45);
display(pIL6_2_peak_Bscan);

savefig(pIL6_2_peak_Tscan, "deliv/figure/hosseini-simplification/IL6-v6-sens-pb-T.png");
savefig(pIL6_2_peak_Bscan, "deliv/figure/hosseini-simplification/IL6-v6-sens-pb-B.png");
