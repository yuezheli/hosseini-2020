# date: 11/1/2024
# purpose of this code: to capture blinatumomab PK and IL6 production 

using Pkg; Pkg.activate(".")

using DifferentialEquations
using DataFrames, CSV
using Plots
using ComponentArrays

include("blin_ode.jl");
include("init.jl");
include("param_blin.jl");

BSA = 1.7*1.7; 

function blinatumomab_cIV(dose_amt, dose_time; BSA = 1.7*1.7)
    global cbs = [];
    if length(dose_time) >= 1
        for i in 1:length(dose_time)
            function affect!(integrator)
                integrator.p.blin_infusion = dose_amt[i] * BSA;
            end
            cb = PresetTimeCallback(dose_time[i],affect!);
            global cbs = push!(cbs, cb);
        end
    end
    cbset = CallbackSet(cbs...);
    return cbset
end

prob0 = ODEProblem(blin_ode!, u0_0, (0, 21.), p_homo_blin); 

p_blin_1 = deepcopy(p_homo_blin); p_blin_1.blin_infusion = 5*BSA; 
prob1 = remake(prob0, p = p_blin_1);
cb1 = blinatumomab_cIV([15, 60], [7, 14]);
sol1 = solve(prob1, saveat = 0.2, reltol = 1E-18, callback = cb1);
sdf1 = DataFrame(sol1); rename!(sdf1, [ Symbol(names(sdf1)[i+1]) => keys(sol1.u[1])[i] for i in 1:length(sol1.u[1])]);

p_blin_2 = deepcopy(p_homo_blin); p_blin_2.blin_infusion = 60*BSA; 
prob2 = remake(prob0, p = p_blin_2);
sol2 = solve(prob2, saveat = 0.2, reltol = 1E-18);
sdf2 = DataFrame(sol2); rename!(sdf2, [ Symbol(names(sdf2)[i+1]) => keys(sol2.u[1])[i] for i in 1:length(sol2.u[1])]);

p_pk = plot(xlabel = "Time (day)", ylabel = "Blincyto amount (ug)", legend = :bottomleft, yaxis = :log10, ylims = [0.1, 100], xticks = [0, 7, 14, 21]);
plot!(xguidefontsize=8,yguidefontsize=8, dpi = 300);
plot!(sdf1.timestamp, sdf1.Blinc_ug, label = "5/15/60 μg/m2/day week 1/2/3"); 
plot!(sdf2.timestamp, sdf2.Blinc_ug, label = "60 μg/m2/day week 1/2/3"); 

p_actT = plot(xlabel = "Time (day)", ylabel = "Plasma activated T cell (#/uL)", legend = :bottomleft, yaxis = :log10, ylims = [1, 1000], xticks = [0, 7, 14, 21]);
plot!(xguidefontsize=8,yguidefontsize=8, dpi = 300);
plot!(sdf1.timestamp, sdf1.actTpb/(p_homo_blin.Vpb*1E3), label = "5/15/60 μg/m2/day week 1/2/3"); 
plot!(sdf2.timestamp, sdf2.actTpb/(p_homo_blin.Vpb*1E3), label = "60 μg/m2/day week 1/2/3"); 

p_B = plot(xlabel = "Time (day)", ylabel = "Plasma B cell (#/uL)", legend = :topright, yaxis = :log10, ylims = [0.01, 1000], xticks = [0, 7, 14, 21]);
plot!(xguidefontsize=8,yguidefontsize=8, dpi = 300);
plot!(sdf1.timestamp, sdf1.Bpb/(p_homo_blin.Vpb*1E3), label = "5/15/60 μg/m2/day week 1/2/3"); 
plot!(sdf2.timestamp, sdf2.Bpb/(p_homo_blin.Vpb*1E3), label = "60 μg/m2/day week 1/2/3"); 

p_cd19neg_B = plot(xlabel = "Time (day)", ylabel = "CD19+CD20- B cell in BM", legend = :topright, yaxis = :log10, ylims = [1E7, 1E10], xticks = [0, 7, 14, 21]);
plot!(xguidefontsize=8,yguidefontsize=8, dpi = 300);
plot!(sdf1.timestamp, sdf1.B19no20tiss3, label = "5/15/60 μg/m2/day week 1/2/3"); 
plot!(sdf2.timestamp, sdf2.B19no20tiss3, label = "60 μg/m2/day week 1/2/3"); 

savefig(plot(p_pk, p_actT, p_B, p_cd19neg_B), "../figure/blincyto-pk.png")

# read on data from Hosseini 2020 Fig 4
obs_il6_1 = CSV.read("../data/hosseini-blincyto-il6/il6-dose-up.csv", DataFrame, header=true);
obs_il6_2 = CSV.read("../data/hosseini-blincyto-il6/il6-60ug-m2-day.csv", DataFrame, header=true);

p_IL6_1 = plot(xlabel = "Time (day)", ylabel = "Plasma IL6 (pg/mL)", legend = :bottomleft, yaxis = :log10, ylims = [1E-1, 1E4], xticks = [0, 7, 14, 21]);
plot!(xguidefontsize=8,yguidefontsize=8, dpi = 300);
plot!(sdf1.timestamp, sdf1.IL6pb, label = "5/15/60 μg/m2/day week 1/2/3"); 
plot!(obs_il6_1.time_day, obs_il6_1.il6_pgmL, label  = "obs", seriestype = :scatter, alpha = 0.7);

p_IL6_2 = plot(xlabel = "Time (day)", ylabel = "Plasma IL6 (pg/mL)", legend = :bottomleft, yaxis = :log10, ylims = [1E-1, 1E4], xticks = [0, 7, 14, 21]);
plot!(xguidefontsize=8,yguidefontsize=8, dpi = 300);
plot!(sdf2.timestamp, sdf2.IL6pb, label = "60 μg/m2/day week 1/2/3"); 
plot!(obs_il6_2.time_day, obs_il6_2.il6_pgmL, label  = "obs", seriestype = :scatter, alpha = 0.7);

savefig(plot(p_IL6_1, p_IL6_2), "../figure/blincyto-IL6.png")
