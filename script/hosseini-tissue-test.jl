# date: 9/12/24
# author: Yuezhe Li 
# purpose of this code: to test if other organs can be turned off

using Pkg; Pkg.activate("")

using DifferentialEquations, ModelingToolkit 
using ModelingToolkit: getdefault
using Plots

# build model 
include("model/tdb_cyno_2.jl");
@mtkbuild tdb = TDB_homo();
prob0 = ODEProblem(tdb, [tdb.TDBc_ugperml => 1E3/getdefault(tdb.V1_TDB), tdb.TCEinjection_effect => 10.], (0., 42.), [tdb.tumor_on => 1.]);
sol0 = solve(prob0, reltol=1e-18, saveat = 0.2); 

# check B cells leaving BM 
# getdefault(tdb.kBtiss3exit)*getdefault(tdb.Vpb)*max.(0,sol0[:B1920tiss3]/getdefault(tdb.Vtissue3) .- sol0[:Bpb]/getdefault(tdb.Vpb)*getdefault(tdb.KBp3))[end]
# 2.055E9

#======================== on and off of all organs ========================#
prob1 = remake(prob0, p = Dict([tdb.tissue1on => 0., tdb.tissue2on => 0., tdb.tissue3on => 0.]));
sol1 = solve(prob1, reltol=1e-18, saveat = 0.2); 

prob2 = remake(prob0, p = Dict([tdb.tissue1on => 0., tdb.tissue2on => 0.]));
sol2 = solve(prob2, reltol=1e-18, saveat = 0.2); 

prob3 = remake(prob0, p = Dict([tdb.tissue1on => 0., tdb.tissue3on => 0.]));
sol3 = solve(prob3, reltol=1e-18, saveat = 0.2); 

prob4 = remake(prob0, p = Dict([tdb.tissue2on => 0., tdb.tissue3on => 0.]));
sol4 = solve(prob4, reltol=1e-18, saveat = 0.2); 

prob5 = remake(prob0, p = Dict([tdb.tissue1on => 0.]));
sol5 = solve(prob5, reltol=1e-18, saveat = 0.2); 

prob6 = remake(prob0, p = Dict([tdb.tissue2on => 0.]));
sol6 = solve(prob6, reltol=1e-18, saveat = 0.2); 

prob7 = remake(prob0, p = Dict([tdb.tissue3on => 0.]));
sol7 = solve(prob7, reltol=1e-18, saveat = 0.2); 

p_comp = plot(title = "");
plot!(sol0.t, sol0[:Btumor], label = "all on", linestyle = :dash, lw = 2);
plot!(sol1.t, sol1[:Btumor], label = "PB only", linestyle = :dash, lw = 1.5);
plot!(sol2.t, sol2[:Btumor], label = "Spleen & LN off");
plot!(sol3.t, sol3[:Btumor], label = "Spleen & BM off");
plot!(sol4.t, sol4[:Btumor], label = "LN & BM off");
plot!(sol5.t, sol5[:Btumor], label = "spleen off", linestyle = :dashdotdot);
plot!(sol6.t, sol6[:Btumor], label = "LN off", linestyle = :dashdotdot);
plot!(sol7.t, sol7[:Btumor], label = "BM off", linestyle = :dashdotdot);
display(p_comp);

#======================== remove T cells but in PB ========================#
include("model/tdb_cyno_3.jl");
@mtkbuild tdb3 = TDB_homo3();
prob3 = ODEProblem(tdb3, [tdb3.TDBc_ugperml => 1E3/getdefault(tdb3.V1_TDB), tdb3.TCEinjection_effect => 10.], (0., 42.), [tdb3.tumor_on => 1.]);
sol_3 = solve(prob3, reltol=1e-18, saveat = 0.2); 

p3 = plot(sol0.t, sol0[:Btumor], label = "all on", linestyle = :dash, lw = 2);
plot!(sol_3.t, sol_3[:Btumor], label = "T cells only in PB & tumor")

savefig(p3, "figure/hosseini-simplification/T-in-PB-TME.png");

#======================== remove T cells in spleen ========================#
include("model/tdb_cyno_4.jl");
@mtkbuild tdb4 = TDB_homo4();
prob4 = ODEProblem(tdb4, [tdb4.TDBc_ugperml => 1E3/getdefault(tdb4.V1_TDB), tdb4.TCEinjection_effect => 10.], (0., 42.), [tdb4.tumor_on => 1.]);
sol_4 = solve(prob4, reltol=1e-18, saveat = 0.2); 

plot(sol0.t, sol0[:Btumor], label = "all on", linestyle = :dash, lw = 2);
plot!(sol_4.t, sol_4[:Btumor], label = "T cells not in spleen")

#======================== remove T cells in spleen & BM ========================#
include("model/tdb_cyno_5.jl");
@mtkbuild tdb5 = TDB_homo5();

prob5 = ODEProblem(tdb5, [tdb5.TDBc_ugperml => 1E3/getdefault(tdb5.V1_TDB), tdb5.TCEinjection_effect => 10.], (0., 42.), [tdb5.tumor_on => 1.]);
sol_5 = solve(prob5, reltol=1e-18, saveat = 0.2); 

plot(sol0.t, sol0[:Btumor], label = "all on", linestyle = :dash, lw = 2);
plot!(sol_5.t, sol_5[:Btumor], label = "T cells not in spleen & BM")

#======================== remove T cells in spleen & BM, remove B cells in spleen ========================#
include("model/tdb_cyno_6.jl");
@mtkbuild tdb6 = TDB_homo6();

prob6 = ODEProblem(tdb6, [tdb6.TDBc_ugperml => 1E3/getdefault(tdb6.V1_TDB), tdb6.TCEinjection_effect => 10.], (0., 42.), [tdb6.tumor_on => 1.]);
sol_6 = solve(prob6, reltol=1e-18, saveat = 0.2); 

plot(sol0.t, sol0[:Btumor], label = "all on", linestyle = :dash, lw = 2);
plot!(sol_6.t, sol_6[:Btumor], label = "T cells not in spleen & BM, B cells not in spleen")

#======================== remove T cells in spleen & BM, remove B cells in spleen & LN ========================#
include("model/tdb_cyno_7.jl");
@mtkbuild tdb7 = TDB_homo7();

prob7 = ODEProblem(tdb7, [tdb7.TDBc_ugperml => 1E3/getdefault(tdb7.V1_TDB), tdb7.TCEinjection_effect => 10.], (0., 42.), [tdb7.tumor_on => 1.]);
sol_7 = solve(prob7, reltol=1e-18, saveat = 0.2); 

plot(sol0.t, sol0[:Btumor], label = "all on", linestyle = :dash, lw = 2);
plot!(sol_7.t, sol_7[:Btumor], label = "T cells not in spleen & BM, B cells not in spleen & LN")

#======================== remove T cells in spleen & BM, remove B cells in spleen; replace B cell output from bone marrow with steady state value ========================#
include("model/tdb_cyno_8.jl");
@mtkbuild tdb8 = TDB_homo8();

prob8 = ODEProblem(tdb8, [tdb8.TDBc_ugperml => 1E3/getdefault(tdb8.V1_TDB), tdb8.TCEinjection_effect => 10.], (0., 42.), [tdb8.tumor_on => 1.]);
sol_8 = solve(prob8, reltol=1e-18, saveat = 0.2); 

plot(sol0.t, sol0[:Btumor], label = "all on", linestyle = :dash, lw = 2);
plot!(sol_8.t, sol_8[:Btumor], label = "T cells not in spleen & BM, B cells not in spleen, steady state B cell output")

#======================== additional visualization ========================#
using DataFrames
using DataFramesMeta
sdf0 = DataFrame(sol0);
sdf6 = DataFrame(sol_6);

@select!(sdf0, :timestamp, :TDBc_ugperml, :restTpb, :actTpb, :act0Tpb, :act0Ttiss2, :restTtiss2, :actTtiss2, :restTtumor, :actTtumor, :act0Ttumor, :IL6pb, :Bpb, :Btiss2, :B1920tiss3, :B19no20tiss3, :Btumor); 
sdf0.Model .= "Hosseini";
@select!(sdf6, :timestamp, :TDBc_ugperml, :restTpb, :actTpb, :act0Tpb, :act0Ttiss2, :restTtiss2, :actTtiss2, :restTtumor, :actTtumor, :act0Ttumor, :IL6pb, :Bpb, :Btiss2, :B1920tiss3, :B19no20tiss3, :Btumor); 
sdf6.Model .= "no T cells in BM & spleen, no B cells in the spleen";

df0 = stack(sdf0, Not([:Model,:timestamp]));
df6 = stack(sdf6, Not([:Model,:timestamp]));

df_long = vcat(df0, df6);

using AlgebraOfGraphics, CairoMakie
layer = visual(Lines,linewidth=6.0, alpha = 0.5);
xy = data(df_long) * mapping(:timestamp, :value, color=:Model,layout=:variable, linestyle = :Model);
# with_theme(theme_default())
fg = draw(layer * xy, axis=(; titlesize=32,xlabelsize=48,ylabelsize=48,xticklabelsize=32,yticklabelsize=32,xlabel="Time (days)",
                      ylabel="Value (Mixed Units)"),figure=(; resolution=(3840, 2160)),
                      legend = (; labelsize=32,titlesize=48),
                      facet = (; linkyaxes = :none));
save("figure/hosseini-simplification/no_T_in_BM_spleen_no_B_in_spleen.png", fg);

#======================== additional dosing schemes ========================#

function tdb_iv_dosing(dose_amt, dose_time, V1_TDB; injection_effect_init = 10.)
    global cbs_human = [];
    if length(dose_time) > 1
        for i in 1:length(dose_time)
            function affect!(integrator)
                integrator[:TDBc_ugperml] += dose_amt[i]/ V1_TDB;
                integrator[:TCEinjection_effect] += injection_effect_init;
            end
            cb = PresetTimeCallback(dose_time[i],affect!);
            global cbs_human = push!(cbs_human, cb);
        end
    end
    cbset_human = CallbackSet(cbs_human...);
    return cbset_human
end

# dosing event 
dose_amt = [2., 60., 60., 30.] .* 1e3;  # [ug]
dose_time = [7, 14, 21, 42];                 # [day] 
cbset_human = tdb_iv_dosing(dose_amt, dose_time, getdefault(tdb.V1_TDB)); 

prob0 = ODEProblem(tdb, [tdb.TDBc_ugperml => 1E3/getdefault(tdb.V1_TDB), tdb.TCEinjection_effect => 10.], (0., 84.), [tdb.tumor_on => 1.]);
sol0 = solve(prob0, reltol=1e-18, saveat = 0.2, callback = cbset_human); 

prob6 = ODEProblem(tdb6, [tdb6.TDBc_ugperml => 1E3/getdefault(tdb6.V1_TDB), tdb6.TCEinjection_effect => 10.], (0., 84.), [tdb6.tumor_on => 1.]);
sol_6 = solve(prob6, reltol=1e-18, saveat = 0.2, callback = cbset_human); 

sdf0 = DataFrame(sol0);
sdf6 = DataFrame(sol_6);

@select!(sdf0, :timestamp, :TDBc_ugperml, :restTpb, :actTpb, :act0Tpb, :act0Ttiss2, :restTtiss2, :actTtiss2, :restTtumor, :actTtumor, :act0Ttumor, :IL6pb, :Bpb, :Btiss2, :B1920tiss3, :B19no20tiss3, :Btumor); 
sdf0.Model .= "Hosseini";
@select!(sdf6, :timestamp, :TDBc_ugperml, :restTpb, :actTpb, :act0Tpb, :act0Ttiss2, :restTtiss2, :actTtiss2, :restTtumor, :actTtumor, :act0Ttumor, :IL6pb, :Bpb, :Btiss2, :B1920tiss3, :B19no20tiss3, :Btumor); 
sdf6.Model .= "no T cells in BM & spleen, no B cells in the spleen";

df0 = stack(sdf0, Not([:Model,:timestamp]));
df6 = stack(sdf6, Not([:Model,:timestamp]));

df_long = vcat(df0, df6);

layer = visual(Lines,linewidth=6.0, alpha = 0.5);
xy = data(df_long) * mapping(:timestamp, :value, color=:Model,layout=:variable, linestyle = :Model);
# with_theme(theme_default())
fg = draw(layer * xy, axis=(; titlesize=32,xlabelsize=48,ylabelsize=48,xticklabelsize=32,yticklabelsize=32,xlabel="Time (days)",
                      ylabel="Value (Mixed Units)"),figure=(; resolution=(3840, 2160)),
                      legend = (; labelsize=32,titlesize=48),
                      facet = (; linkyaxes = :none));
save("figure/hosseini-simplification/no_T_in_BM_spleen_no_B_in_spleen-label-dosing.png", fg);


#======================== additional comparison, more simplified model ========================#
prob8 = ODEProblem(tdb8, [tdb8.TDBc_ugperml => 1E3/getdefault(tdb8.V1_TDB), tdb8.TCEinjection_effect => 10.], (0., 84.), [tdb8.tumor_on => 1.]);
sol_8 = solve(prob8, reltol=1e-18, saveat = 0.2, callback = cbset_human); 

sdf0 = DataFrame(sol0);
sdf8 = DataFrame(sol_8);

@select!(sdf0, :timestamp, :TDBc_ugperml, :restTpb, :actTpb, :act0Tpb, :act0Ttiss2, :restTtiss2, :actTtiss2, :restTtumor, :actTtumor, :act0Ttumor, :IL6pb, :Bpb, :Btiss2, :Btumor); 
sdf0.Model .= "Hosseini";
@select!(sdf8, :timestamp, :TDBc_ugperml, :restTpb, :actTpb, :act0Tpb, :act0Ttiss2, :restTtiss2, :actTtiss2, :restTtumor, :actTtumor, :act0Ttumor, :IL6pb, :Bpb, :Btiss2, :Btumor); 
sdf8.Model .= "no T cells in BM & spleen, no B cells in the spleen; steady state B cell output from BM";

df0 = stack(sdf0, Not([:Model,:timestamp]));
df8 = stack(sdf8, Not([:Model,:timestamp]));

df_long = vcat(df0, df8);

layer = visual(Lines,linewidth=6.0, alpha = 0.5);
xy = data(df_long) * mapping(:timestamp, :value, color=:Model,layout=:variable, linestyle = :Model);
# with_theme(theme_default())
fg = draw(layer * xy, axis=(; titlesize=32,xlabelsize=48,ylabelsize=48,xticklabelsize=32,yticklabelsize=32,xlabel="Time (days)",
                      ylabel="Value (Mixed Units)"),figure=(; resolution=(3840, 2160)),
                      legend = (; labelsize=32,titlesize=48),
                      facet = (; linkyaxes = :none));
save("figure/hosseini-simplification/no_T_in_BM_spleen_no_B_in_spleen_BM-label-dosing.png", fg);

