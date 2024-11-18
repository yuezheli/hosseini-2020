using DataFrames, DataFramesMeta, CSV
using DifferentialEquations
using Plots, PDFmerger
using Random, Statistics

include("Hosseini_params.jl");
include("Hosseini_init.jl");
include("Hosseini_ode.jl")

Hosseini_params.end_time = 80.0

u0 = Hosseini_init(Hosseini_params);

prob = ODEProblem(Hosseini_ode!,u0, (0.0, Hosseini_params.end_time), Hosseini_params);

sol = solve(prob, alg=AutoTsit5(Rosenbrock23()), reltol=1e-9, abstol=1e-9, saveat = 1.);

#=
using MAT

vars = matread("data/MatlabResults.mat")["data"];
df_names = matread("data/MatlabNames.mat")["names"];
mt = matread("data/MatlabTime.mat")["time"][:,1];
df_names = Symbol.(df_names)[:,1];
mdf = DataFrame(vars,:auto);
rename!(mdf, df_names);

# Get rid of repeated 0 time point in MATLAB data (i.e remove the before dose time)
mt = mt[2:end];
mdf = mdf[2:end,:];
mdf.Time .= mt;

# Interpolate the julia simulation solution at the same time points
sol_mtime = sol(mt);
mdf.Source .= "MATLAB";

# convert julia solultion into dataframe
sdf = DataFrame(sol_mtime);
rename!(sdf, [:value1, :value2, :value3, :value4, :value5, :value6, :value7, :value8, :value9, :value10, :value11, :value12, :value13, :value14, :value15, :value16, 
        :value17, :value18, :value19, :value20, :value21, :value22, :value23, :value24, :value25, :value26, :value27, :value28, :value29, :value30, 
        :value31, :value32, :value33, :value34, :value35, :value36] .=> [:actTtiss, :Btiss, :TDBc_ugperkg, :actTpb, :restTtiss, :TDBp_ugperkg, :restTpb,
        :Bpb, :BAFF, :act0Tpb, :act0Ttiss, :injection_effect, :Btiss2, :act0Ttiss2, :restTtiss2, 
        :actTtiss2, :RTXc_ugperkg, :RTXp_ugperkg, :drug_effect, :B1920tiss3, :B19no20tiss3, 
        :restTtiss3, :act0Ttiss3, :actTtiss3, :Blinc_ug, :restTtumor, :actTtumor, :Btumor, :act0Ttumor, :TDBsc_ugperkg, 
        :TDBc_ugperml_AUC, :IL6pb, :IL6tiss, :IL6tiss2, :IL6tiss3, :IL6tumor] );
select!(sdf, Not(:timestamp));
sdf.Time .= mt;
sdf.Source .= "Julia";

mdf_long = stack(mdf, Not([:Source,:Time]));
sdf_long = stack(sdf, Not([:Source,:Time]));

df_long = vcat(mdf_long, sdf_long);

# Plot everything, MATLAB comparison
using AlgebraOfGraphics, CairoMakie
layer = visual(Lines,linewidth=6.0, alpha = 0.5);
xy = data(df_long) * mapping(:Time, :value, color=:Source,layout=:variable, linestyle = :Source);
# with_theme(theme_default())
fg = draw(layer * xy, axis=(; titlesize=32,xlabelsize=48,ylabelsize=48,xticklabelsize=32,yticklabelsize=32,xlabel="Time (days)",
                      ylabel="Value (Mixed Units)"),figure=(; resolution=(3840, 2160)),
                      legend = (; labelsize=32,titlesize=48),
                      facet = (; linkyaxes = :none))
save("validation-MATLAB.png", fg)
=#


#------------------------------------- REPLACE PK_v26 with UNIT CONVERSION -------------------------------------#
sdf0 = DataFrame(sol);
rename!(sdf0, [:value1, :value2, :value3, :value4, :value5, :value6, :value7, :value8, :value9, :value10, :value11, :value12, :value13, :value14, :value15, :value16, 
        :value17, :value18, :value19, :value20, :value21, :value22, :value23, :value24, :value25, :value26, :value27, :value28, :value29, :value30, 
        :value31, :value32, :value33, :value34, :value35, :value36] .=> [:actTtiss, :Btiss, :TDBc_ugperkg, :actTpb, :restTtiss, :TDBp_ugperkg, :restTpb,
        :Bpb, :BAFF, :act0Tpb, :act0Ttiss, :injection_effect, :Btiss2, :act0Ttiss2, :restTtiss2, 
        :actTtiss2, :RTXc_ugperkg, :RTXp_ugperkg, :drug_effect, :B1920tiss3, :B19no20tiss3, 
        :restTtiss3, :act0Ttiss3, :actTtiss3, :Blinc_ug, :restTtumor, :actTtumor, :Btumor, :act0Ttumor, :TDBsc_ugperkg, 
        :TDBc_ugperml_AUC, :IL6pb, :IL6tiss, :IL6tiss2, :IL6tiss3, :IL6tumor] );

include("Hosseini_ode_2.jl");

prob2 = ODEProblem(Hosseini_ode_2!,u0, (0.0, Hosseini_params.end_time), Hosseini_params);
sol2 = solve(prob2, alg=AutoTsit5(Rosenbrock23()), reltol=1e-9, abstol=1e-9, saveat = 1.);
sdf2 = DataFrame(sol2);
rename!(sdf2, [:value1, :value2, :value3, :value4, :value5, :value6, :value7, :value8, :value9, :value10, :value11, :value12, :value13, :value14, :value15, :value16, 
        :value17, :value18, :value19, :value20, :value21, :value22, :value23, :value24, :value25, :value26, :value27, :value28, :value29, :value30, 
        :value31, :value32, :value33, :value34, :value35, :value36] .=> [:actTtiss, :Btiss, :TDBc_ugperkg, :actTpb, :restTtiss, :TDBp_ugperkg, :restTpb,
        :Bpb, :BAFF, :act0Tpb, :act0Ttiss, :injection_effect, :Btiss2, :act0Ttiss2, :restTtiss2, 
        :actTtiss2, :RTXc_ugperkg, :RTXp_ugperkg, :drug_effect, :B1920tiss3, :B19no20tiss3, 
        :restTtiss3, :act0Ttiss3, :actTtiss3, :Blinc_ug, :restTtumor, :actTtumor, :Btumor, :act0Ttumor, :TDBsc_ugperkg, 
        :TDBc_ugperml_AUC, :IL6pb, :IL6tiss, :IL6tiss2, :IL6tiss3, :IL6tumor] );


Plots.plot(title = "T-cell dependent bispecific, PB");
Plots.plot!(sdf0.timestamp, sdf0.TDBc_ugperkg, label = "PK_v26", alpha = 0.5, linestyle = :dash, lw = 4);
Plots.plot!(sdf2.timestamp, sdf2.TDBc_ugperkg, label = "unit conversion", alpha = 0.5, linestyle = :solid, lw = 2);
Plots.xlabel!("Time (day)"); Plots.ylabel!("T-cell dependent bispecific (ug/kg)")

# check additional 
Plots.plot(title = "CD8+ T-cells, PB");
Plots.plot!(sdf0.timestamp, sdf0.actTpb/ Hosseini_params.Vpb * 1e-3, label = false);
Plots.xlabel!("Time (day)"); Plots.ylabel!("CD8+ T-cells (/uL)")

Plots.plot(title = "CD8+ T-cells, PB");
Plots.plot!(sdf0.timestamp, sdf0.actTpb ./ (sdf.actTpb .+ sdf.restTpb .+ sdf.act0Tpb ) * 100, label = false);
Plots.xlabel!("Time (day)"); Plots.ylabel!("CD69+ CD8+ T-cells (%)")

Plots.plot(title = "B-cells, PB");
Plots.plot!(sdf0.timestamp, sdf0.Bpb/ Hosseini_params.Vpb * 1e-3, label = false);
Plots.xlabel!("Time (day)"); Plots.ylabel!("B-cells (/uL)")

Plots.plot(title = "IL6, PB");
Plots.plot!(sdf0.timestamp, sdf0.IL6pb, label = false);
Plots.xlabel!("Time (day)"); Plots.ylabel!("IL6 (pg/mL)")

#------------------------------------- FULL BLOW 2-comp PK -------------------------------------#
include("Hosseini_ode_homo.jl");

p_cyno.end_time = 40.  

p_cyno_0 = ComponentArray(p_cyno, 
        CL_mosun = Hosseini_params.Cl_tdb * p_cyno.BW,         # [mL/d]
        Q_mosun  = Hosseini_params.Cld_tdb * p_cyno.BW,        # [mL/d]
        V1_mosun = Hosseini_params.Vc_tdb * p_cyno.BW,         # [mL]
        V2_mosun = Hosseini_params.Vp_tdb * p_cyno.BW,         # [mL]
        infusion_mosun = 0.0)

p_cyno_0.Cl_tdb = Hosseini_params.Cl_tdb

Dose = u0.TDBc_ugperkg * p_cyno.BW / p_cyno_0.V1_mosun         # [ug/mL]

# initialization 
u02 = ComponentArray(actTtiss=0., actTpb=0., act0Tpb=0., act0Ttiss=0., act0Ttiss2=0., actTtiss2=0., act0Ttiss3=0., actTtiss3=0., act0Ttumor=0., actTtumor=0., 
IL6pb=0., IL6tiss=0., IL6tiss2=0., IL6tiss3=0., IL6tumor=0., 
injection_effect=0., drug_effect=0., RTXc_ugperkg=0., RTXp_ugperkg=0., Blinc_ug=0., TDBc_ugperml=Dose, TDBp_ugperml=0. ,
Bpb = (1-p_cyno.depleteBpb) * p_cyno.Bpbo_perml * p_cyno.Vpb, 
Btiss = p_cyno.KBp * p_cyno.Bpbo_perml * p_cyno.Vtissue, 
Btiss2 = p_cyno.KBp2 * p_cyno.Bpbo_perml * p_cyno.Vtissue2, 
B1920tiss3 = p_cyno.Bpbref_perml * p_cyno.KBp3 * p_cyno.Vtissue3, 
B19no20tiss3 = p_cyno.Bpbref_perml * p_cyno.KBp3 * p_cyno.B19no20_B1920_ratio * p_cyno.Vtissue3, 
Btumor = p_cyno.KBptumor * p_cyno.Bpbo_perml * p_cyno.Vtumor, 
restTpb = (1-p_cyno.depleteTpb) * p_cyno.Trpbo_perml * p_cyno.Vpb,
restTtiss = p_cyno.KTrp * p_cyno.Trpbo_perml * p_cyno.Vtissue, 
restTtiss2 = p_cyno.KTrp2 * p_cyno.Trpbo_perml * p_cyno.Vtissue2, 
restTtiss3 = p_cyno.KTrp3 * p_cyno.Trpbo_perml * p_cyno.Vtissue3,
restTtumor = p_cyno.KTrptumor * p_cyno.Trpbo_perml * p_cyno.Vtumor,  
BAFF = p_cyno.BAFFo );

sol02 = solve(ODEProblem(Hosseini_ode_homo!, u02, (0., Hosseini_params.end_time), p_cyno_0), alg=AutoTsit5(Rosenbrock23()), reltol=1e-18, saveat = 1.);
sdf02 = DataFrame(sol02);
rename!(sdf02, [:value1, :value2, :value3, :value4, :value5, :value6, :value7, :value8, :value9, :value10, :value11, :value12, :value13, :value14, :value15, 
        :value16, :value17, :value18, :value19, :value20, :value21, :value22, :value23, :value24, :value25, :value26, :value27, :value28, :value29, :value30, 
        :value31, :value32, :value33, :value34,] .=> [:actTtiss, :actTpb, :act0Tpb, :act0Ttiss, :act0Ttiss2, :actTtiss2, :act0Ttiss3, :actTtiss3, :act0Ttumor, 
        :actTtumor, :IL6pb, :IL6tiss, :IL6tiss2, :IL6tiss3, :IL6tumor, :injection_effect, :drug_effect, :RTXc_ugperkg, :RTXp_ugperkg, :Blinc_ug, :TDBc_ugperml, 
        :TDBp_ugperml, :Bpb, :Btiss, :Btiss2, :B1920tiss3, :B19no20tiss3, :Btumor, :restTpb, :restTtiss, :restTtiss2, :restTtiss3, :restTtumor, :BAFF] );


plot(title = "T-cell dependent bispecific, PB", titlefontsize = 8, xticks = (0:7:70, string.(0:7:70)) );
plot!(sdf0.timestamp, sdf0.TDBc_ugperkg * p_cyno_0.BW / p_cyno_0.V1_mosun,label = "MATLAB", alpha = 0.5, lw = 4, linestyle = :dash);
plot!(sdf02.timestamp, sdf02.TDBc_ugperml, label = "2-comp PK", alpha = 0.8, lw = 2, linestyle = :solid);
xlabel!("Time (day)"); xlims!(0,42); ylabel!("T-cell dependent bispecific (ug/mL)")

#------------------------------------- CYNO PARAMS from Supp Table 2 vs. MATLAB PARAMS -------------------------------------#
# purpose of this chunk of code: to update params listed in Supp Table 2, and check simulation results
p_cyno_2 = ComponentArray(p_cyno, 
        CL_mosun = p_cyno.Cl_tdb * p_cyno.BW,         # [mL/d]
        Q_mosun  = p_cyno.Cld_tdb * p_cyno.BW,        # [mL/d]
        V1_mosun = p_cyno.Vc_tdb * p_cyno.BW,         # [mL]
        V2_mosun = p_cyno.Vp_tdb * p_cyno.BW,         # [mL]
        infusion_mosun = 0.0);

u03 = ComponentArray(actTtiss=0., actTpb=0., act0Tpb=0., act0Ttiss=0., act0Ttiss2=0., actTtiss2=0., act0Ttiss3=0., actTtiss3=0., act0Ttumor=0., actTtumor=0., 
IL6pb=0., IL6tiss=0., IL6tiss2=0., IL6tiss3=0., IL6tumor=0., 
injection_effect=0., drug_effect=0., RTXc_ugperkg=0., RTXp_ugperkg=0., Blinc_ug=0., TDBc_ugperml=u0.TDBc_ugperkg * p_cyno.BW / p_cyno_2.V1_mosun, TDBp_ugperml=0. ,
Bpb = (1-p_cyno.depleteBpb) * p_cyno.Bpbo_perml * p_cyno.Vpb, 
Btiss = p_cyno.KBp * p_cyno.Bpbo_perml * p_cyno.Vtissue, 
Btiss2 = p_cyno.KBp2 * p_cyno.Bpbo_perml * p_cyno.Vtissue2, 
B1920tiss3 = p_cyno.Bpbref_perml * p_cyno.KBp3 * p_cyno.Vtissue3, 
B19no20tiss3 = p_cyno.Bpbref_perml * p_cyno.KBp3 * p_cyno.B19no20_B1920_ratio * p_cyno.Vtissue3, 
Btumor = p_cyno.KBptumor * p_cyno.Bpbo_perml * p_cyno.Vtumor, 
restTpb = (1-p_cyno.depleteTpb) * p_cyno.Trpbo_perml * p_cyno.Vpb,
restTtiss = p_cyno.KTrp * p_cyno.Trpbo_perml * p_cyno.Vtissue, 
restTtiss2 = p_cyno.KTrp2 * p_cyno.Trpbo_perml * p_cyno.Vtissue2, 
restTtiss3 = p_cyno.KTrp3 * p_cyno.Trpbo_perml * p_cyno.Vtissue3,
restTtumor = p_cyno.KTrptumor * p_cyno.Trpbo_perml * p_cyno.Vtumor,  
BAFF = p_cyno.BAFFo );

sol3 = solve(ODEProblem(Hosseini_ode_homo!, u03, (0., p_cyno_2.end_time), p_cyno_2), alg=AutoTsit5(Rosenbrock23()), reltol=1e-18, saveat = 0.1);
sdf3 = DataFrame(sol3);
rename!(sdf3, [:value1, :value2, :value3, :value4, :value5, :value6, :value7, :value8, :value9, :value10, :value11, :value12, :value13, :value14, :value15, 
        :value16, :value17, :value18, :value19, :value20, :value21, :value22, :value23, :value24, :value25, :value26, :value27, :value28, :value29, :value30, 
        :value31, :value32, :value33, :value34,] .=> [:actTtiss, :actTpb, :act0Tpb, :act0Ttiss, :act0Ttiss2, :actTtiss2, :act0Ttiss3, :actTtiss3, :act0Ttumor, 
        :actTtumor, :IL6pb, :IL6tiss, :IL6tiss2, :IL6tiss3, :IL6tumor, :injection_effect, :drug_effect, :RTXc_ugperkg, :RTXp_ugperkg, :Blinc_ug, :TDBc_ugperml, 
        :TDBp_ugperml, :Bpb, :Btiss, :Btiss2, :B1920tiss3, :B19no20tiss3, :Btumor, :restTpb, :restTtiss, :restTtiss2, :restTtiss3, :restTtumor, :BAFF] );


plot(title = "T-cell dependent bispecific, PB", titlefontsize = 8, xticks = (0:7:70, string.(0:7:70)) );
plot!(sdf02.timestamp, sdf02.TDBc_ugperml, label = "MATLAB params", alpha = 0.5, lw = 4, linestyle = :dash);
plot!(sdf3.timestamp, sdf3.TDBc_ugperml, label = "cyno params", alpha = 0.7);
xlabel!("Time (day)"); xlims!(0,42); ylabel!("T-cell dependent bispecific (ug/mL)")

