# author: Yuezhe Li
# date: June 21, 2023
# purpose: to recap cyno simulation from Hosseini et al., 2020
# https://www.nature.com/articles/s41540-020-00145-7

using Pkg; Pkg.activate()

using DataFrames, DataFramesMeta, CSV
using DifferentialEquations
using Plots, PDFmerger

include("Hosseini_params.jl");
include("Hosseini_ode_homo.jl");  # two-compartmental PK for mosun

# cyno params from Supp Table 2
p_cyno_2 = ComponentArray(p_cyno, 
        CL_mosun = p_cyno.Cl_tdb * p_cyno.BW,         # [mL/d]
        Q_mosun  = p_cyno.Cld_tdb * p_cyno.BW,        # [mL/d]
        V1_mosun = p_cyno.Vc_tdb * p_cyno.BW,         # [mL]
        V2_mosun = p_cyno.Vp_tdb * p_cyno.BW,         # [mL]
        infusion_mosun = 0.0);

p_cyno_2.tissue3on = 1.
p_cyno_2.tissue2on = 1.
p_cyno_2.act0on = 1.
injection_effect_init = 10.0

TDBc_ugperkg = 10.  # [ug]
subject_dose_time = 7. *[0 1 2 3];
tspan = (0., 80.);
#subject_dose_amt = p_cyno.BW * ones(length(subject_dose_time)) * TDBc_ugperkg; # [ug]
subject_dose_amt = p_cyno.BW * [1, 0.6, 1e-9, 0.025] * TDBc_ugperkg; # [ug]

u03 = ComponentArray(actTtiss=0., actTpb=0., act0Tpb=0., act0Ttiss=0., act0Ttiss2=0., actTtiss2=0., act0Ttiss3=0., actTtiss3=0., act0Ttumor=0., actTtumor=0., 
IL6pb=0., IL6tiss=0., IL6tiss2=0., IL6tiss3=0., IL6tumor=0., 
injection_effect=injection_effect_init, drug_effect=0., RTXc_ugperkg=0., RTXp_ugperkg=0., Blinc_ug=0., TDBc_ugperml= TDBc_ugperkg * p_cyno.BW / p_cyno_2.V1_mosun, TDBp_ugperml=0. ,
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

cbs = [];
if length(subject_dose_time) > 1
    for i in 2:length(subject_dose_time)
        function affect!(integrator)
            integrator.u.TDBc_ugperml += subject_dose_amt[i]/p_cyno_2.V1_mosun;
            integrator.u.injection_effect += injection_effect_init;
        end
        cb = PresetTimeCallback(subject_dose_time[i],affect!);
        global cbs = push!(cbs, cb);
    end
end
cbset = CallbackSet(cbs...);

sol4 = solve(ODEProblem(Hosseini_ode_homo!, u03, tspan, p_cyno_2), alg=AutoTsit5(Rosenbrock23()), reltol=1e-18, saveat = 0.2, callback = cbset);
sdf4 = DataFrame(sol4);
rename!(sdf4, [:value1, :value2, :value3, :value4, :value5, :value6, :value7, :value8, :value9, :value10, :value11, :value12, :value13, :value14, :value15, 
        :value16, :value17, :value18, :value19, :value20, :value21, :value22, :value23, :value24, :value25, :value26, :value27, :value28, :value29, :value30, 
        :value31, :value32, :value33, :value34,] .=> [:actTtiss, :actTpb, :act0Tpb, :act0Ttiss, :act0Ttiss2, :actTtiss2, :act0Ttiss3, :actTtiss3, :act0Ttumor, 
        :actTtumor, :IL6pb, :IL6tiss, :IL6tiss2, :IL6tiss3, :IL6tumor, :injection_effect, :drug_effect, :RTXc_ugperkg, :RTXp_ugperkg, :Blinc_ug, :TDBc_ugperml, 
        :TDBp_ugperml, :Bpb, :Btiss, :Btiss2, :B1920tiss3, :B19no20tiss3, :Btumor, :restTpb, :restTtiss, :restTtiss2, :restTtiss3, :restTtumor, :BAFF] );

# save the results
CSV.write("old-sims.csv", sdf4);


# plasma PK (supp fig 1)
p1 = plot(title = "mosun, PB", titlefontsize = 8, xticks = (0:7:70, string.(0:7:70)) );
plot!(sdf4.timestamp, sdf4.TDBc_ugperml ,label = "cyno param (supp table 2)", alpha = 0.5, lw = 4);
xlabel!("Time (day)"); ylabel!("T-cell dependent bispecific (ug/mL)"); ylims!(1E-5, 1E2); plot!(yaxis=:log); xlims!(0, 35)

# B/T ratio in spleen (supp fig 3)
p2 = plot(title = "B/T ratio in spleen", titlefontsize = 8, size = (400, 300));
plot!(sdf4.timestamp, sdf4.Btiss ./ (sdf4.actTtiss .+ sdf4.restTtiss) ,label = false, alpha = 0.5, lw = 4);
xlabel!("Time (day)"); ylabel!("spleen B:T ratio");
xlims!(0, 30)

# B/T ratio in LN (supp fig 3)
p3 = plot(title = "B/T ratio in LNs", titlefontsize = 8, size = (400, 300));
plot!(sdf4.timestamp, sdf4.Btiss2 ./ (sdf4.actTtiss2 .+ sdf4.restTtiss2) ,label = false, alpha = 0.5, lw = 4);
xlabel!("Time (day)"); ylabel!("LN B:T ratio");
xlims!(0, 30)

# B cell in BM
p4 = plot(title = "B cells, PB", titlefontsize = 8, xticks = (0:7:70, string.(0:7:70)), size = (400, 300));
plot!(sdf4.timestamp, sdf4.B1920tiss3 / (p_cyno_2.Vtissue3*1E3) ,label = "CD19+CD20+ cells in BM", alpha = 0.5, lw = 4);
plot!(sdf4.timestamp, sdf4.B19no20tiss3 / (p_cyno_2.Vtissue3*1E3) ,label = "CD19+CD20- cells in BM", alpha = 0.5, lw = 4);
xlabel!("Time (day)"); ylabel!("B-cell (#/uL)"); ylims!(1E3, 1E5); plot!(yaxis = :log)

# injection effect
p5 = plot(title = "injection effect", titlefontsize = 8, xticks = (0:7:70, string.(0:7:70)), size = (400, 300));
plot!(sdf4.timestamp, sdf4.injection_effect,label = false, alpha = 0.5, lw = 4);
xlabel!("Time (day)"); ylabel!("injection effect")

# CD8+ T cells in PB (fig 3)
p6 = plot(title = "CD8+ T cells, PB", titlefontsize = 8, xticks = (0:7:70, string.(0:7:70)), size = (400, 300));
plot!(sdf4.timestamp, (sdf4.actTpb .+  sdf4.act0Tpb .+ sdf4.restTpb) / (p_cyno_2.Vpb*1E3) ,label = false, alpha = 0.5, lw = 4);
xlabel!("Time (day)"); ylabel!("CD8+ T-cell (#/uL)"); ylims!(0, 8000); xlims!(0, 30)

# CD69+ CD8+ T cell % in PB (fig 3)
p7 = plot(title = "CD8+ T cells, PB", titlefontsize = 8, xticks = (0:7:70, string.(0:7:70)), size = (400, 300));
plot!(sdf4.timestamp, sdf4.actTpb ./ (sdf4.actTpb .+  sdf4.act0Tpb .+ sdf4.restTpb) * 100 ,label = false, alpha = 0.5, lw = 4);
xlabel!("Time (day)"); ylabel!("CD69+ CD8+ T-cell (%)"); ylims!(-1, 100); xlims!(0, 30)

# IL6 concentration (supp fig 4)
p8 = plot(title = "IL6, PB", titlefontsize = 8, xticks = (0:7:70, string.(0:7:70)), size = (400, 300));
plot!(sdf4.timestamp, sdf4.IL6pb,label = false, alpha = 0.5, lw = 4);
xlabel!("Time (day)"); ylabel!("plasma IL6 (pg/mL)");
ylims!(0.1, 1e4); xlims!(0, 30); plot!(yaxis=:log)

# B cell in PB (fig 3)
p9 = plot(title = "B cells, PB", titlefontsize = 8, xticks = (0:7:70, string.(0:7:70)), size = (400, 300));
plot!(sdf4.timestamp, sdf4.Bpb / (p_cyno_2.Vpb*1E3) ,label = false, alpha = 0.5, lw = 4);
xlabel!("Time (day)"); ylabel!("B-cell (#/uL)"); ylims!(1, 2000); xlims!(0, 30)  #; plot!(yaxis=:log)

# save figures 
savefig(p1, "cyno.pdf")
savefig(p2, "temp.pdf"); append_pdf!("cyno.pdf", "temp.pdf", cleanup=true);
savefig(p3, "temp.pdf"); append_pdf!("cyno.pdf", "temp.pdf", cleanup=true);
savefig(p4, "temp.pdf"); append_pdf!("cyno.pdf", "temp.pdf", cleanup=true);
savefig(p5, "temp.pdf"); append_pdf!("cyno.pdf", "temp.pdf", cleanup=true);
savefig(p6, "temp.pdf"); append_pdf!("cyno.pdf", "temp.pdf", cleanup=true);
savefig(p7, "temp.pdf"); append_pdf!("cyno.pdf", "temp.pdf", cleanup=true);
savefig(p8, "temp.pdf"); append_pdf!("cyno.pdf", "temp.pdf", cleanup=true);
savefig(p9, "temp.pdf"); append_pdf!("cyno.pdf", "temp.pdf", cleanup=true);

