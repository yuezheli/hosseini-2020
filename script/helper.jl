# author: Yuezhe LI
# date: 10/3/2023
# purpose of this code: to host helper functions for simulation using bs-lonca-3.jl

include("../model/param.jl");
include("../model/init.jl");
include("../model/bs-lonca-3.jl");

function preequlibrium3(tumor_vol,  # [mL]
    cd19_neg_frac, p_homo, injection_effect_init = 10., u0_eql = u0_0)
    tspan = (0., 7.);
    p_01 = deepcopy(p_homo); 
    p_01.kBtumorprolif = 0.
    u0_tmp = deepcopy(u0_eql);
    u0_tmp.Btumor = tumor_vol * 0.375/ Vc * (1 - cd19_neg_frac); 
    u0_tmp.Btumor_cd19neg = tumor_vol * 0.375/ Vc * cd19_neg_frac; 
    # solve ODE 
    sol01 = solve(ODEProblem(bs_lonca_ode_3!, u0_tmp, tspan, p_01), alg=AutoTsit5(Rosenbrock23()), saveat = 1.);
    u_end = sol01.u[end]
    for i = 1:length(u_end)
        u_end[i] = max(u_end[i], 0.0);
    end
    u_end.injection_effect = injection_effect_init
    return u_end;
end

