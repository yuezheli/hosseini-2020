# date: 8/30/24
# author: Yuezhe Li 

# default dose following mosunetuzumab dosing scheme in LOTIS-7 trial
function tdb_dosing_sc(subject_dose_amt_TDB = [5., 45., 45., 45., 45.] * 1E3, subject_dose_time_TDB = 7. *[0 1 2 3 6], injection_effect_init = 10., p_homo_0 = p_homo)
    global cbs01 = [];
    if subject_dose_time_TDB[1] == 0
        for i in 2:length(subject_dose_time_TDB)
            function affect_TDB!(integrator)
                integrator.u.TDBdepot_ug += subject_dose_amt_TDB[i] * p_homo_0.fbio_TDB;
                integrator.u.injection_effect += injection_effect_init;
            end
            cb01 = PresetTimeCallback(subject_dose_time_TDB[i],affect_TDB!);
            global cbs01 = push!(cbs01, cb01);
        end
    else 
        for i in 1:length(subject_dose_time_TDB)
            function affect_TDB!(integrator)
                integrator.u.TDBdepot_ug += subject_dose_amt_TDB[i] * p_homo_0.fbio_TDB;
                integrator.u.injection_effect += injection_effect_init;
            end
            cb01 = PresetTimeCallback(subject_dose_time_TDB[i],affect_TDB!);
            global cbs01 = push!(cbs01, cb01);
        end
    end
    cbset01 = CallbackSet(cbs01...);
    return cbset01; 
end
