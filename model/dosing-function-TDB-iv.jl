# date: 8/30/24

# default dose following mosunetuzumab label dosing scheme
function tdb_dosing_iv(subject_dose_amt_TDB = [1., 2., 60., 60., 30.] * 1E3, subject_dose_time_TDB = 7. *[0 1 2 3 6], injection_effect_init = 10., p_homo_0 = p_homo)
    global cbs01 = [];
    if subject_dose_time_TDB[1] == 0
        for i in 2:length(subject_dose_time_TDB)
            function affect_TDB!(integrator)
                integrator.u.TDBc_ugperml += subject_dose_amt_TDB[i] / p_homo_0.V1_TDB ;
                integrator.u.injection_effect += injection_effect_init;
            end
            cb01 = PresetTimeCallback(subject_dose_time_TDB[i],affect_TDB!);
            global cbs01 = push!(cbs01, cb01);
        end
    else
        for i in 1:length(subject_dose_time_TDB)
            function affect_TDB!(integrator)
                integrator.u.TDBc_ugperml += subject_dose_amt_TDB[i] / p_homo_0.V1_TDB ;
                integrator.u.injection_effect += injection_effect_init;
            end
            cb01 = PresetTimeCallback(subject_dose_time_TDB[i],affect_TDB!);
            global cbs01 = push!(cbs01, cb01);
        end
    end
    cbset01 = CallbackSet(cbs01...);
    return cbset01; 
end
