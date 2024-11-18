# date: 8/30/24
# author: Yuezhe LI 
# purpose of this code: to simulate IV dosing events for ADCs

function adc_iv_dosing(dose_amt, dose_time, V1_ADC)
    global cbs = [];
    if length(dose_time) > 1
        for i in 1:length(dose_time)
            function affect!(integrator)
                integrator[:ADCc_ugperml] += dose_amt[i]/ V1_ADC;
            end
            cb = PresetTimeCallback(dose_time[i],affect!);
            global cbs = push!(cbs, cb);
        end
    end
    cbset = CallbackSet(cbs...);
    return cbset
end

