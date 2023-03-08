# F0AM-extension
Extension for F0AM that allows for iterating through combinations of initial conditions and calculating a variety of parameters related to RO2 chemistry for each simulation.

There are two relevant files:

1. F0AM_span_initial_conditions
This is the script that generates the initial condition matrix, runs a F0AM simulation for each of the initial conditions in that matrix, and then tabulates the calculated RO2-chemistry-related parameters.  This is the file where you specify which initial conditions you want to vary, and over what ranges.  It’s also where you specify the usual met variables and chem files to use.  And, at the end of the file, you can specify what to call the tabulated output for saving.

2. calculate_RO2_parameters
This is a function that calculates a variety of RO2-chemistry-related parameters (e.g., beta, tau, etc.).  You shouldn’t have to modify anything in this file unless you want to change how parameters are calculated or introduce new parameters. 

There are two important outputs that will be saved when you run F0AM_span_initial conditions:

1. cond_table: table with a row for each experiment tabulating all of the calculated parameters for each experiment

2. time_series: structure with a field for each experiment, each of which contains a couple of RO2-related time series

The parameters tabulated in cond_table are as follows:

half_int_ind: time index at which the intermediate species selected reaches half of its maximum

x90_VOC_ind: time index at which 90% of the initial VOC has been consumed

beta_8hr: beta for the entire 8-hr experiment

beta_half_int: beta until the intermediate species selected reaches half of its maximum

beta_90: beta until 90% of the initial VOC has been consumed

beta_gen1: beta for the first generation (until the intermediate species selected hits its maximum)

beta_gen2: beta for the second generation (starting from when the intermediate species selected hits its maximum, ending when the intermediate species selected reaches half of its maximum)

gamma_8hr: gamma for the entire 8-hr experiment

gamma_half_int: gamma until the intermediate species selected reaches half of its maximum

gamma_90: gamma until 90% of the initial VOC has been consumed

gamma_gen1: gamma for the first generation (until the intermediate species selected hits its maximum)

gamma_gen2: gamma for the second generation (starting from when the intermediate species selected hits its maximum, ending when the intermediate species selected reaches half of its maximum)

tau_8hr: tau for the entire 8-hr experiment

tau_half_int: tau until the intermediate species selected reaches half of its maximum

tau_90: tau until 90% of the initial VOC has been consumed

tau_gen1: tau for the first generation (until the intermediate species selected hits its maximum)

tau_gen2: tau for the second generation (starting from when the intermediate species selected hits its maximum, ending when the intermediate species selected reaches half of its maximum)

tau_acyl_8hr: tau for acyl RO2 for the entire 8-hr experiment

tau_acyl_half_int: tau for acyl RO2 until the intermediate species selected reaches half of its maximum

tau_acyl_90: tau for acyl RO2 until 90% of the initial VOC has been consumed

tau_acyl_gen1: tau for acyl RO2 for the first generation (until the intermediate species selected hits its maximum)

tau_acyl_gen2: tau for acyl RO2 for the second generation (starting from when the intermediate species selected hits its maximum, ending when the intermediate species selected reaches half of its maximum)

OH: median OH over the experiment

HO2: median HO2 over the experiment

NO_NO2: median NO/NO2 ratio over the experiment

max_VOC: maximum VOC concentration

end_VOC: final VOC concentration

max_int: maximum intermediate species concentration

end_int: final intermediate species concentration

frac_OH: fraction of VOC that reacts with OH

frac_O3: fraction of VOC that reacts with O3

frac_NO3: fraction of VOC that reacts with NO3

initial_HONO: initial HONO concentration

initial_H2O2: initial H2O2 concentration

initial_CH3ONO: initial CH3ONO concentration

initial_NO: initial NO concentration

initial_NO2: initial NO2 concentration

initial_VOC: initial VOC concentration

The parameters saved in time_series are as follow:

time: on a one-second time base

interp_int: intermediate species concentration interpolated on one-second time base

interp_VOC: VOC concentration interpolated on one-second time base

interp_RO2: total RO2 concentration interpolated on one-second time base

beta_t: instantaneous beta interpolated on one-second time base

gamma_t: instantaneous gamma interpolated on one-second time base

tau_t: instantaneous tau interpolated on one-second time base

tau_acyl_t: instantaneous tau for acyl RO2 interpolated on one-second time base

frac_OH_t: instantaneous fraction of VOC oxidized by OH interpolated on one-second time base

frac_O3_t: instantaneous fraction of VOC oxidized by O3 interpolated on one-second time base

frac_NO3_t: instantaneous fraction of VOC oxidized by NO3 interpolated on one-second time base

General definitions:

Beta definition: fraction of RO2 that reacts with NO vs HO2 (beta = 1 means all reaction with NO, beta = 0 means all reaction with HO2)

Tau definition: RO2 bimolecular lifetime (taking into account reaction with NO and HO2)

Gamma definition: fraction of RO2 that reacts with NO vs NO2 (gamma = 1 means all reaction with NO, gamma = 0 means all reaction with NO2)
