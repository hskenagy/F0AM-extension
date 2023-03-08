function [cond_table, time_series] = calculate_RO2_parameters(S, timemin, VOC_name, int_name)

% This is a function that calculates a variety of RO2-chemistry-related 
% parameters (e.g., beta, tau, etc.).  You shouldn’t have to modify 
% anything in this file unless you want to change how parameters are 
% calculated or introduce new parameters. 


% INPUTS
% S: F0AM output for a single simulation
% timemin: simulation time (in seconds) to begin integration at
% VOC_name: name of VOC precursor used
% int_name: name of intermediate product to use when defining generations
% (recommend MVK for C5H8, APINBOOH for APINENE)

% OUTPUTS
% cond_table: table containing all of the important RO2-related parameters
% for each experiment
% time_series: structure containing a couple important RO2-related time series

cond_table = table;

time = timemin:1:max(S.Time);
time_series.time = time;

% calculate cutoff times between each generation
% use intermediate compound to distinguish gen1 and gen2 (e.g. MVK for
% isoprene)
% assign gen1 as everything before intermediate compound peaks
% assign gen2 ending when intermediate compound reaches half its maximum
interp_int = interp1(S.Time(S.Time>timemin), ...
    S.Conc.(int_name)(S.Time>timemin, :), time);
[~, int_max] = max(interp_int(:));
k = find(interp_int < 0.5 * max(interp_int));
half_int_ind = min(k(k>int_max));
if isempty(half_int_ind)
    half_int_ind = length(interp_int);
end
cond_table.half_int_ind = half_int_ind;
time_series.interp_int = interp_int;

% calculate cutoff times for single generation situations (e.g. a-pinene)
interp_VOC = interp1(S.Time(S.Time>timemin), ...
    S.Conc.(VOC_name)(S.Time>timemin, :), time);
l = find(interp_VOC < 0.1 * max(interp_VOC));
x90_VOC_ind = min(l);
if isempty(x90_VOC_ind)
    x90_VOC_ind = length(interp_VOC);
end
cond_table.x90_VOC_ind = x90_VOC_ind;
time_series.interp_VOC = interp_VOC;

        
rates = sumRates_ROx(S);
nox_rates = sumRates_NOx(S);
PANs_prod = nox_rates.Loss(:,3);
losses = [rates.Loss(:,[3,4,5,6]), PANs_prod];

% interpolate values so they are on a linear timescale (otherwise
% medians and means don't make sense)
interp_losses = interp1(S.Time(S.Time>timemin),...
    losses(S.Time>timemin, :), time);
cumul_losses = cumsum(interp_losses, 'omitnan');
gen1_cumul_losses = cumsum(interp_losses(1:int_max,:), 1, 'omitnan');
gen2_cumul_losses = cumsum(interp_losses(int_max:half_int_ind,:), 1, 'omitnan');

RO2 = struct2cell(ExtractSpecies(S.iRO2,S));
RO2 = sum([RO2{:}],2);
interp_RO2 = interp1(S.Time(S.Time>timemin),...
    RO2(S.Time>timemin, :), time);
time_series.interp_RO2 = interp_RO2;

% beta time series
time_series.beta_t = interp_losses(:,3) ./ (interp_losses(:,3) + interp_losses(:,1));
% beta integrated over the entire 8-hr experiment
cond_table.beta_8hr = cumul_losses(end,3) ./ (cumul_losses(end,3) + cumul_losses(end,1));
% beta integrated until intermediate species reaches half its max
cond_table.beta_half_int = cumul_losses(half_int_ind,3) ./ (cumul_losses(half_int_ind, 3) + cumul_losses(half_int_ind,1));
% beta integrated until 90% of the initial VOC has been oxidized
cond_table.beta_90 = cumul_losses(x90_VOC_ind,3) ./ (cumul_losses(x90_VOC_ind, 3) + cumul_losses(x90_VOC_ind,1));
% betas integrated for each generation (defined based on intermeidate
% concentrations, e.g. MVK)
cond_table.beta_gen1 = gen1_cumul_losses(end,3) ./ ...
    (gen1_cumul_losses(end,3) + gen1_cumul_losses(end,1));
cond_table.beta_gen2 = gen2_cumul_losses(end,3) ./ ...
    (gen2_cumul_losses(end,3) + gen2_cumul_losses(end,1));

% gamma time series
time_series.gamma_t = interp_losses(:,3) ./ (interp_losses(:,3) + interp_losses(:,5));
% gamma integrated over the entire 8-hr experiment
cond_table.gamma_8hr = cumul_losses(end,3) ./ (cumul_losses(end,3) + cumul_losses(end,5));
% beta integrated until intermediate species reaches half its max
cond_table.gamma_half_int = cumul_losses(half_int_ind,3) ./ (cumul_losses(half_int_ind, 3) + cumul_losses(half_int_ind,5));
% gamma integrated until 90% of the initial VOC has been oxidized
cond_table.gamma_90 = cumul_losses(x90_VOC_ind,3) ./ (cumul_losses(x90_VOC_ind, 3) + cumul_losses(x90_VOC_ind,5));
% gammas integrated for each generation (defined based on intermeidate
% concentrations, e.g. MVK)
cond_table.gamma_gen1 = gen1_cumul_losses(end,3) ./ ...
    (gen1_cumul_losses(end,3) + gen1_cumul_losses(end,5));
cond_table.gamma_gen2 = gen2_cumul_losses(end,3) ./ ...
    (gen2_cumul_losses(end,3) + gen2_cumul_losses(end,5));


% tau
time_series.tau_t = (-1 * interp_RO2' ./ ...
    (interp_losses(:,1) + interp_losses(:,3)));
cond_table.tau_8hr = median(-1 * interp_RO2' ./ ...
    (interp_losses(:,1) + interp_losses(:,3)), 'omitnan'); % tau_RO2
cond_table.tau_half_int = median(-1 * interp_RO2(:,1:half_int_ind)' ./ ...
    (interp_losses(1:half_int_ind,1) + interp_losses(1:half_int_ind,3)), 'omitnan');
cond_table.tau_90 = median(-1 * interp_RO2(:,1:x90_VOC_ind)' ./ ...
    (interp_losses(1:x90_VOC_ind,1) + interp_losses(1:x90_VOC_ind,3)), 'omitnan');% tau_RO2
cond_table.tau_gen1 = median(-1 * interp_RO2(:,1:int_max)' ./ ...
    (interp_losses(1:int_max,1) + interp_losses(1:int_max,3)), 'omitnan');
cond_table.tau_gen2 = median(-1 * interp_RO2(:,int_max:half_int_ind)' ./ ...
    (interp_losses(int_max:half_int_ind,1) + interp_losses(int_max:half_int_ind,3)), 'omitnan');

% tau acyl
time_series.tau_acyl_t = (-1 * interp_RO2' ./ ...
    (interp_losses(:,1) + interp_losses(:,3) + interp_losses(:,4)));
cond_table.tau_acyl_8hr = median(-1 * interp_RO2' ./ ...
    (interp_losses(:,1) + interp_losses(:,3) + interp_losses(:,4)), 'omitnan'); % tau_RO2_acyl
cond_table.tau_acyl_half_int = median(-1 * interp_RO2(:,1:half_int_ind)' ./ ...
    (interp_losses(1:half_int_ind,1) + interp_losses(1:half_int_ind,3) + interp_losses(1:half_int_ind,4)), 'omitnan');
cond_table.tau_acyl_90 = median(-1 * interp_RO2(:,1:x90_VOC_ind)' ./ ...
    (interp_losses(1:x90_VOC_ind,1) + interp_losses(1:x90_VOC_ind,3) + interp_losses(1:x90_VOC_ind,4)), 'omitnan');
cond_table.tau_acyl_gen1 = median(-1 * interp_RO2(:,1:int_max)' ./ ...
    (interp_losses(1:int_max,1) + interp_losses(1:int_max,3) + interp_losses(1:int_max,4)), 'omitnan');
cond_table.tau_acyl_gen2 = median(-1 * interp_RO2(:,int_max:half_int_ind)' ./ ...
    (interp_losses(int_max:half_int_ind,1) + interp_losses(int_max:half_int_ind,3) + interp_losses(int_max:half_int_ind,4)), 'omitnan');

OH_interp = interp1(S.Time(S.Time>timemin), ...
    S.Conc.OH(S.Time>timemin), time);
cond_table.OH = median(OH_interp, 'omitnan');
HO2_interp = interp1(S.Time(S.Time>timemin), ...
    S.Conc.HO2(S.Time>timemin), time);
cond_table.HO2 = median(HO2_interp, 'omitnan');
NO_interp = interp1(S.Time(S.Time>timemin), ...
    S.Conc.NO(S.Time>timemin), time);
NO2_interp = interp1(S.Time(S.Time>timemin), ...
    S.Conc.NO2(S.Time>timemin), time);
cond_table.NO_NO2 = median((NO_interp./NO2_interp), 'omitnan');
cond_table.max_VOC = max(S.Conc.(VOC_name)(:));
cond_table.end_VOC = S.Conc.(VOC_name)(end);
cond_table.max_int = max(S.Conc.(int_name)(:));
cond_table.end_int = S.Conc.(int_name)(end);

Sprates = PlotRates(VOC_name,S,3,'plotme',0);
interp_VOC_losses = interp1(S.Time(S.Time>timemin),...
    Sprates.Loss(S.Time>timemin, :), time);
if strcmp(VOC_name, 'APINENE')
    z = 11;
elseif strcmp(VOC_name, "C5H8")
    z = 1;
else
    disp('unsupported VOC')
end
OH_ind = cellfun(@(c) isequal(c,z), regexp(Sprates.Lnames, 'OH'));
OH_loss = interp_VOC_losses(:,find(OH_ind==1));
O3_ind = cellfun(@(c) isequal(c,z), regexp(Sprates.Lnames, 'O3'));
O3_loss = interp_VOC_losses(:,find(O3_ind==1));
NO3_ind = cellfun(@(c) isequal(c,z), regexp(Sprates.Lnames, 'NO3'));
NO3_loss = interp_VOC_losses(:,find(NO3_ind==1));
sum_loss = OH_loss + O3_loss + NO3_loss;
total_loss = sum(sum_loss, 'omitnan');
cond_table.frac_OH = sum(OH_loss, 'omitnan') / total_loss;
cond_table.frac_O3 = sum(O3_loss, 'omitnan') / total_loss;
cond_table.frac_NO3 = sum(NO3_loss, 'omitnan') / total_loss;
time_series.frac_OH_t = OH_loss ./ sum_loss;
time_series.frac_O3_t = O3_loss ./ sum_loss;
time_series.frac_NO3_t = NO3_loss ./ sum_loss;

if isfield(S.Conc, 'HONO')
    cond_table.inital_HONO = S.Conc.HONO(1);
else
    cond_table.initial_HONO = 0;
end
if isfield(S.Conc, 'H2O2')
    cond_table.inital_H2O2 = S.Conc.H2O2(1);
else
    cond_table.initial_H2O2 = 0;
end
if isfield(S.Conc, 'CH3ONO')
    cond_table.inital_CH3ONO = S.Conc.CH3ONO(1);
else
    cond_table.initial_CH3ONO = 0;
end
if isfield(S.Conc, 'NO')
    cond_table.inital_NO = S.Conc.NO(1);
else
    cond_table.initial_NO = 0;
end
if isfield(S.Conc, 'NO2')
    cond_table.inital_NO2 = S.Conc.NO2(1);
else
    cond_table.initial_NO2 = 0;
end
if isfield(S.Conc, VOC_name)
    cond_table.inital_VOC = S.Conc.(VOC_name)(1);
else
    cond_table.initial_VOC = 0;
end

end