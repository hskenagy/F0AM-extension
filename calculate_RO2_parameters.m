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
time_series.VOC = interp_VOC;

        
% rates = sumRates_ROx(S);
% nox_rates = sumRates_NOx(S);
% PANs_prod = nox_rates.Loss(:,3);
% losses = [rates.Loss(:,[3,4,5,6]), PANs_prod];
losses = sumRates_RO2loss(S).Loss;
% 1 - RO2 + HO2
% 2 - RO2 + RO2
% 3 - RO2 + NO
% 4 - RO2 + NO2
% 5 - acylRO2 + NO2
% 6 - acylRO2 + NO
% 7 - RO2 + NO3
% 8 - other RO2 loss reactions

% interpolate values so they are on a linear timescale (otherwise
% means and medians don't make sense)
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
k_NO_RO2 = 2.7e-12 * exp(360./S.Met.T);
k_HO2_RO2 = 2.91e-13 * exp(1300./S.Met.T);
time_series.beta_nosum = (k_NO_RO2 .* S.Conc.NO) ./ ...
    (k_NO_RO2 .* S.Conc.NO + k_HO2_RO2 .* S.Conc.HO2);


% tau
time_series.tau_t = (-1 * interp_RO2' ./ ...
    (interp_losses(:,1) + interp_losses(:,3)));
cond_table.tau_8hr = mean(-1 * interp_RO2' ./ ...
    (interp_losses(:,1) + interp_losses(:,3)), 'omitnan'); % tau_RO2
cond_table.tau_half_int = mean(-1 * interp_RO2(:,1:half_int_ind)' ./ ...
    (interp_losses(1:half_int_ind,1) + interp_losses(1:half_int_ind,3)), 'omitnan');
cond_table.tau_90 = mean(-1 * interp_RO2(:,1:x90_VOC_ind)' ./ ...
    (interp_losses(1:x90_VOC_ind,1) + interp_losses(1:x90_VOC_ind,3)), 'omitnan');% tau_RO2
cond_table.tau_gen1 = mean(-1 * interp_RO2(:,1:int_max)' ./ ...
    (interp_losses(1:int_max,1) + interp_losses(1:int_max,3)), 'omitnan');
cond_table.tau_gen2 = mean(-1 * interp_RO2(:,int_max:half_int_ind)' ./ ...
    (interp_losses(int_max:half_int_ind,1) + interp_losses(int_max:half_int_ind,3)), 'omitnan');
time_series.tau_nosum = 1 ./ (k_NO_RO2 * S.Conc.NO * NumberDensity(S.Met.P, S.Met.T) * 1e-9 + ...
    k_HO2_RO2 * S.Conc.HO2 * NumberDensity(S.Met.P, S.Met.T) * 1e-9);
cond_table.tau_90_withRO2 = mean(-1 * interp_RO2(:,1:x90_VOC_ind)' ./ ...
    (interp_losses(1:x90_VOC_ind,1) + ...
    interp_losses(1:x90_VOC_ind,3) + ...
    interp_losses(1:x90_VOC_ind,2)), 'omitnan');% tau_RO2


% gamma time series
time_series.gamma_t = interp_losses(:,6) ./ (interp_losses(:,6) + interp_losses(:,5));
% gamma integrated over the entire 8-hr experiment
cond_table.gamma_8hr = cumul_losses(end,6) ./ (cumul_losses(end,6) + cumul_losses(end,5));
% beta integrated until intermediate species reaches half its max
cond_table.gamma_half_int = cumul_losses(half_int_ind,6) ./ (cumul_losses(half_int_ind, 6) + cumul_losses(half_int_ind,5));
% gamma integrated until 90% of the initial VOC has been oxidized
cond_table.gamma_90 = cumul_losses(x90_VOC_ind,6) ./ (cumul_losses(x90_VOC_ind, 6) + cumul_losses(x90_VOC_ind,5));
% gammas integrated for each generation (defined based on intermeidate
% concentrations, e.g. MVK)
cond_table.gamma_gen1 = gen1_cumul_losses(end,6) ./ ...
    (gen1_cumul_losses(end,6) + gen1_cumul_losses(end,5));
cond_table.gamma_gen2 = gen2_cumul_losses(end,6) ./ ...
    (gen2_cumul_losses(end,6) + gen2_cumul_losses(end,5));


% rho time series
time_series.rho_t = interp_losses(:,2) ./ (interp_losses(:,2)+ interp_losses(:,1));
% rho integrated over the entire 8-hr experiment
cond_table.rho_8hr = cumul_losses(end,2) ./ (cumul_losses(end,2) + cumul_losses(end,1));
% rho integrated until intermediate species reaches half its max
cond_table.rho_half_int = cumul_losses(half_int_ind,2) ./ (cumul_losses(half_int_ind,2) + cumul_losses(half_int_ind,1));
% rho integrated until 90% of the initial VOC has been oxidized
cond_table.rho_90 = cumul_losses(x90_VOC_ind,2) ./ (cumul_losses(x90_VOC_ind,2) + cumul_losses(x90_VOC_ind,1));
% rhos integrated for each generation (defined based on intermeidate
% concentrations, e.g. MVK)
cond_table.rho_gen1 = gen1_cumul_losses(end,2) ./ ...
    (gen1_cumul_losses(end,2) + gen1_cumul_losses(end,1));
cond_table.rho_gen2 = gen2_cumul_losses(end,2) ./ ...
    (gen2_cumul_losses(end,2) + gen2_cumul_losses(end,1));

OH_interp = interp1(S.Time(S.Time>timemin), ...
    S.Conc.OH(S.Time>timemin), time);
cond_table.OH = mean(OH_interp, 'omitnan');
time_series.OH = OH_interp;
HO2_interp = interp1(S.Time(S.Time>timemin), ...
    S.Conc.HO2(S.Time>timemin), time);
cond_table.HO2 = mean(HO2_interp, 'omitnan');
NO_interp = interp1(S.Time(S.Time>timemin), ...
    S.Conc.NO(S.Time>timemin), time);
cond_table.NO = mean(NO_interp, 'omitnan');
NO2_interp = interp1(S.Time(S.Time>timemin), ...
    S.Conc.NO2(S.Time>timemin), time);
cond_table.NO2 = mean(NO2_interp, 'omitnan');
RO2_interp = interp1(S.Time(S.Time>timemin), ...
    S.Conc.RO2(S.Time>timemin), time);
cond_table.RO2 = mean(RO2_interp, 'omitnan');
O3_interp = interp1(S.Time(S.Time>timemin), ...
    S.Conc.O3(S.Time>timemin), time);
cond_table.O3 = mean(O3_interp, 'omitnan');
NO3_interp = interp1(S.Time(S.Time>timemin), ...
    S.Conc.NO3(S.Time>timemin), time);
cond_table.NO3 = mean(NO3_interp, 'omitnan');
CH3O2_interp = interp1(S.Time(S.Time>timemin), ...
    S.Conc.CH3O2(S.Time>timemin), time);
cond_table.CH3O2 = mean(CH3O2_interp, 'omitnan');
cond_table.max_VOC = max(S.Conc.(VOC_name)(:));
cond_table.end_VOC = S.Conc.(VOC_name)(end);
cond_table.max_int = max(S.Conc.(int_name)(:));
cond_table.end_int = S.Conc.(int_name)(end);

time_series.NO_NO2_t = NO_interp ./ NO2_interp;
cond_table.NO_NO2_8hr = mean((NO_interp./NO2_interp), 'omitnan');
cond_table.NO_NO2_half_int = mean((NO_interp(1:half_int_ind)./...
    NO2_interp(1:half_int_ind)), 'omitnan');
cond_table.NO_NO2_half_int = mean((NO_interp(1:x90_VOC_ind)./...
    NO2_interp(1:x90_VOC_ind)), 'omitnan');
cond_table.NO_NO2_gen1 = mean((NO_interp(1:int_max)./...
    NO2_interp(1:int_max)), 'omitnan');
cond_table.NO_NO2_gen2 = mean((NO_interp(int_max:half_int_ind)./...
    NO2_interp(int_max:half_int_ind)), 'omitnan');
cond_table.NO_NO2_90 = mean((NO_interp(1:x90_VOC_ind)./...
    NO2_interp(1:x90_VOC_ind)), 'omitnan');

time_series.RO2_HO2_t = RO2_interp ./ HO2_interp;
cond_table.RO2_HO2_8hr = mean((RO2_interp./HO2_interp), 'omitnan');
cond_table.RO2_HO2_half_int = mean((RO2_interp(1:half_int_ind)./...
    HO2_interp(1:half_int_ind)), 'omitnan');
cond_table.RO2_HO2_half_int = mean((RO2_interp(1:x90_VOC_ind)./...
    HO2_interp(1:x90_VOC_ind)), 'omitnan');
cond_table.RO2_HO2_gen1 = mean((RO2_interp(1:int_max)./...
    HO2_interp(1:int_max)), 'omitnan');
cond_table.RO2_HO2_gen2 = mean((RO2_interp(int_max:half_int_ind)./...
    HO2_interp(int_max:half_int_ind)), 'omitnan');
cond_table.RO2_HO2_90 = mean((RO2_interp(1:x90_VOC_ind)./...
    HO2_interp(1:x90_VOC_ind)), 'omitnan');

time_series.CH3O2_RO2_t = CH3O2_interp ./ RO2_interp;
cond_table.CH3O2_RO2_8hr = mean((CH3O2_interp./RO2_interp), 'omitnan');
cond_table.CH3O2_RO2_half_int = mean((CH3O2_interp(1:half_int_ind)./...
    RO2_interp(1:half_int_ind)), 'omitnan');
cond_table.CH3O2_RO2_half_int = mean((CH3O2_interp(1:x90_VOC_ind)./...
    RO2_interp(1:x90_VOC_ind)), 'omitnan');
cond_table.CH3O2_RO2_gen1 = mean((CH3O2_interp(1:int_max)./...
    RO2_interp(1:int_max)), 'omitnan');
cond_table.CH3O2_RO2_gen2 = mean((CH3O2_interp(int_max:half_int_ind)./...
    RO2_interp(int_max:half_int_ind)), 'omitnan');


Sprates = PlotRates(VOC_name,S,3,'plotme',0);
if ~isempty(Sprates.Lnames)
interp_VOC_losses = interp1(S.Time(S.Time>timemin),...
    Sprates.Loss(S.Time>timemin, :), time);


[~,iOH] = ismember('OH',S.Cnames);
[~,iO3] = ismember('O3',S.Cnames);
[~,iNO3] = ismember('NO3',S.Cnames);
iG = S.Chem.iG(Sprates.iRx_Loss);
jOH = any(iG == iOH, 2);
jO3 = any(iG == iO3, 2);
jNO3 = any(iG == iNO3, 2);
OH_loss = interp_VOC_losses(:,jOH==1);
O3_loss = interp_VOC_losses(:,jO3==1);
NO3_loss = interp_VOC_losses(:,jNO3==1);

if ~isempty(O3_loss)
    if ~isempty(NO3_loss)
        sum_loss = OH_loss + O3_loss + NO3_loss;
    else
        sum_loss = OH_loss + O3_loss;
    end
else
    sum_loss = OH_loss;
end
total_loss = sum(sum_loss, 'omitnan');
cond_table.frac_OH = sum(OH_loss, 'omitnan') / total_loss;
cond_table.frac_O3 = sum(O3_loss, 'omitnan') / total_loss;
cond_table.frac_NO3 = sum(NO3_loss, 'omitnan') / total_loss;
if isempty(cond_table.frac_NO3)
    cond_table.frac_NO3 = NaN;
end
time_series.frac_OH_t = OH_loss ./ sum_loss;
time_series.frac_O3_t = O3_loss ./ sum_loss;
time_series.frac_NO3_t = NO3_loss ./ sum_loss;

else
    cond_table.frac_OH = NaN;
    cond_table.frac_O3 = NaN;
    cond_table.frac_NO3 = NaN;
end

if isfield(S.Conc, 'HONO')
    cond_table.initial_HONO = S.Conc.HONO(1);
else
    cond_table.initial_HONO = 0;
end
if isfield(S.Conc, 'H2O2')
    cond_table.initial_H2O2 = S.Conc.H2O2(1);
else
    cond_table.initial_H2O2 = 0;
end
if isfield(S.Conc, 'CH3ONO')
    cond_table.initial_CH3ONO = S.Conc.CH3ONO(1);
else
    cond_table.initial_CH3ONO = 0;
end
if isfield(S.Conc, 'NO')
    cond_table.initial_NO = S.Conc.NO(1);
else
    cond_table.initial_NO = 0;
end
if isfield(S.Conc, 'NO2')
    cond_table.initial_NO2 = S.Conc.NO2(1);
else
    cond_table.initial_NO2 = 0;
end
if isfield(S.Conc, VOC_name)
    cond_table.initial_VOC = S.Conc.(VOC_name)(1);
else
    cond_table.initial_VOC = 0;
end

if strcmp(VOC_name, 'C5H8')

delta_VOC = cond_table.initial_VOC - interp_VOC;
% m/z = 118.1311
mz118 = SearchMW(118.1311,S.Cnames,'v331');
isopooh = 0;
for i = 1:length(mz118.index)
    isopooh = isopooh + interp1(S.Time(S.Time>timemin), ...
        S.Conc.(mz118.names{i})(S.Time>timemin, :), time);
end
cond_table.isopooh_yield = mean(isopooh(:,2:x90_VOC_ind) ./...
    delta_VOC(:,2:x90_VOC_ind));
time_series.isopooh = isopooh;

mz148 = SearchMW(147.1293,S.Cnames,'v331');
ihn = 0;
for i = 1:length(mz148.index)
    ihn = ihn + interp1(S.Time(S.Time>timemin), ...
        S.Conc.(mz148.names{i})(S.Time>timemin, :), time);
end
cond_table.ihn_yield = mean(ihn(:,2:x90_VOC_ind) ./...
    delta_VOC(:,2:x90_VOC_ind));
time_series.ihn = ihn;

mz116 = SearchMW(116.1152,S.Cnames,'v331');
hpald = 0;
for i = 1:length(mz116.index)
    hpald = hpald + interp1(S.Time(S.Time>timemin), ...
        S.Conc.(mz116.names{i})(S.Time>timemin, :), time);
end
cond_table.hpald_yield = mean(hpald(:,2:x90_VOC_ind) ./...
    delta_VOC(:,2:x90_VOC_ind));
time_series.hpald = hpald;

end

end