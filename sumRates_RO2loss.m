function rates = sumRates_RO2loss(S)
% function rates = sumRates_ROx(S)
% Calculates rates loss for RO2

% LOSS REACTION DEFINITIONS
% Lnames = {'HO2 + RO2';'RO2 + RO2';'RO2 + NO';'RO2 + NO2';
%           'acylRO2 + NO2';'RO2 + NO3';'other'};

% 
% INPUTS
% S is a F0AM model structure.
%
% OUTPUTS
% rates is structure containing the following fields.
% rates.Loss is a matrix of summed production rates with dimensions of (# model stepx) x (# reactions)
% rates.Lnames is a cell array containing the names of all reaction families.
% rates.iRx_Loss is a cell array containing the indices of all reactions included in a family.
%   You can use this to access those reactions in S.Chem, e.g. 
%   Rnames_oVOChv = S.Chem.Rnames(rates.iRx_Prod{4});
%   This might be handy if you want to look at major contributors to a family.


% Get net P/L rates
RO2names = S.Cnames(S.Chem.iRO2);
if ismember('CISOPA', S.Cnames)
    Fnames = [{'ROx'}; RO2names; {'CISOPA';'CISOPC';'TISOPA';'TISOPC';'HCOCO'}];
else
    Fnames = [{'ROx'}; RO2names];
end
Frates = PlotRates(Fnames,S,5000,'sumEq',0,'plotme',0);

%% Loss

% get species indices
[~,iHO2] = ismember('HO2',S.Cnames); %location of HO2
[~,iRO2] = ismember(RO2names,S.Cnames); %location of RO2 (speciated)
[~,iNO] = ismember('NO',S.Cnames); % NO index
[~,iNO2] = ismember('NO2',S.Cnames);
APinfo = SearchSMILES('peroxyacylRadical',S.Cnames,'v331');
[~,iNO3] = ismember('NO3',S.Cnames);

% index reactions for each species
iG = S.Chem.iG(Frates.iRx_Loss,:); %reactant indices
jHO2 = any(iG == iHO2,2);
jRO2sum = any(iG == 2,2); %sum of RO2 is always second species
jRO2speciated = any(ismember(iG,iRO2),2);
jNO = any(iG == iNO,2);
jNO2 = any(iG == iNO2,2);
jAP = any(ismember(iG,APinfo.index),2); %peroxyacyl radicals
jNO3 = any(iG == iNO3,2);

% index full reactions
Lnames = {'HO2 + RO2';'RO2 + RO2';'RO2 + NO';'RO2 + NO2';...
    'acylRO2 + NO2';'RO2 + NO3';'other'};
jHO2RO2 = jHO2 & jRO2speciated;
jRO2RO2 = jRO2sum & jRO2speciated;
jRO2NO  = jRO2speciated & jNO;
jRO2NO2 = jRO2speciated & jNO2;
jacylRO2NO2 = jAP & jNO2;
jRO2NO3 = jRO2speciated & jNO3;

jother  = ~(jHO2RO2 | jRO2RO2 | jRO2NO | jRO2NO2 | jacylRO2NO2 | jRO2NO3);

% order indices
jall = {jHO2RO2;jRO2RO2;jRO2NO;jRO2NO2;jacylRO2NO2;jRO2NO3;jother};

% loop through, sum and store reaction index
rates.Lnames = Lnames;
rates.Loss = zeros(length(S.Time),length(Lnames));
rates.iRx_Loss = cell(length(Lnames),1);
for i = 1:length(Lnames)
    if ~isempty(jall{i})
        rates.Loss(:,i) = sum(Frates.Loss(:,jall{i}),2);
    end
    rates.iRx_Loss{i} = Frates.iRx_Loss(jall{i});
end



% that's all

%% Old ideas
% other alternate, use of names
    % get reactant names
%     [Lnames1,Lnames2] = strtok(Frates.Lnames,'+'); %extra step here protects against typos in reaction names
%     Lnames2 = strtok(Lnames2,'+');
%     Lnames1 = strtrim(Lnames1);
%     Lnames2 = strtrim(Lnames2);
    
    % alternate use of f's
%     fHO2 = full(Snow.Chem.f(Frates.iRx_Loss,iHO2)); %stoich coeff for HO2
%     fRO2 = Snow.Chem.f(Frates.iRx_Loss,iRO2); %stoich coeff for RO2
%     fRO2 = full(sum(fRO2,2)); %sum along species
%     jHO2RO2 = fHO2 == -1 & fRO2 == -1;

