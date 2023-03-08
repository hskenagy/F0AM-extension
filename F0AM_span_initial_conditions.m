% This is the script that generates the initial condition matrix, runs a 
% F0AM simulation for each of the initial conditions in that matrix, and 
% then tabulates the calculated RO2-chemistry-related parameters.  This is 
% the file where you specify which initial conditions you want to vary, 
% and over what ranges.  It’s also where you specify the usual met 
% variables and chem files to use.  And, at the end of the file, you can 
% specify what to call the tabulated output for saving.


% start by clearing all variables in the workspace (so there isn't any
% confusion and also because this thing is about to use a lot of memory)
clear all

% set the dimensions for your initial condition matrix
m = 10;
n = 10;
o = 4;

% specify the ranges for the initial condition matrix
HONO_span = 5*logspace(0,3,m);
H2O2_span =  5*logspace(0,4,n);
NO_span = logspace(-3,3,o);

% generate a matrix of initial conditions, then vectorize it
[HONO_matrix, H2O2_matrix, NO_matrix]=ndgrid(HONO_span, H2O2_span, NO_span);
HONO_concs = reshape(HONO_matrix,[],1);
H2O2_concs = reshape(H2O2_matrix,[],1);
NO_concs = reshape(NO_matrix,[],1);


%% METEOROLOGY

Met = {...
%   names       values          
    'P'         1013                       ; %Pressure, mbar
    'T'         298                        ; %Temperature, K
    'RH'        1                         ; %Relative Humidity, percent
    'LFlux'     '340nm_update.txt'     ; %Text file for radiation spectrum
    'jcorr'     2*3.3e-4                      ; %Adjusted for Fall 2020 experiments based off of jNO2 of 0.06 min^-1
    'kdil'      0;%1/(7500*60)                        ; %dilution rate in s^-1
    };

%% CHEMICAL CONCENTRATIONS

InitConc = {...
%   names       conc(ppb)           HoldMe
   'C5H8'    100                     0;
%    'APINENE' 100                     0;
   'NO'      NO_concs                0;
   'NO2'     0                       0;
   'HONO'    HONO_concs              0;
   'H2O2'    H2O2_concs              0;
    };
%% CHEMISTRY

ChemFiles = {...
    'MCMv331_K(Met)';
    'MCMv331_J(Met,1)'; %Jmethod flag of 1 specifies using "BottomUp" J-value method.
%     'MCMv331_AllRxns';
     'MCMv331_Inorg_Isoprene';
%     'MCMv331_Inorg_apinene'
%    'CH3ONO_hv'; %not included in MCM (doesn't exist in the atmosphere)
    };

%% DILUTION CONCENTRATIONS

BkgdConc = {...
%   names           values
    'DEFAULT'       0;   %0 for all zeros, 1 to use InitConc
   };

%% OPTIONS

ModelOptions.Verbose       = 1;
ModelOptions.EndPointsOnly = 0;
ModelOptions.LinkSteps     = 0;
%ModelOptions.Repeat        = 1;
ModelOptions.IntTime       = 3600*4; %Total model time in seconds
ModelOptions.SavePath      = 'IterateInitialConditionsOutput.mat';
ModelOptions.GoParallel    = 0;

%% MODEL RUN

S = F0AM_ModelCore(Met,InitConc,ChemFiles,BkgdConc,ModelOptions);

%% PLOTTING

SplitRun(S,'step')
Splot = who('-regexp','S\d');
%%
labels = cell(length(Splot),1);
cond_table = table;
time_series = struct;

% here you must specify which VOC precursor you're using and which
% intermediate species you want to use for the calculations (e.g. MVK for
% C5H8)
for i = 1:length(Splot)
    [cond_table(i,:), time_series.(strcat(['S' num2str(i)]))] = calculate_RO2_parameters(eval(['S' num2str(i)]), 0, 'C5H8', 'MVK');
%     [cond_table(i,:), time_series.(strcat(['S' num2str(i)]))] = calculate_RO2_parameters(eval(['S' num2str(i)]), 0, 'APINENE', 'APINBOOH');
end

% the cond_table is a table with a row for each experiment and tabulates all of
% the important RO2-related parameters for each experiment (as specified in
% calculate_RO2_parameters)
writetable(cond_table);
% the time_series structure contains a field for each experiment, each of
% which contains a couple important RO2-related time series (as specified
% in calculate_RO2_parameters)
save('time_series.mat', 'time_series')
   
% remove the F0AM output from your workspace so your computer's memory
% doesn't explode
clearvars S*
