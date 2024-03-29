function Q = SearchMW(mw,Cnames,MCMversion)
% function Q = SearchMW(mw,Cnames, MCMversion)

% NOTE: This will NOT work for
%       1) Inorganic species (e.g. OH, NO, etc.)
%       2) any non-MCM molecules added to your species list (e.g. extra monoterpenes)
%
% INPUTS:
% mw: molecular weight
% Cnames:   a cell array of MCM names to search.
% MCMversion: string specifying which version of MCM to use, 'v32' or 'v331'.
%
% OUTPUTS:
% Q: a structure with the following fields:
%   names:  cell array of MCM names for species
%   index:  vector of indices that map Q.names to the input Cnames
%   MW:     vector of molecular weights
%   nX:     number of a given atom, where X = C,O,N,H,S,F,Cl,Br,I
%
% HSK adapted to search for molecular weight (based on searchSMILES)
% 20130104 GMW
% 20210303 GMW  Added alkoxyRadical option
% 20210914 GMW  Added peroxyacylRadical option (same as APRadical)
%               Added aryloxyRadical group
% 20210915 GMW  Updated alkoxyRadical pattern to catch non-primary radicals ((C)[O] and similar).
%
%%%%%DIRECTIONS FOR CREATING MCMSpeciesInfo.mat%%%%%
% 1) Go to MCM website and select ALL VOCs
% 2) Go to "Extract", select "Molecular Weights for all species in a subset" and download
% 3) Read textfile into matlab with the following (replacing filename and number of header
%   lines as needed):
%   [MCMnames,SMILES,InChI,MolWeight] = textread(flnm,'%s %s %s %f','headerlines',nhead);
%    save MCMvxxx_SpeciesInfo.mat

%%%%%FUNCTIONAL GROUP REGULAR EXPRESSIONS%%%%%


%%%%%GRAB SMILES STRINGS%%%%%
% Variables are:
% MCMnames
% SMILES
% MolWeight
switch MCMversion
    case 'v32'
        load MCMv32SpeciesInfo.mat
    case 'v331'
        load MCMv331SpeciesInfo.mat
end


[tf,loc] = ismember(Cnames,MCMnames);
S = cell(size(Cnames));
S(~tf) = {''};
S(tf) = SMILES(loc(tf));
M = nan(size(MolWeight));
M(tf) = MolWeight(loc(tf));

%%%%%DO THE FINDING%%%%%
match = regexp(S,'.*');
Q.index = find(~cellfun('isempty',match)); %index for matching names
Q.names = Cnames(Q.index);
Q.MW = M(Q.index);

%%%%%COUNT ATOMS%%%%%
Q.nC = cellfun('length',regexp(S(Q.index),'C(?!l)'));
Q.nO = cellfun('length',regexp(S(Q.index),'O'));
Q.nN = cellfun('length',regexp(S(Q.index),'N'));
Q.nS = cellfun('length',regexp(S(Q.index),'S'));
Q.nF = cellfun('length',regexp(S(Q.index),'F'));
Q.nCl = cellfun('length',regexp(S(Q.index),'Cl'));
Q.nBr = cellfun('length',regexp(S(Q.index),'Br'));
Q.nI = cellfun('length',regexp(S(Q.index),'I'));

%estimate hydrogen as residual mass
mH = 1.00794;   mC = 12.0107;   mO = 15.9994;   mN = 14.0067;   mS = 32.065;
mF = 18.998403; mCl = 35.453;   mBr = 79.904;   mI = 126.90447;
Q.nH = round((Q.MW - (Q.nC*mC + Q.nO*mO + Q.nN*mN + Q.nS*mS + ...
    Q.nF*mF + Q.nCl*mCl + Q.nBr*mBr + Q.nI*mI))/mH);

Qfields = fields(Q);

for i = 1:length(Qfields)
    Q.(Qfields{i}) = Q.(Qfields{i})(abs(Q.MW-mw)<0.0001);
end

