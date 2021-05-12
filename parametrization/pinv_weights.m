%%
% Import tables containing information about individual experiments and
% yields of measured analytes
conditions = importdata('conditions_matrix_updated.csv');
yields = importdata('yields_matrix_updated.csv');
%% 
% Store the parameters 
parameters = {'2_bromophenylpyrroline','methyl acrylate','decamethylcobaltocene','Lewis acid','Ni_cod','water_workup','bicarb_workup','ammonium_chloride_workup','sodium_carbonate_workup','naoh_workup','BEt3','AlMe3','AlEt3','toluene','benzene','trifluorotoluene','thf','ether', 'xylene'}

%% 
% Extract the numbers from data tables
conditions = conditions.data;
yields = yields.data;
%% 
% Remove experiments where mass balance was less than 75 and greater than
% 120%
yields_trimmed = yields(yields(:,6)>75 & yields(:,6)<120, :);
conds_trimmed = conditions(yields(:,6)>75 & yields(:,6)<120, :);
parameters(sum(conds_trimmed)==0)=[];
%% 
% Center the data by subtracting the mean and dividing by the mean of each
% column
conds_trimmed(:, sum(conds_trimmed)==0) = [];
conds_centered = conds_trimmed - ones(size(conds_trimmed,1),1)*sum(conds_trimmed)/size(conds_trimmed, 1);
yields_centered = yields_trimmed - ones(size(yields_trimmed,1),1)*sum(yields_trimmed)/size(yields_trimmed, 1);
%% 
% Calculate pseudoinverse to get contributions of each parameter to yields
centered_weights = pinv(conds_centered)*yields_centered;
%% 
% Plot for quick visual
plot(centered_weights)

legend('anti', 'syn', 'phenyl', '2-bromophenyl', '2-ethyl')
%% 
fileID = fopen('centered_weights.txt','w');
fprintf(fileID,'%f,%f,%f,%f,%f,%f\n',centered_weights')
fclose(fileID)

