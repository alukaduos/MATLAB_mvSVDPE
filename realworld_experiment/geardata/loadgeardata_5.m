clc;clear;
filepath = "D:\MATLAB\Work\MTS_Entropy\MATLAB_mvSVDPE\realworld_experiment\geardata\";% change the directory of geardata;
filelistset = {'Chipped_20_0.csv';'Health_20_0.csv';'Miss_20_0.csv';'Root_20_0.csv';'Surface_20_0.csv';...
    'Chipped_30_2.csv';'Health_30_2.csv';'Miss_30_2.csv';'Root_30_2.csv';'Surface_30_2.csv'};
no_test=5;
filelist = cell(no_test,1);
for i = 1:no_test
    filelist{i} = strcat(filepath,filelistset{i});
end
length_sample = 3e3;num_samples = 200; num_groups = 5;
X_mts_cell = cell(num_groups*num_samples,1);gearlabel = zeros(num_groups*num_samples,1);
startpoints = 1:length_sample:num_samples*length_sample;
countsamples = 0;
for i=1:no_test
    tempdata = readtable(filelist{i});
    tempdata = table2array(tempdata(12:end,1:3));
    for j = 1:num_samples
        countsamples = countsamples+1;
        gearlabel(countsamples) = mod(i,5);
        X_mts_cell{countsamples} = tempdata(startpoints(j):startpoints(j)+length_sample-1,:)';
        % if mod(i,5)~=2
        %     if j>50
        %         break;
        %     end
        % end
    end
end
num_totalsamples = countsamples;
% save geardata_5.mat