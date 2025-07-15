clc;clear;
filepath = "D:\MATLAB\Work\MTS_Entropy\MATLAB_mvSVDPE\realworld_experiment\CWRU\cwrudataset\";
listname = 'filenamelist.xlsx';
filelist = readtable(strcat(filepath,listname));
length_sample = 1e3;num_samples = 120; num_groups = 4;
X_mts_cell = cell(num_groups*num_samples,1);cwrulabel = zeros(num_groups*num_samples,1);
startpoints = 1:length_sample:num_samples*length_sample;

countsamples = 0;
for i = 1:8
    temp_matname = filelist.matname{i};
    tempfilename = strcat(filepath,temp_matname);
    load(tempfilename);
    tempdata = [eval(filelist.datamatname2{i}),eval(filelist.datamatname3{i})];
        
    for j = 1:num_samples
        countsamples = countsamples+1;
        cwrulabel(countsamples) = mod(i,4);
        X_mts_cell{countsamples} = tempdata(startpoints(j):startpoints(j)+length_sample-1,:)';

    end
end
num_totalsamples = countsamples;
