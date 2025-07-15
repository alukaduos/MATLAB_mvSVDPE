% clc;clear;close all;
if 1
    filepath_h = 'D:\Matlab\Work\Dongfeng2024\rollingdataset\Healthy\';
    filepath_f = 'D:\Matlab\Work\Dongfeng2024\rollingdataset\Faulty\';
    num_hf = 80;
    vib_data = zeros(5000,3,num_hf*2);flag_wptorvmd=1;fs=1000;
    for i = 1:num_hf
        filename_h = strcat(filepath_h,'H',num2str(i),'.mat');
        filename_f = strcat(filepath_f,'F',num2str(i),'.mat');
        load(filename_h);
        vib_data(:,:,i) = H;
        load(filename_f);
        vib_data(:,:,i+num_hf) = H;
    end
    vib_data = permute(vib_data,[2,1,3]);
    n_layers = size(vib_data,3);
    X_mts_cell = cell(n_layers,1);
    for i = 1:n_layers
        X_mts_cell{i} = squeeze(vib_data(:,:,i));
    end
    vib_label = [zeros(num_hf,1);ones(num_hf,1)];
    % waveletcoef_data = zeros(5000,2^decomp_level,200);
else
    filepath = "F:\MATLAB\syndata\";
    matname = "chaos_syndata_35_100group_10000.mat";
    filename = strcat(filepath,matname);
    load(filename);
    chaos_style_label = [0*ones(6,1);ones(6,1);2*ones(6,1);3*ones(6,1);4*ones(6,1);5*ones(4,1)];
    chaos_num_label = 1:34;
    for i = 1:34
        for j = 1:100
            temp_layer_indx = (i-1)*100+j;
            vib_data(:,1,temp_layer_indx) = 0
        end
    end

end

