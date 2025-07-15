load chaos_syndata_5_300000.mat
length_sample_each = 3000;
num_chaos = size(chaos_syndata,1);
length_total = size(chaos_syndata{1},1);
num_samples = fix(length_total/length_sample_each);
num_totalsamples = num_chaos*num_samples;
dim_mts = zeros(num_totalsamples,1);

X_mts_cell = cell(num_totalsamples,1);labels = zeros(num_totalsamples,1);
for i_c = 1:num_chaos
    temp_totaldata = chaos_syndata{i_c};
    
    for i_s = 1:num_samples
        X_mts_cell{(i_c-1)*num_samples+i_s} = (temp_totaldata((i_s-1)*length_sample_each+1:i_s*length_sample_each,:))';
        labels((i_c-1)*num_samples+i_s) = i_c;
        dim_mts((i_c-1)*num_samples+i_s)  = size(chaos_syndata{i_c},2);
    end
end

%save chaos_syndata_5_3000_100.mat X_mts_cell labels num_totalsamples num_samples