length_sample_each = n;
num_signals = 9;
num_samples = 100;
indx_hurst = 0.7; coef_r = [0 0 0 0 0 -0.8 -0.3 0.3 0.8];
total_samples = num_signals*num_samples;
X_mts_cell = cell(total_samples,1);mtslable = zeros(total_samples,1); indx_signal = 0;
for i_signal = 1:num_signals
    for i_samples = 1:num_samples
        indx_signal = indx_signal + 1;
        tempmts = createmanystoc(length_sample_each,i_signal,indx_hurst,coef_r(i_signal));
        X_mts_cell{indx_signal} = tempmts;
        mtslable(indx_signal) = i_signal;
    end
end
num_totalsamples = total_samples;

function tempmts = createmanystoc(length_samples,i_signal,indx_hurst,coef_r)
tempmts = zeros(2,length_samples);
if i_signal ==1
    tempmts = rand(2,length_samples);
elseif i_signal == 2
    tempmts = normrnd(0,1,2,length_samples);
elseif i_signal == 3
    tempmts(1,:) = rand(1,length_samples);
    tempmts(2,:) = normrnd(0,1,1,length_samples);
elseif i_signal == 4
    tempmts(1,:) = rand(1,length_samples);
    tempmts(2,:) = wfbm(indx_hurst,length_samples);
elseif i_signal == 5
    tempmts(1,:) = normrnd(0,1,1,length_samples);
    tempmts(2,:) = wfbm(indx_hurst,length_samples);
else
    mu_mat = zeros(1,2); sigma_mat = [1 coef_r;coef_r,2];
    temp = mvnrnd(mu_mat,sigma_mat,length_samples);
    tempmts = temp';
end
end
