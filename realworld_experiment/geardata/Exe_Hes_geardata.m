clc;clear;close all;
loadgeardata_5;
tau_1_set = 1:3;         nn1 = length(tau_1_set);
tau_2_set = 3:2:13;        nn2 = length(tau_2_set);
dvalue_set = 5:2:20;       nn3 = length(dvalue_set);
d_pvalue_set = 3:2:10;     nn4 = length(d_pvalue_set);
Traindataset = cell(nn1,nn2,nn3,nn4);
countcompletes = 0;     nnall = nn1*nn2*nn3*nn4;
% it takes a very long duration to compute the following 4-layers cycles
for i1 = 1:nn1
    for i2 = 1:nn2
        for i3 = 1:nn3
            for i4 = 1:nn4
                % tic
                tau_1 = tau_1_set(i1); tau_2 = tau_2_set(i2);
                d = dvalue_set(i3); d_p = d_pvalue_set(i4);
                parametersval = [tau_1,tau_2,d,d_p];
                Traindataset{i1,i2,i3,i4} = calsvdtrain(parametersval);
                countcompletes = countcompletes+1;                
                countcompletes/nnall*100
                % toc
            end
        end
    end
end
save traindata_gear_5;

function Traindata = calsvdtrain(parametersval)
load geardata_5.mat
pe_mts_cell = cell(num_totalsamples,1);
dim_pe_mts = zeros(num_totalsamples,1);
for i = 1:num_totalsamples
    X_mts = X_mts_cell{i};
    pe_mts_cell{i} = calsvdpe_mts(X_mts,parametersval);
    dim_pe_mts(i) = length(pe_mts_cell{i});
    %i/num_totalsamples/2
end
dim_pe_vec = max(dim_pe_mts);
pe_mts = zeros(dim_pe_vec,num_totalsamples);
for i = 1:num_totalsamples
    if dim_pe_mts(i)<dim_pe_vec
        pe_mts(:,i) = padarray(pe_mts_cell{i},dim_pe_vec-dim_pe_mts(i),0,'post');
    else
        pe_mts(:,i) = pe_mts_cell{i};
    end
    %i/num_totalsamples/2+.5
end
Traindata = [pe_mts',gearlabel];

end

