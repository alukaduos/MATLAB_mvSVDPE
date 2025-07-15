function classLoss = svm_no_opt(parametersval)
Traindata = buildtraindata(parametersval);
Train_X = Traindata(:,1:end-1); Train_y = Traindata(:,end);
%SVMModel = fitcsvm(Train_X,Train_y,'Standardize',true,'KernelFunction','RBF',...    'KernelScale','auto','OptimizeHyperparameters','all');
% SVMModel = fitcsvm(Train_X,Train_y,'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',struct('AcquisitionFunctionName', 'expected-improvement-plus','Verbose',0,'ShowPlots',0));
SVMModel = fitcsvm(Train_X,Train_y);
CVSVMModel = crossval(SVMModel);
classLoss = kfoldLoss(CVSVMModel);
end


function Traindata = buildtraindata(parametersval)
loadhfdata;
num_totalsamples = length(X_mts_cell);
% simutimes = 10;
% pe_mts = zeros(m+2,simutimes);
pe_mts_cell = cell(num_totalsamples,1);

dim_pe_mts = zeros(num_totalsamples,1);



for i = 1:num_totalsamples
    % X_mts = normrnd(0,2,[m,n]);
    X_mts = X_mts_cell{i};
    % X_mts(end,:) = 0*X_mts(1,:)+0.02*rand([1,n]);
    pe_mts_cell{i} = calsvdpe_mts(X_mts,parametersval);
    dim_pe_mts(i) = length(pe_mts_cell{i});
    %i/num_totalsamples/2
end
dim_pe_vec = max(dim_pe_mts);
pe_mts = zeros(dim_pe_vec,num_totalsamples);
pe_mts_ori = pe_mts;
for i = 1:num_totalsamples
    if dim_pe_mts(i)<dim_pe_vec
        pe_mts(:,i) = padarray(pe_mts_cell{i},dim_pe_vec-dim_pe_mts(i),0,'post');
    else
        pe_mts(:,i) = pe_mts_cell{i};
    end
    % i/num_totalsamples/2+.5
end

Traindata = [pe_mts',vib_label];
end