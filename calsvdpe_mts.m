function pe_mts = calsvdpe_mts(X_mts,parametersval)
tau_1 = parametersval(1);tau_2 = parametersval(2);d = parametersval(3); d_p = parametersval(4);
[m,n] = size(X_mts);
k_max = fix((n - (d-1)*tau_2)/tau_1);
Y = zeros(m,k_max,d);%3D-rearrangements of X_mts, dimension*k_max*emd_dim
for k = 1:k_max
    for j = 0:d-1
        Y(:,k,j+1) = X_mts(:,k*tau_1+j*tau_2);
    end
end
Y_p = permute(Y,[1,3,2]); %dimension*emd_dim*k_max
s_v = squeeze(pagesvd(Y_p)); %singular val of Y_p
[U_v,~,V_v] = pagesvd(Y_p); %orthogonal U and V
%V_v = V_v';
U_dim = size(U_v,1); V_dim = size(V_v,1);dis_cos_u = zeros(1,k_max); dis_cos_v = dis_cos_u;
for k = 1:k_max
    temp_u = U_v(:,:,k);
    mean_u = mean(temp_u,2);
    mean_u_eye = sum(eye(U_dim),2)/U_dim;
    temp_v = V_v(:,:,k);
    mean_v = mean(temp_v,1)';
    mean_v_eye = sum(eye(V_dim),2)/V_dim;
    dis_cos_u(k) = dot(mean_u, mean_u_eye)/norm(mean_u)/norm(mean_u_eye);
    dis_cos_v(k) = dot(mean_v, mean_v_eye)/norm(mean_v)/norm(mean_v_eye);

end
S_v = [dis_cos_u;dis_cos_v;s_v];
i_max = size(S_v,2)-d_p+1;%num_samples_of_permutations
m_p = min(m,d);
perm = zeros(m_p+2,d_p,i_max);rank_perm = perm;
for i = 1:i_max
    perm(:,:,i) = S_v(:,i:i+d_p-1);
   % rank_perm(:,:,i) = tiedrank(perm(:,:,i)')';
    rank_perm(:,:,i) = normalrank(perm(:,:,i));
end
pe_mts = calpe_permat(rank_perm);
end
function pe_val = calpe_permat(x)%each row a permutation, rows=(dis_U;dis_V;S_1;S_2;......); each layer an embedded mat
num_feature = size(x,1);
pe_val = zeros(num_feature,1);
for i = 1:num_feature
    temp_x = squeeze(x(i,:,:));
    pe_val(i) = calpe_permarray(temp_x);
end
end
function pe_val = calpe_permarray(x)
% [len_perm,num_sample] = size(x);
[~,~,idx_perm] = unique(x','rows');
idx_counts = accumarray(idx_perm,1);
perm_ratio = idx_counts/size(x,2);
pe_val = -sum(log2(perm_ratio).*perm_ratio);
end
