function rankmat = normalrank(x_mat)
% 将 矩阵x_mat 的每一行转化成秩
x_mat = squeeze(x_mat);

x_mat = x_mat';
[m,n] = size(x_mat);
[~,rposimat] = sort(x_mat);
rankmat = zeros(m,n);
tempr = zeros(m,1);
for i = 1:n
    tempr(rposimat(:,i)) = 1:m;
    rankmat(:,i) = tempr;
end
rankmat = rankmat';
end