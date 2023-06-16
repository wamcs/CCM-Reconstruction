function Q_list = get_mub(NA, NP, Q0)

% This function is get the Q_list which is used in the baseline algorithm. 
% Input: 
%   Q0 is the initial Q 
% Q = Q0 * Mub matrix 

load('mub.mat'); 
K8 = zeros(8, 8, 15);
[~, ~, num_k2] = size(K2); 
[~, ~, num_k4] = size(K4);

for i = 1:num_k2
    for j = 1:num_k4 
        K8(:, :, (i-1)*num_k4+j) = kron(K2(:, :, i), K4(:, :, j));  
    end 
end

Q_list = zeros(NA, NP, 15); 
for i = 1:15
    Q_list(:, :, i) = Q0 * K8(:, :, i); 
end 

end