function [v, cqi] = get_feedback(q_t, C, codebook)

% This function is used to get the feedback
% Input: 
%       C is the ground-truth channel convariance matrix
% Output: 
%       v is the precoder matrix indictor 

[~, M] = size(codebook);
max_value = 0;
index = 1; 
for m = 1:M
    v_m = codebook(:, m);
    value = v_m'*q_t'*C*q_t*v_m;
    if value > max_value
        max_value = value;
        index = m; 
    end 
end
v = index; 
cqi = abs(max_value);  
end 