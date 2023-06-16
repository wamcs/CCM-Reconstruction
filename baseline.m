function [ ratioc,rmses,C_list] = baseline(max_it, C)

% This function is the implementation of baseline algorithm which is
% described in the patent. 
% Input: 
%       type = 1, use codebook_i 
%       type =2, use codebook_ii 
% Output: 
%       rmses: RMSE lists 
%       ratioc: ratio lists 


% load C
C = C ./ trace(C);
[u1, d1]=eigs(C,1);

% codebook 
codebook = get_codebook(8); 

% parameter setting 
NA = 32; NP = 8; 
% max_it = 50; 

% Q part 
Q0 = init_Q(C); 
Q_list = get_mub(NA, NP, Q0);
 
v_list = zeros(max_it, 1); 
C_list_cqi = zeros(NA, NA, max_it); % cqi*Q*vm0*vm0'*q' 
rmses = [];  % rmses lists 
ratioc = [];
C_list = zeros(NA,NA,max_it);

for it = 1:max_it 
    idx = mod(it, 15); 
    if idx == 0
        idx = 15;
    end
    Q_t = Q_list(:,:, idx);
    [v_list(it), cqi] = get_feedback(Q_t, C, codebook); % feedback
    v_m0 = codebook(:, v_list(it)); 
    
    C_list_cqi(:, :, it) = cqi*Q_t*v_m0*v_m0'*Q_t'; 
    C_est = 1/it * sum(C_list_cqi, 3); 
    C_list(:,:,it) = C_est;
    [w, d_est] = eigs(C_est, 1);
    ratio =  real(w'*C*w/(d1));
    ratioc = [ratioc, ratio];
    C_est = C_est./norm(C_est,'fro');
    rmse = norm(C_est - C./ norm(C, 'fro'), 'fro');
    rmses = [rmses, rmse]; 
    if rmse <= 1e-4
        break; 
    end 
    fprintf('the ratio is %f at iteration %d \n', ratio, it);
    fprintf('the rmse is %f at iteration %d \n', rmse, it);
    
end 
end 