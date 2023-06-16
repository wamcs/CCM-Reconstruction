function Q0 = init_Q(C)
    
% This function is used to get the initial Q 

data_real = xlsread('Q0.xlsx', 2);
data_img = xlsread('Q0.xlsx', 3);
Qs = zeros(32, 8, 7);
max_rsrp = 0; 
Q0 = zeros(32, 8); 
for i=1:7
   start_idx = 2+9*(i-1); 
   end_idx = start_idx + 7; 
   Q = data_real(start_idx:end_idx, :) + 1i*data_img(start_idx:end_idx, :); 
   Q = Q'; 
   Qs(:, :, i) = Q;
   rsrp = trace(Q'*C*Q); 
   if rsrp > max_rsrp 
      Q0 = Q;
      max_rsrp = rsrp;
   end 
end
end 