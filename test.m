core_number = 10;
parpool('local',core_number);


N = 100;
T = 35;
rt = zeros(N,T);
rm = zeros(N,T);


NA = 32;
NP = 8;

rt_patent = zeros(N,T);
rm_patent = zeros(N,T);

 
load('data_2_ports.mat');
C_list = zeros(NA,NA,N);
C_est_list = zeros(N,NA, NA,T);
vio_list = zeros(N,T);
C_est_list_patent = zeros(N,NA,NA,T); 

for n =1:N
 
    H=Hall(:,:,1,n,1);
    C = H'*H;
    C_list(:,:,n) = C;
     [t1,t2] = benchmark(H,1);
     Type_ii(n) = t2;
     Type_i(n) = t1;
end

parfor i = 1:N
    [r1,r2,vios,c_ls] = main(C_list(:,:,i),T,2);
    rt(i,:) = r1;
    rm(i,:) = r2;
    C_est_list(i,:,:,:) = c_ls;
    vio_list(i,:) = vios;
    [r1_p, r2_p,c_ls] = baseline(T, C_list(:,:,i));
    rt_patent(i,:) = r1_p;
    rm_patent(i,:) = r2_p;
     C_est_list_patent(i,:,:,:) = c_ls;

end
delete(gcp('nocreate'));
