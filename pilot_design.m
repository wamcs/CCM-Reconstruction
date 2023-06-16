function Q = pilot_design(it, C_est, Q_list, codebook)

[NA, NP, ~] = size(Q_list);
[~, M] = size(codebook); 

B = zeros(NA,NP);
for i=1:it
    B = B + 1/exp(10*(it-i))*Q_list(:,:,i);
end


basis  = get_basis(codebook, NP, M);

rank_C = rank(C_est); 
sigmas = zeros(NP,1);

if rank_C >= NP
    sigmas(1:NP) = 1;
else
    sigmas(1:rank_C) = 1;
end


sigmas = sort(sigmas, 'descend');

D = diag(sigmas);
V = basis;
R = V*D*V'; 

[U,D_1,~] = svd(C_est^(0.5));

U_1 = U(1:end,1:rank_C);
U_2 = U(1:end,rank_C+1:end);
U_2 = sum(U_2,2);
D_2 = D_1(1:rank_C,1:rank_C);

Sigma = zeros(NA,NP);
Sigma(1:NP,1:NP) = sqrt(D);

X = zeros(NA,NA);
X(1:end,1:rank_C) = U_1;


if rank_C <= NP
D = V'*B';
D = D(1:rank_C,1:rank_C);
D = D/U_1(1:rank_C,1:end)'/D_2;
[D_U,~,D_V] = svd(D);

I = zeros(NA,NA);
I(1:rank_C,1:rank_C) = D_U*D_V';
% I(1:rank_C,1:rank_C) = orth(rand(rank_C, rank_C) + 1i*rand(rank_C, rank_C));
else

I = zeros(NA,NA);
I(1:rank_C,1:rank_C) = orth(rand(rank_C, rank_C) + 1i*rand(rank_C, rank_C));
% I(rank_C+1:end,rank_C+1:end) = orth(rand(NA-rank_C, NA-rank_C) + 1i*rand(NA-rank_C, NA-rank_C));
end

X = X*I;
F = X*Sigma*V';


Q_bar = D_2\(U_1'*F);
Q = zeros(NA,NP);
Q_addtion = repmat(U_2,1,NP);
Q(1:rank_C,1:end) = U_1(1:rank_C,1:end)'\Q_bar;


Q = Q+Q_addtion;

if norm(Q'*C_est*Q-R,'fro') <= 1e-8
    disp('OK'); 
end


end
