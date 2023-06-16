function basis = get_basis(codebook,NP,M)

% This function is used to get the basis which constructed the matrix R. 

anchor = randi(M);
B = codebook'*codebook;

index = find(abs(B(anchor,:))<1e-10);
if rank(B(:,[anchor,index])) ~= NP
    disp('there is an error in basis')
    disp(anchor)
end

A = codebook(:,[anchor,index]);

[~,~,E] = qr(A);
basis = A*E;
basis = basis(:,1:NP);
% fprintf('the rank of basis is %d.\n',rank(basis))
end 