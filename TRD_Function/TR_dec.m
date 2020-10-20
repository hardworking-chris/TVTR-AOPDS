function X  =  TR_dec(G,u,Omega)
N = ndims(G);
for n = 1:N
	B = tenmat_sb(Z_neq(G,n),2);  B=B';% B is the right part of the right part of the relation equation
    Omega = tenmat_sb(Omega,2);
    Gn = Gunfold(G{n},2);
    u = tenmat_sb(u,2);
    I_n = size(G{n},2);
    parfor i=1:I_n        
        Gn(i,:) = (u(i,:).*Omega(i,:)) * pinv(B*repmat(Omega(i,:),size(B,2),1));
    end
	G{n} = Gfold(Gn,size(G{n}),2);
end

% update X
X = G;
