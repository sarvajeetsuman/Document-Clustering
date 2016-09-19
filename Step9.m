%It decomposes X matrix into 3 matrices U , sigma , Vt having sigma dimension = 64 
size(X);
 [U,Sigma,Vt]=svds(X,64,'L');
% tp=transpose(A);
% t2p=inv(S);
V=transpose(X)*U*inv(Sigma);
