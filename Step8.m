%Matlab file for generating the X matrix with dimension 21578*8514 
X=zeros(21578,8514);
[rowp,colp]=size(X);
for i=1:rowp
	base1='C:\Users\ccuser\Documents\MATLAB\tfidffiles\tfidf-';
	extension1='.dat';
   	filenumber1=sprintf('%05d',i);
    fe=strcat(filenumber1,extension1);
   	base=strcat(base1,fe)
    file=fopen(base,'r');
    sizeA = [8514 Inf];
   A=fscanf(file,'%d %f',sizeA);
   A=transpose(A);
   j=1;
   jp=0;
   [row,col]=size(A);

    while j <=col
   
        jp=j+1;
    X(i,A(1,j))=A(1,jp);
    j=j+2;
    end
   fclose(file);
end
X=transpose(X);
%[U,Sigma,Vt]=svd(X,0);
% tp=transpose(A);
% t2p=inv(S);
%V=transpose(X)*U*inv(S);
% fileID = fopen('type1.txt','r');
% formatSpec = '%f %d';
% sizeA = [2 Inf];
% A = fscanf(fileID,formatSpec,sizeA)
% fprintf(fileID,'%d %4.4f\n',y);
% transpose(A)
% fclose(fileID);
