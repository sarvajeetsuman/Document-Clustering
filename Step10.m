%Write V matrix in to myfile_final.txt
fid = fopen('myfile_final16.txt', 'wt'); % Open for writing
[rowh,colp]=size(V)
for i=1:rowh
   fprintf(fid, '%f ', V(i,:));
   fprintf(fid, '\n');
end
fclose(fid);
