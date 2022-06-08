function  OutputMatrixCSR(fileID, ia, ja, Kc, Mc, Cc)
% Output of the system matrices in CSR (Compressed Sparse Row) format

fprintf(fileID,'%s\n', "  n   nnz  nm");
fprintf(fileID,'%d %d %d \n', length(ia)-1, length(ja), nargin-3);
fprintf(fileID,'%s\n', "ia");
fprintf(fileID,'%d\n', ia);
fprintf(fileID,'%s\n', "ja");
fprintf(fileID,'%d\n', ja);

% Stiffness matrix
fprintf(fileID,'%s\n', "Kc");
fprintf(fileID,'%-23.16g\n',Kc);

if nargin == 5
	% Mass matrix
    fprintf(fileID,'%s\n', "Mc");
    fprintf(fileID,'%-23.16g\n',Mc);
end

if nargin == 6
	% Mass and damping matrices
    fprintf(fileID,'%s\n', "Mc");
    fprintf(fileID,'%-23.16g\n',Mc);
    fprintf(fileID,'%s\n', "Cc");
    fprintf(fileID,'%-23.16g\n',Cc);
end

end