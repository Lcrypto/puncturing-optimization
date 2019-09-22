function CWeff = CWeff(g, R, H)
%CWeff compute effective column weight of g's column over R subset of rows
%   Detailed explanation goes here
col = H(:, g);
Lambda = find(col); % nonzero positions in column
Nonzeros = intersect(Lambda, R); % effective nonzero positions
CWeff = length(Nonzeros);

end

