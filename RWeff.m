function RWeff = RWeff(p, C, H)
%RWeff compute effective row weight of p's row over C subset of columns
%   Detailed explanation goes here
row = H(p, :);
Gamma = find(row); % nonzero positions in row
Nonzeros = intersect(Gamma, C); % effective nonzero positions
RWeff = length(Nonzeros);
end

