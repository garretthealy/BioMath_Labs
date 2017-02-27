function [lat,objvalue] = lattice(seq,index)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Solve the lattice structure problem associated with Lab 5
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2
    error('Two inputs are needed, see function code for more info')
end

s = size(index);

m = s(1); 
n = s(2);

lat = char(m,n);

for i = 1:m
    for j = 1:n
        ind = index(i,j);
        lat(i,j) = seq(ind);
    end
end

objvalue=0;
tlat = lat;

for i = 1:m
    for j = 1:n
        if lat(i,j)=='H'
            if i~=m && lat(i+1,j)=='H'
                objvalue = objvalue - 1;
            end
            if j~= n && lat(i,j+1)=='H'
                objvalue = objvalue - 1;
            end
        end
    end
end


end