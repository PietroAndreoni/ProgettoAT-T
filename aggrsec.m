function [mat_aggr] = aggrsec(mat,n1,n2,type)
%aggrsec Summary of this function goes here
%   Detailed explanation goes here
global nTOT kTOT


if type == "sec"
    aggr = kTOT;
elseif type == "co"
    aggr = nTOT;
else
    print("Error! Type sec for sectoral aggregation or co for country aggregation");
end

if n1 ~= 0 && n2 == 0
    mat_aggr = zeros(aggr, size(mat,2));
elseif n1 == 0 && n2 ~= 0 
    mat_aggr = zeros(size(mat,1), aggr);
elseif n1 ~= 0 && n2 ~= 0 
    mat_aggr = zeros(aggr,aggr);
else
    print("Error! pick the right indicators");
end

for i = 1:size(mat,1)
    for j = 1:size(mat,2)
        
        for nn = 1: aggr
        
            if n1 ~= 0 && n2 == 0 && i < size(mat,1)/aggr + 1
                mat_aggr(nn,j) = mat_aggr(nn,j) + mat( (nn-1)* size(mat,1)/aggr + i, j);
            end
            
            if n1 == 0 && n2 ~= 0 && j < size(mat,2)/aggr + 1
                mat_aggr(i,nn) = mat_aggr(i,nn) + mat(i, (nn-1)* size(mat,2)/aggr + j);
            end
            
            for nnn = 1: aggr            
                  
                if n1 ~= 0 && n2 ~= 0 && i < size(mat,1)/aggr + 1 && j < size(mat,2)/aggr + 1
                    mat_aggr(nn,nnn) = mat_aggr(nn,nnn) + mat((nn-1)* size(mat,1)/aggr + i, (nnn-1)* size(mat,2)/aggr + j);
                end

            end
        end       
    end
end

end

