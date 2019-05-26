function [infos] = delta_analysis(Dmat,tol,n1,n2)
%delta_analysis Seek major changes in a Delta matrix
%   Requires n1 and n2 as logical that state if raw and column of the
%   matrix respectively are by country or not (ie A -> n1,n2 = 1, B ->
%   n1=0, n2=1)
%   tol is the threshold below witch changes are considered irrelevant

global nTOT kTOT fTOT trTOT vTOT

for i=1:size(Dmat,1)
    for j=1:size(Dmat,2)
        if abs(Dmat(i,j)) <= tol
            Dmat(i,j) = 0;
        end
    end
end

[a,b] = find(Dmat);
infos = [a b];

if n1 ~= 0
    country1 = fix( a./kTOT ) + 1; % + 1 is necessary to access correct country
    sec1 = mod( a,kTOT);
    infos = [infos country1 sec1];  
end
    
if n2 ~= 0
    country2 = fix( b./kTOT ) + 1;
    sec2 = mod( b,kTOT);  
    infos = [infos country2 sec2];  
end
    
end

