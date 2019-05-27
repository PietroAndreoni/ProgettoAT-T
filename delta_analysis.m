function [infos] = delta_analysis(Dmat,tol,n1,n2)
%delta_analysis Seek major changes in a Delta matrix
%   Requires n1 and n2 as logical that state if raw and column of the
%   matrix respectively are by country or not (ie A -> n1,n2 = 1, B ->
%   n1=0, n2=1)
%   tol is the threshold below witch changes are considered irrelevant
%   output is [index1 index2 value countryrow secrow countrycol seccol]
%   if matrix is of type n1,n2 = (1,1) and 
%   [index1 index2 value countrycol seccol] if matrix is of type n1,n2 =
%   (0,1) (index1 in this case equals secrow)

global nTOT kTOT fTOT trTOT vTOT

for i=1:size(Dmat,1)
    for j=1:size(Dmat,2)
        if abs(Dmat(i,j)) <= tol
            Dmat(i,j) = 0;
        end
    end
end

[a,b] = find(Dmat);
c = zeros(size(a,1),1);
for i=1:size(a)
    c(i) = Dmat(a(i),b(i));
end
infos = [a b c];

if n1 ~= 0
    country1 = fix( a./ ( size(Dmat,1)/nTOT ) ) + 1; % + 1 is necessary to access correct country 
    sec1 = mod( a,( size(Dmat,1)/nTOT ) );
    for i=1:size(sec1,1)
        if sec1(i) == 0
            sec1(i) = 163;
        end
    end
    infos = [infos country1 sec1];  
end
    
if n2 ~= 0
    country2 = fix( b./ ( size(Dmat,1)/nTOT ) ) + 1;
    sec2 = mod( b,( size(Dmat,1)/nTOT ));  
    for i=1:size(sec2,1)
        if sec2(i) == 0
            sec2(i) = 163;
        end
    end
    infos = [infos country2 sec2];  
end
    
end

