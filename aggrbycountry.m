function [Zaggr,Aaggr,faggr,Baggr,Raggr,Eaggr,vaggr] = aggrbycountry(Z,A,f,B,R,E,v)
%AGGRBYCOUNTRY Create a country I/O with all sectors aggregated
%   Requires dimensions defined as global variables
%   Requires all input insertion

global nTOT kTOT fTOT trTOT vTOT

%{
if size(Z,1) ~= size(Z,2) || size(Z,1)  ~= kTOT*nTOT
    disp("(at least) matrix A is incorrect in size (square matrix required)")
    return
end

if size(B,2) ~= size(Z,2) || size(B,1)  ~= trTOT
        disp("(at least) matrix B is incorrect in size")
    return
end

if size(f,1) ~= size(Z,1) || size(f,2) ~= fTOT*nTOT
        disp("(at least) matrix f is incorrect in size")
    return
end

if size(v,2) ~= size(Z,1) || size(v,2) ~= vTOT
        disp("matrix v is incorrect in size")
    return
end    
%}

Zaggr = zeros(nTOT,nTOT);
Aaggr = zeros(nTOT,nTOT);
faggr = zeros(nTOT,nTOT);
Baggr = zeros(trTOT,nTOT);
Raggr = zeros(trTOT,nTOT);
Eaggr = zeros(trTOT,nTOT);
vaggr = zeros(fTOT,nTOT);

for n1=1:nTOT
    for n2=1:nTOT
        for i=1:kTOT
            for j=1:kTOT
                
                Zaggr(n1,n2) = Zaggr(n1,n2) + Z((n1-1)*kTOT + i, (n2-1)*kTOT + j);
                Aaggr(n1,n2) = Aaggr(n1,n2) + A((n1-1)*kTOT + i, (n2-1)*kTOT + j);
                
                if j < fTOT+1
                    faggr(n1,n2) = faggr(n1,n2) + f((n1-1)*kTOT + i, (n2-1)*fTOT + j);
                end
                
                if n1 == 1
                    
                if i < trTOT+1
                    Baggr(i,n2) = Baggr(i,n2) + B(i, (n2-1)*kTOT + j);                   
                    Raggr(i,n2) = Raggr(i,n2) + R(i, (n2-1)*kTOT + j);             
                    Eaggr(i,n2) = Eaggr(i,n2) + E(i, (n2-1)*kTOT + j);
                    
                end
                
                if i < vTOT+1
                    vaggr(i,n2) = vaggr(i,n2) + v(i, (n2-1)*vTOT + j);
                end 
                
                end
                
            end
        end
    end
end
            

end

