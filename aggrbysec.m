function [Zaggr,Aaggr,faggr,Baggr,Raggr,Eaggr,vaggr] = aggrbysec(Z,A,f,B,R,E,v)
%AGGRBYSEC Create a world I/O by sector
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

Zaggr = zeros(kTOT,kTOT);
Aaggr = zeros(kTOT,kTOT);
faggr = zeros(kTOT,fTOT);
Baggr = zeros(trTOT,kTOT);
Raggr = zeros(trTOT,kTOT);
Eaggr = zeros(trTOT,kTOT);
vaggr = zeros(fTOT,kTOT);

for i=1:kTOT
    for j=1:kTOT
        for n1=1:nTOT
            
            if i<trTOT+1
                Baggr(i,j)= Baggr(i,j) + B(i,(n1-1)*trTOT + j);
                Raggr(i,j)= Raggr(i,j) + R(i,(n1-1)*trTOT + j);
                Eaggr(i,j)= Eaggr(i,j) + E(i,(n1-1)*trTOT + j);              
            end
            
            if i<fTOT+1
                vaggr(i,j)= vaggr(i,j) + v(i,(n1-1)*fTOT + j);
            end
            
            for n2=1:nTOT
            
                if j<fTOT+1
                    faggr(i,j)= faggr(i,j) + f((n1-1)*kTOT + i,(n2-1)*fTOT + j);
                end
                
                Zaggr(i,j)=Zaggr(i,j) + Z((n1-1)*kTOT + i,(n2-1)*kTOT + j);
                Aaggr(i,j)=Aaggr(i,j) + A((n1-1)*kTOT + i,(n2-1)*kTOT + j);
                
            end
        end
    end
end

end

