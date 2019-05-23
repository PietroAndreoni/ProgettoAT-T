%% Basic calculations
clear ;
load("base_students.mat");

% number of countries, sectors, final demand types and exogenous resources
% in the model
nTOT = 48;
kTOT = 163;
fTOT = 7;
trTOT = 12;

%calculate Leontieff coefficients matrix and overall production
L = inv(diag(ones(size(A,1),1))-A);
x0 = L*f;
x = sum(x0,2);
Z = A*diag(x);

%calculate aggregate final demand for sector i in country A
f_tot=zeros(nTOT*kTOT,1);
for i = 1:nTOT*kTOT
    f_tot(i) = sum( f(i,:) );
end

%calculate exogenous transaction matrix and embodied resource matrix
R = B* ( diag(L*f_tot) );
E = (B*L) * diag(f_tot);

%% Some extra calculations about final demand

%calculate aggregate final demand for sector i in country A FROM country B
%(supply side)
f_totS = zeros(nTOT*kTOT,nTOT);
for i = 1:nTOT*kTOT
    for n = 1:nTOT 
        a = (n-1)*fTOT+1 : (n*fTOT);
        f_totS(i,n) = sum( f(i,a) );
    end
end

%calculate aggregate final demand per country (GDP) supply side
GDP = zeros(nTOT,1);
for n = 1:nTOT
    GDP(n) = sum(f_tot((n-1)*kTOT+1:n*kTOT));
end


%calculate aggregate final demand for sector i in country A FROM country B
%(demand side)
f_totD = zeros(nTOT,fTOT*nTOT);
for n = 1:nTOT 
    for i = 1:nTOT*fTOT

        a = (n-1)*kTOT+1 : (n*kTOT);
        f_totD(n,i) = sum( f(a,i) );
    end
end

GDP2 = zeros(nTOT,1);
for n = 1:nTOT 
    GDP2(n) = sum (f_totD(n,:));
end

%calculate aggregate final demand TO country A FROM country B (import
%export relationships)
f_totCC = zeros(nTOT,nTOT);
for n = 1:nTOT 
    for j = 1:nTOT
    f_totCC(n,j) = sum(f_totD(n,(n-1)*fTOT+1 : n*fTOT));
    end
end

%% access to data
% support variable for accessing raws: type number of country and sector
% you want the information about; calculate index1 (conceptually access 
% exchanges FROM sector k1 in contry n1 TO index2 or index3)
port = 1;
k1 = 1;

% support variables for accessing columns: type number of country and
% sector you want the information about; calculate index2 (conceptually 
% access exchanges TO sector k2 in contry n2 FROM index1, ex_tr or f_dem)
spain = 1;
k2 = 1;

% support variables for accessing a given type of exogenous flow and a certain
% component of the final demand in country n3 (if necessary)
ex_tr=1;
f_dem= 1; 
n3= 1;

%calculate index to access matrixes
index1 = (port-1)*kTOT + k1;
index2 = (spain-1)*kTOT + k2;
index3 = (n3-1)*fTOT + f_dem;

%access matrixes
A (index1,index2)
L (index1,index2)
E (ex_tr,index1)
R (ex_tr,index1)
v (f_dem,index1)
x0 (index1,index3)
f  (index1,index3)
f_tot (index1)

%% other calculations and support informations
%index values for portugal, trasport sector and households demand
n = 22;
sec = 57;
fd = 1; 


%sup stands for support variable

sup_n = Z((n-1)*kTOT+1:n*kTOT,:); 

j = 1;
for i=1:kTOT*nTOT
    if mod(i-sec,kTOT) == 0
        sup_nsec(:,j) = sup_n(:,i);
        j = j + 1;
    end
end

intersecn = sum(sup_nsec,1); % intermediate demand in country n all sectors for product sec in all countries

%extract the n country fd final demand sector (portf) and select only
%raws related to sec sector (porttr). This way we have a clear 
%rapresentation of who n country buy their sec product from.
supf_n = f(:,(n-1)*fTOT + fd); 

j = 1;
for i=1:kTOT*nTOT
    if mod(i-sec,kTOT) == 0
        finalsecn(j) = supf_n(i); %final demand (from fd final demand sector) for sec product from all countries in country n
        j = j + 1;
    end
end

%% electricity sector analysis

all_sec=1:163;
ee_ind = 96:107;
port = 22; 
spain = 9;
fd = 1:7;

supfinalee = f((port-1)*kTOT + ee_ind, (port-1)*fTOT + fd);
finalee = sum(supfinalee,2)';
supinteree = Z((port-1)*kTOT + all_sec, (port-1)*kTOT + ee_ind);
interee = sum(supinteree,1);
%searching for imports
importee = Z((port-1)*kTOT + ee_ind, (port-1)*kTOT + ee_ind);



%% WORLD I/O
%More interesting! Calculate a world I/O (in order to have an idea of gross
%relationship, regardless of commercial relationship between states)

Zaggr = zeros(kTOT,kTOT);
faggr = zeros(kTOT,fTOT);
Baggr = zeros(trTOT,kTOT);
vaggr = zeros(fTOT,kTOT);

for i=1:kTOT
    for j=1:kTOT
        for port=1:nTOT
            
            if i<trTOT+1
                Baggr(i,j)= Baggr(i,j) + B(i,(port-1)*trTOT + j);
            end
            
            if i<fTOT+1
                vaggr(i,j)= vaggr(i,j) + v(i,(port-1)*fTOT + j);
            end
            
            for spain=1:nTOT
            
                if j<fTOT+1
                    faggr(i,j)= faggr(i,j) + f((port-1)*kTOT + i,(spain-1)*fTOT + j);
                end
                
                Zaggr(i,j)=Zaggr(i,j) + Z((port-1)*kTOT + i,(spain-1)*kTOT + j);
                
            end
        end
    end
end