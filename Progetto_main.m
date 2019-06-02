%% Basic calculations
clear ;
load("base_students.mat");

% number of countries, sectors, final demand types and exogenous resources
% in the model
global nTOT kTOT fTOT trTOT vTOT

nTOT = 48;
kTOT = 163;
fTOT = 7;
trTOT = 12;
vTOT = 7;

%%
port = 22; 

total_prod = sum( f( (port-1)*kTOT + (96:107), (port-1)*fTOT + 1 )); 
share_new = [0.23, 0.125, 0, 0.3, 0.233, 0.032, 0.064, 0.012, 0, 0, 0.004, 0];

for el = 96:107
    f( (port-1)*kTOT + el, (port-1)*fTOT + 1 ) =  share_new(el-95)*total_prod;
end

%wip
a = sum ( Z ( :, (UK-1)*fTOT + (96:107) ), 1);
tot = sum (a);
share = a/tot;

a = sum ( sum( Z,1));
b = sum ( sum( Z,2));

clear el;
%%
%calculate Leontieff coefficients matrix and overall production
L = inv(diag(ones(size(A,1),1))-A);
Lnon = diag(ones(size(A,1),1))-A;
x0 = Lnon\f;
x = sum(x0,2);
Z = A*diag(x+0.0001);

%calculate aggregate final demand for sector i in country A
f_tot=zeros(nTOT*kTOT,1);
for i = 1:nTOT*kTOT
    f_tot(i) = sum( f(i,:) );
end

%calculate exogenous transaction matrix and embodied resource matrix
R = B* ( diag(Lnon\f_tot) );
E = (B/Lnon) * diag(f_tot);
 

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
sec = 91;
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

intersecn = sum(sup_nsec,1); % intermediate demand in country n all sectors for product sec by country

%extract the n country fd final demand sector (portf) and select only
%raws related to sec sector (porttr). This way we have a clear 
%rapresentation of who n country buy their sec product from.
supf_n = f(:,(n-1)*fTOT + fd); 

j = 1;
for i=1:kTOT*nTOT
    if mod(i-sec,kTOT) == 0
        finalsecn(j) = supf_n(i); %final demand in country n (from fd final demand sector) for product sec by country
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
importee = Z((port-1)*kTOT + ee_ind, (spain-1)*kTOT + ee_ind);
