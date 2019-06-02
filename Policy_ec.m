global nTOT kTOT fTOT trTOT vTOT
load("Baseline.mat");

%% data

count = 0;

target_ec = 75000;  %target for electric car penetration
dummy = 1;         % if target=750000(all years) dummy = 10, else dummy = 1
avcost_ec = 0.018;  %average cost of electric cars (in M€)
ref_ec = 0;         %electric cars circulating in reference year

avcost_nec = 0.018;

chargdiff = 6;      %number of charging station per 1000 vehicles (source: https://www.sciencedirect.com/science/article/pii/S1361920917305643)
chargcost = 0.03;   %cost of a single charging station (in M€)
ref_charg = 0;      %charging stations operating in reference year

aveff_ec = 17.9;    %average efficiency for electric cars in kWh/100km
aveff_nec = 7.5;    %average efficiency for non electric cars in l/100km
avdist = 14000;     %average distance routed in a year

price_ff = 0.00000148;   %average substitute fuels price (gasoline+gas) in M€/l
price_ee = 0.00000025;    %average price of electricity in M€/kWh

litbat = 30;          %litium necessary to build a 30 kWh car battery in kg (1kg/kWh)
price_lit = 0.00002;  %price of litium in M€/kg

price_bat = 0.003; %battery price (30 KWh) assuming a target price of 100€/kWh - ie perfect competitivity 

av_ageec = 12;          %average useful life of an electric car 
av_agenec = 12;          %average useful life of non electric car 
av_agecc = 25;          %average useful life of a charging station
av_agew = 20;           %average useful life of a wind farm

totinv_h = (target_ec - ref_ec) * avcost_ec;
totinv_g = ( (target_ec - ref_ec) / 1000 * chargdiff - ref_charg) * chargcost;
totinv_wind = 95.56;
dem_ff = (target_ec - ref_ec) * avdist /100 * aveff_nec * price_ff;
dem_ee = (target_ec - ref_ec) * avdist /100 * aveff_ec * price_ee;
dem_lit = (target_ec - ref_ec) * litbat * price_lit;

%useful indexes: countries
port = 22; 
spain = 9;
france = 11;
germany = 6;
poland = 21;
hungary = 13;
UK = 27;
japan = 29;
china = 30;
southam = 45;

%useful indexes: sectors
ffref = 57;
tr = 91;
ir_dir = 72;
ir_ind = 73;
nmet_dir = 82;
nmet_ind = 83;
macheq = 86;
elman = 88;
costr = 113;

%useful indexes: electricity sector
coal = 96;
gas = 97;
hydro = 99;
wind = 100;
oil = 101;
biomass = 102;
solarth = 103;
solarPV = 104;

%useful indexes: final demand
h = 1;
g = 3;

%matrixes for the policy
f1 = f;
B1 = B;
v1 = v;
Z1 = Z;

%electric auto share import between import countries
ger_share = 0.33;
fr_share = 0.33;
UK_share = 0.33;

%investment share for charging station between transport 
costr_sharecs = 0.6;
elmac_sharecs = 0.4;

%investment share for wind turbines
costr_sharew = 0.1;
macheq_sharew = 0.77;
elmac_sharew = 0.13;


%% UNA TANTUM investments (to be divided for the useful life of the product)
count = count +1;

% final demand investment increase associated to electric car purchase by
% households (bought from France, Germany and UK)
f1((france-1)*kTOT + tr,(port-1)*fTOT + h) = f1((france-1)*kTOT + tr,(port-1)*fTOT + h) + fr_share*totinv_h/av_ageec;
f1((germany-1)*kTOT + tr,(port-1)*fTOT + h) = f1((germany-1)*kTOT + tr,(port-1)*fTOT + h) + ger_share*totinv_h/av_ageec;
f1((UK-1)*kTOT + tr,(port-1)*fTOT + h) = f1((UK-1)*kTOT + tr,(port-1)*fTOT + h) + UK_share*totinv_h/av_ageec;

%final demand decrease associated to less non electric car purchase by
%households (uniformly divided respecting actual import share)
supf_a = f(:,(port-1)*fTOT + h); 

j = 1;
for i=1:kTOT*nTOT
    if mod(i-tr,kTOT) == 0
        supfinaltransa(j) = supf_a(i); 
        j = j + 1;
    end
end

finaltransa = supfinaltransa/sum(supfinaltransa);

for n=1:nTOT
f1((n-1)*kTOT + tr,(port-1)*fTOT + h ) = f1((n-1)*kTOT + tr,(port-1)*fTOT + h ) - finaltransa(n)*avcost_nec*(target_ec-ref_ec)/av_agenec;
end

%final demand investment by government due to installation of charging
%station (electrical compon-ents are assumed to be bought in Portugal)
%cost is assumed to be splitted into 70% costruction, 30% electric
%components
f1((port-1)*kTOT + elman,(port-1)*fTOT + g) = f1((port-1)*kTOT + elman,(port-1)*fTOT + g) + costr_sharecs*totinv_g/av_agecc;
f1((port-1)*kTOT + costr,(port-1)*fTOT + g) = f1((port-1)*kTOT + costr,(port-1)*fTOT + g) + elmac_sharecs*totinv_g/av_agecc;

%% technology changes
count = count +1;
% technology change due to transformation of transport sector caused by
% increase in demand for electric cars for these countries
% france import batteries from poland, germany from hungary and UK self
% produces
Z1( (poland-1)*kTOT + elman, (france-1)*kTOT + tr  ) = Z1( (poland-1)*kTOT + elman, (france-1)*kTOT + tr ) + fr_share*price_bat*target_ec/dummy; %increase in battery demand from transport sector
Z1( (southam-1)*kTOT + nmet_dir, (poland-1)*kTOT + elman ) = Z1( (southam-1)*kTOT + nmet_dir, (poland-1)*kTOT + elman ) + fr_share*dem_lit/dummy; %increase in litium demand from battery sector

Z1( (hungary-1)*kTOT + elman, (germany-1)*kTOT + tr ) = Z1( (hungary-1)*kTOT + elman, (germany-1)*kTOT + tr ) + ger_share*price_bat*target_ec/dummy;
Z1( (southam-1)*kTOT + nmet_dir, (hungary-1)*kTOT + elman) = Z1( (southam-1)*kTOT + nmet_dir, (hungary-1)*kTOT + elman) + ger_share*dem_lit/dummy;

Z1( (UK-1)*kTOT + elman, (UK-1)*kTOT + tr ) = Z1( (UK-1)*kTOT + elman, (UK-1)*kTOT + tr ) + UK_share*price_bat*target_ec/dummy;
Z1( (southam-1)*kTOT + nmet_dir, (UK-1)*kTOT + elman) = Z1( (southam-1)*kTOT + nmet_dir, (UK-1)*kTOT + elman) + UK_share*dem_lit/dummy;

%% annual demand shock (not to be divided for useful life)
count = count +1; 
%demand shift from fossil fuel to electricity
%decrease in fossil fuel demand
supf_p = f(:,(port-1)*fTOT + h); 

j = 1;
for i=1:kTOT*nTOT
    if mod(i-ffref,kTOT) == 0
        supfinalrefp(j) = supf_p(i); %final demand (from fd final demand sector) for sec product from all countries in country n
        j = j + 1;
    end
end

finalrefp = supfinalrefp/sum(supfinalrefp);

for n=1:nTOT
f1((n-1)*kTOT + ffref,(port-1)*fTOT + h ) = f1((n-1)*kTOT + ffref,(port-1)*fTOT + h ) - finalrefp(n)*dem_ff;
end

%% electricity: SCENARIO 1
count = count +1;
%increase in electricity demand covered by programmable plants (biomass)
f1((port-1)*kTOT + biomass,(port-1)*fTOT + h ) = f1((port-1)*kTOT + biomass,(port-1)*fTOT + h ) + dem_ee;

%% electricity: SCENARIO 2
count = count +1;
%increase in total demand is supported by added capacity in wind 
f1((port-1)*kTOT + wind,(port-1)*fTOT + h ) = f1((port-1)*kTOT + wind,(port-1)*fTOT + h ) + dem_ee;

%investment for plants costruction (fix proportions)
f1((port-1)*kTOT + macheq,(port-1)*fTOT + g ) = f1((port-1)*kTOT + macheq,(port-1)*fTOT + g ) + macheq_sharew*totinv_wind/av_agew;
f1((port-1)*kTOT + costr,(port-1)*fTOT + g ) = f1((port-1)*kTOT + costr,(port-1)*fTOT + g ) + costr_sharew*totinv_wind/av_agew;
f1((port-1)*kTOT + elman,(port-1)*fTOT + g ) = f1((port-1)*kTOT + elman,(port-1)*fTOT + g ) + elmac_sharew*totinv_wind/av_agew;

%% electricity: SCENARIO 3 
count = count +1;
%extra demand is covered omogenously by the existing mix
ttotal_prod = sum( f1( (port-1)*kTOT + (96:107), (port-1)*fTOT + h ) ) ; 
share_prod = f1( (port-1)*kTOT + (96:107), (port-1)*fTOT + h ) ./ ttotal_prod;

for i=96:107
    f1((port-1)*kTOT + i,(port-1)*fTOT + h ) = f1((port-1)*kTOT + i,(port-1)*fTOT + h ) + share_prod(i-95)*dem_ee;
end

%% Result calculations
xnon = diag(x+0.0001);
A1 = Z1 / xnon;
Lnon1 = diag(ones(size(A1,1),1))-A1;
L1 = inv(diag(ones(size(A1,1),1))-A1);
x1 = Lnon1\f1;

f_tot1=zeros(nTOT*kTOT,1);
for i = 1:nTOT*kTOT
    f_tot1(i) = sum( f1(i,:) );
end

R1 = B1* ( diag(Lnon1\f_tot1) );
E1 = (B1/Lnon1) * diag(f_tot1);

DZ = Z1 - Z;
Df = f1 - f;
DB = B1 - B;
DA = A1 - A;
DR = R1 - R;
DE = E1 - E;
Dv = v1 - v;

%% Result analysis

info = delta_analysis(f,1000,1,1);
[DZaggr_sec,DAaggr_sec,Dfaggr_sec,DBaggr_sec,DRaggr_sec,DEaggr_sec,Dvaggr_sec] = aggrbysec (DZ,DA,Df,DB,DR,DE,Dv);
[DZaggr_co,DAaggr_co,Dfaggr_co,DBaggr_co,DRaggr_co,DEaggr_co,Dvaggr_co] = aggrbycountry (DZ,DA,Df,DB,DR,DE,Dv);

[Zaggr_co,Aaggr_co,faggr_co,Baggr_co,Raggr_co,Eaggr_co,vaggr_co] = aggrbycountry (Z1,A1,f1,B1,R1,E1,v1);
[Zaggr_sec,Aaggr_sec,faggr_sec,Baggr_sec,Raggr_sec,Eaggr_sec,vaggr_sec] = aggrbysec (Z1,A1,f1,B1,R1,E1,v1);

fsupp = zeros(7824,366);
for i=1:size(Df,1)
    for j=1:size(Df,2)
        
        if Df(i,j) < 0 
            fsupp(i,j) = 1;
        end
    end
end

xlswrite("resultbio.xlsx",DA,'DA');
xlswrite("resultbio.xlsx",DAaggr_co,'DAco');
xlswrite("resultbio.xlsx",DAaggr_sec,'DAsec');

xlswrite("BIGFIXwind.xlsx",DE,'DE');
xlswrite("BIGFIXwind.xlsx",DEaggr_co,'DEco');
xlswrite("BIGFIXwind.xlsx",DEaggr_sec,'DEsec');

xlswrite("BIGFIXwind.xlsx",DR,'DR');
xlswrite("BIGFIXwind.xlsx",DRaggr_co,'DRco');
xlswrite("BIGFIXwind.xlsx",DRaggr_sec,'DRsec');

