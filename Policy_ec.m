%% data
global nTOT kTOT fTOT trTOT vTOT

target_ec = 75000;  %target for electric car penetration
avcost_ec = 0.04;   %average cost of electric cars (in M€)
ref_ec = 0;         %electric cars circulating in reference year

avcost_nec = 0.015;

chargdiff = 6;      %number of charging station per 1000 vehicles (source: https://www.sciencedirect.com/science/article/pii/S1361920917305643)
chargcost = 0.03;   %cost of a single charging station (in M€)
ref_charg = 0;      %charging stations operating in reference year

aveff_ec = 17.9;    %average efficiency for electric cars in kW/100km
aveff_nec = 7.5;    %average efficiency for non electric cars in l/100km
avdist = 14000;     %average distance routed in a year

price_ff = 0.00000148;   %average substitute fuels price (gasoline+gas) in M€/l
price_ee = 0.000003;     %average price of electricity in M€/kW

motmet = 0;          %iron and hard metals necessary for building an engine
price_ir = 0;        %price of iron

litbat = 14;          %litium necessary to build a - kW car battery in kg
price_lit = 0.00002;  %price of litium in M€/kg

price_bat = 0.006227; %battery price    

%% calculations

totinv_h = (target_ec - ref_ec) * avcost_ec;
totinv_g = ( (target_ec - ref_ec) / 1000 * chargdiff - ref_charg) * chargcost;
totinv_wind = 0;
dem_ff = (target_ec - ref_ec) * avdist /100 * aveff_nec * price_ff;
dem_ee = (target_ec - ref_ec) * avdist /100 * aveff_ec * price_ee;
dem_ir = (target_ec - ref_ec) * motmet * price_ir;
dem_lit = (target_ec - ref_ec) * litbat * price_lit;


%% implementation
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
A1 = A;
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
costr_sharew = 0.7;
macheq_sharew = 0.3;

% final demand investment increase associated to electric car purchase by
% households (bought from France, Germany and UK)
f1((france-1)*kTOT + tr,(port-1)*fTOT + h) = f((france-1)*kTOT + tr,(port-1)*fTOT + h) + fr_share*totinv_h;
f1((germany-1)*kTOT + tr,(port-1)*fTOT + h) = f((germany-1)*kTOT + tr,(port-1)*fTOT + h) + ger_share*totinv_h;
f1((UK-1)*kTOT + tr,(port-1)*fTOT + h) = f((UK-1)*kTOT + tr,(port-1)*fTOT + h) + UK_share*totinv_h;

%final demand decrease associated to less electric car purchase by
%households (uniformly divided respecting actual import share)
supf_a = f(:,(port-1)*fTOT + h); 

j = 1;
for i=1:kTOT*nTOT
    if mod(i-tr,kTOT) == 0
        supfinaltransa(j) = supf_a(i); %final demand (from fd final demand sector) for sec product from all countries in country n
        j = j + 1;
    end
end

finaltransa = supfinaltransa/sum(supfinaltransa);

for n=1:nTOT
f1((n-1)*kTOT + tr,(port-1)*fTOT + h ) = f((n-1)*kTOT + tr,(port-1)*fTOT + h ) - finaltransa(n)*avcost_nec*target_ec;
end



%final demand investment by government due to installation of charging
%station (electrical components are assumed to be bought in Portugal)
%cost is assumed to be splitted into 70% costruction, 30% electric
%components
f1((port-1)*kTOT + elman,(port-1)*fTOT + g) = f((port-1)*kTOT + elman,(port-1)*fTOT + g) + costr_sharecs*totinv_g;
f1((port-1)*kTOT + costr,(port-1)*fTOT + g) = f((port-1)*kTOT + costr,(port-1)*fTOT + g) + elmac_sharecs*totinv_g;

% technology change due to transformation of transport sector caused by
% increase in demand for electric cars for these countries
% france import batteries from poland, germany from hungary and UK self
% produces
Z1( (france-1)*kTOT + tr, (poland-1)*kTOT + elman ) = Z( (france-1)*kTOT + tr, (poland-1)*kTOT + elman ) + fr_share*price_bat; %increase in battery demand from transport sector
Z1( (poland-1)*kTOT + elman, (southam-1)*kTOT + nmet_dir) = Z( (poland-1)*kTOT + elman, (southam-1)*kTOT + nmet_dir) + fr_share*dem_lit; %increase in litium demand from battery sector

Z1( (germany-1)*kTOT + tr, (hungary-1)*kTOT + elman ) = Z( (germany-1)*kTOT + tr, (hungary-1)*kTOT + elman ) + ger_share*price_bat;
Z1( (hungary-1)*kTOT + elman, (southam-1)*kTOT + nmet_dir) = Z( (hungary-1)*kTOT + elman, (southam-1)*kTOT + nmet_dir) + ger_share*dem_lit;

Z1( (UK-1)*kTOT + tr, (UK-1)*kTOT + elman ) = Z( (UK-1)*kTOT + tr, (UK-1)*kTOT + elman ) + UK_share*price_bat;
Z1( (UK-1)*kTOT + elman, (southam-1)*kTOT + nmet_dir) = Z( (UK-1)*kTOT + elman, (southam-1)*kTOT + nmet_dir) + UK_share*dem_lit;

%technology change due to charging stations costruction ?

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
f1((n-1)*kTOT + ffref,(port-1)*fTOT + h ) = f((n-1)*kTOT + ffref,(port-1)*fTOT + h ) - finalrefp(n)*dem_ff;
end

%% electricity: SCENARIO 1
%increase in electricity demand covered by programmable plants (biomass)
f1((port-1)*kTOT + biomass,(port-1)*fTOT + h ) = f((port-1)*kTOT + biomass,(port-1)*fTOT + h ) + dem_ee;

%% electricity: SCENARIO 2
%increase in total demand is supported by added capacity in wind 
f1((port-1)*kTOT + wind,(port-1)*fTOT + h ) = f((port-1)*kTOT + biomass,(port-1)*fTOT + h ) + dem_ee;

%investment for plants costruction (fix proportions)
f1((port-1)*kTOT + macheq,(port-1)*fTOT + g ) = f((port-1)*kTOT + macheq,(port-1)*fTOT + g ) + macheq_sharew*totinv_wind;
f1((port-1)*kTOT + costr,(port-1)*fTOT + g ) = f((port-1)*kTOT + costr,(port-1)*fTOT + g ) + costr_sharew*totinv_wind;

%% electricity: SCENARIO 3 
%extra demand is covered with import via Spain

%% Result calculations

A1 = Z1 * inv(diag(x+0.0001));
B1 = B;
L1 = inv(diag(ones(size(A,1),1))-A1);
x1 = L*f1;

f_tot1=zeros(nTOT*kTOT,1);
for i = 1:nTOT*kTOT
    f_tot1(i) = sum( f1(i,:) );
end

R1 = B1* ( diag(L1*f_tot1) );
E1 = (B1*L1) * diag(f_tot1);

DA = A - A1;
DR = B - R1;
DE = E - E1;

%% Result analysis
info = delta_analysis(DR,1000,0,1);

