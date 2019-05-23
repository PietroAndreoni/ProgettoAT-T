%% data
target_ec = 750000; %target for electric car penetration
avcost_ec = ;       %average cost of electric cars (in M€)
ref_ec = ;          %electric cars circulating in reference year

chargdiff = ;       %number of charging station per 1000 vehicles
chargcost = ;       %cost of a single charging station (in M€)
ref_charg = ;       %charging stations operating in reference year

aveff_ec = ;        %average efficiency for electric cars in kW/100km
aveff_nec = ;       %average efficiency for non electric cars in l/100km
avdist = ;          %average distance routed in a year

price_ff = ;        %average substitute fuels price (gasoline+gas) in M€/l
price_el = ;        %average price of electricity in M€/kW

motmet = ;          %iron and hard metals necessary for building an engine
price_ir = ;        %price of iron

litbat = ;          %litium necessary to build a - kW car battery
price_lit = ;       %price of litium

%% calculations

totinv_h = (target_ec - ref_ec) * avcost_ec;
totinv_g = ( (target_ec - ref_ec) / 1000 * chargdiff - ref_charg) * chargcost;
totinv_wind = ;
dem_ff = (target_ec - ref_ec) * avdist /100 * aveff_nec * price_ff;
dem_ee = (target_ec - ref_ec) * avdist /100 * aveff_ec * price_ee;
dem_ir = (target_ec - ref_ec) * motmet * price_ir;
dem_lit = (target_ec - ref_ec) * litbat * price_lit;

%% analysis


%% implementation
%useful indexes: countries
port = 22; 
spain = 9;
france = 11;
germany = 6;
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

% final demand investment increase associated to electric car purchase by
% households (bought from France, Germany and UK)
f1((france-1)*kTOT + tr,(port-1)*fTOT + h)) = f((france-1)*kTOT + tr,(port-1)*fTOT + h)) + 0.33*totinv_h;
f1((germany-1)*kTOT + tr,(port-1)*fTOT + h)) = f((germany-1)*kTOT + tr,(port-1)*fTOT + h)) + 0.33*totinv_h;
f1((UK-1)*kTOT + tr,(port-1)*fTOT + h)) = f((UK-1)*kTOT + tr,(port-1)*fTOT + h)) + 0.33*totinv_h;

%final demand investment by government due to installation of charging
%station (electrical components are assumed to be bought in Portugal)
%cost is assumed to be splitted into 70% costruction, 30% electric
%components
f1((port-1)*kTOT + elman,(port-1)*fTOT + g)) = f((port-1)*kTOT + elman,(port-1)*fTOT + g)) + 0.7*totinv_g;
f1((port-1)*kTOT + costr,(port-1)*fTOT + g)) = f((port-1)*kTOT + costr,(port-1)*fTOT + g)) + 0.3*totinv_g;

% technology change due to transformation of transport sector caused by
% increase in demand for electric cars for these countries
% france import batteries from poland, germany from hungary and UK self
% produces
Z1( (france-1)*kTOT + tr, (poland-1)*kTOT + elman ) = Z( (france-1)*kTOT + tr, (poland-1)*kTOT + elman ) + ; %increase in battery demand from transport sector
Z1( (poland-1)*kTOT + elman, (southam-1)*kTOT + nmet_dir) = Z( (poland-1)*kTOT + elman, (southam-1)*kTOT + nmet_dir) + ; %increase in litium demand from battery sector

Z1( (germany-1)*kTOT + tr, (hungary-1)*kTOT + elman ) = Z( (germany-1)*kTOT + tr, (hungary-1)*kTOT + elman ) + ;
Z1( (hungary-1)*kTOT + elman, (southam-1)*kTOT + nmet_dir) = Z( (hungary-1)*kTOT + elman, (southam-1)*kTOT + nmet_dir) + ;

Z1( (UK-1)*kTOT + tr, (UK-1)*kTOT + elman ) = Z( (UK-1)*kTOT + tr, (UK-1)*kTOT + elman ) + ;
Z1( (UK-1)*kTOT + elman, (southam-1)*kTOT + nmet_dir) = Z( (UK-1)*kTOT + elman, (southam-1)*kTOT + nmet_dir) + ;

%technology change due to charging stations costruction ?

%demand shift from fossil fuel to electricity
%decrease in fossil fuel demand
for n=1:nTOT
f1((n-1)*kTOT + ffref,(port-1)*fTOT + h ) = f((n-1)*kTOT + ffref,(port-1)*fTOT + h ) - dem_ff/nTOT;
end

%% electricity: SCENARIO 1
%increase in electricity demand covered by programmable plants (biomass)
f1((port-1)*kTOT + biomass,(port-1)*fTOT + h ) = f((port-1)*kTOT + biomass,(port-1)*fTOT + h ) + dem_ee;

%% electricity: SCENARIO 2
%increase in total demand is supported by added capacity in wind 
f1((port-1)*kTOT + wind,(port-1)*fTOT + h ) = f((port-1)*kTOT + biomass,(port-1)*fTOT + h ) + dem_ee;

%investment for plants costruction (fix proportions)
f1((port-1)*kTOT + macheq,(port-1)*fTOT + g ) = f((port-1)*kTOT + macheq,(port-1)*fTOT + g ) + 0.7*totinv_wind;
f1((port-1)*kTOT + costr,(port-1)*fTOT + g ) = f((port-1)*kTOT + costr,(port-1)*fTOT + g ) + 0.3*totinv_wind;

%% electricity: SCENARIO 3 
%extra demand is covered with import via Spain

