%% INPUT SECTION (to be customised depending on the data available for the analysis)
%%%%%%%% MCS for N-1 outages
disp('**** INPUTTING DATA***')
warning off
N = 5; % nr of MCS years to be sampled for N-1 independent unavailabilities
type_of_unavailabilities = 0; % 0 = only N-1, 1 = only N-k, 2 = both
%%%%%%%% plotting options
figu_plot = 0; % 0=no figures, 1 = yes figures
%%%%%%%% input from WEATHER data (here are series of wind speed 10 min average)
if type_of_unavailabilities > 0
infox = ncinfo('fg_ens_mean_0.1deg_reg_1980-1994_v26.0e.nc');
Wind =  ncread('fg_ens_mean_0.1deg_reg_1980-1994_v26.0e.nc','fg');
tempo =  ncread('fg_ens_mean_0.1deg_reg_1980-1994_v26.0e.nc','time');
LONGI =  ncread('fg_ens_mean_0.1deg_reg_1980-1994_v26.0e.nc','longitude');
LATI =  ncread('fg_ens_mean_0.1deg_reg_1980-1994_v26.0e.nc','latitude');
end
% add path to botev tilting algorithm for copula computation
metho =2; % use block maxima, 1 = peak over threshold method to produce GEV distributions
addpath('BOTEV');
%%%%%%%% add data about the GRID AC/DC
% sample of the grid in hvdc-wise format
caso ='grid39'; %grid5
% loading the grid
switch caso
    case 'grid5'
        mpc = case5_2grids_MC;
    case 'grid39'
        mpc = case39_mcdc;
end
%%%%%%%% input TOPOLOGICAL DATA (bus coordinates) for the resilience based
% availability assessment
% in case of missing geo data use the force method for the AC grid
%
if exist('default_coo.mat') > 0
    load default_coo.mat
else

 [ ok TOLLE errorep LUNGH0 distanzaEstim mpc0] = Modello_forze_matp( mpc,1 );

 S.mpc0.bus_coord = mpc0.coord;
% in case of missing geo data use the force method for the DC grid
%
 [ ok TOLLE errorep LUNGH0 distanzaEstim mpc0] = Modello_forze_matp( mpc,2 );

 S.mpc0.busdc_coord = mpc0.busdc_coord;
 save default_coo.mat S
end
mpc.bus_coord=S.mpc0.bus_coord;
mpc.busdc_coord=S.mpc0.busdc_coord;

% consistency check for coordinates
quali_trafi = find(mpc.branch(:,9)>0);

for itr = 1:size(mpc.convdc,1)
    mpc.busdc_coord(mpc.convdc(itr,1),:) = mpc.bus_coord(mpc.convdc(itr,2),:);
end

%%%%%%%% adding data about type of DC branch (1 = overhead, 0= cable)
mpc.branchdc(:,17) = 1; %
% adding percentage of cable type over the total length (not used in the
% release)
mpc.branchdc(:,18) = 0.1; % pu of cable;

% giving names to buses AC and DC
for j = 1:size(mpc.bus,1)
   mpc.bus_name{j} = ['AC-' num2str(mpc.bus(j,1))];
end
for j = 1:size(mpc.busdc,1)
   mpc.busdc_name{j} = ['DC-' num2str(mpc.busdc(j,1))];
end

mpc0 = mpc;
% correspondence between wires N,R,P and DC branches
idexdc = [];poli=[];
for i = 1:size(mpc0.branchdc,1)
    if mpc0.branchdc(i,10)==2
        for ipol = 1:3
            idexdc = [idexdc i];
             poli=[poli ipol-1];
        end
    else
        quali_poli = find(mpc.branchdc(i,14:16)==1); %1 =p,2=n, 3=r
        for ipol = 1:length(quali_poli)
            idexdc = [idexdc i];
            switch quali_poli(ipol)
                case 1
                    poli=[poli 1];
                case 2
                    poli=[poli 2];
                case 3
                    poli=[poli 0];
            end
        end
    end
end

%%%%%%%% definition of the real geo coordinates setting the ZONE UTM where the
% grid is located (if coordinates are available, skip this part)


if type_of_unavailabilities > 0

ZONA = '30U';
[latlim longlim] = utmzone(ZONA);

closesLO1 = find(abs(LONGI-longlim(1))==min(abs(LONGI-longlim(1))));
closesLO2 = find(abs(LONGI-longlim(2))==min(abs(LONGI-longlim(2))));
closesLA1 = find(abs(LATI-latlim(1))==min(abs(LATI-latlim(1))));
closesLA2 = find(abs(LATI-latlim(2))==min(abs(LATI-latlim(2))));


LATI_ = [closesLA1:closesLA2];
LONGI_ = [closesLO1:closesLO2];
LATIc_ = LATI(LATI_);%[closesLA1:closesLA2];
LONGIc_ = LONGI(LONGI_);%[closesLO1:closesLO2];

utmstruct = defaultm('utm');
dczone = utmzone(mean(latlim),mean(longlim));
utmstruct.zone = dczone;
utmstruct.geoid = wgs84Ellipsoid;
utmstruct = defaultm(utmstruct);

%
[x11 y11] = minvtran(utmstruct,5e5+[mpc0.bus_coord(:,1)], 5e6+[mpc0.bus_coord(:,2)]);

mpc0.bus_coord2 = [x11 y11];
[x11d y11d] = minvtran(utmstruct,5e5+[mpc0.busdc_coord(:,1)], 5e6+[mpc0.busdc_coord(:,2)]);

mpc0.busdc_coord2 = [x11d y11d];
end
%%%%%%%% ended coordinates identification based on UTM zone

%%%%%%%% GRID ASSET VULNERABILITY characterization
meanVu = 19; % mean value of resistance to wind speed in m/s
sigmaVu = 0.2; % dev stad for the lognormal fragility curves to wind speed
%%%%%%%% PARAMETers for DEFINITION for line clustering  and ctg identification
% module
SOGLIA_M = 6;% definition of min wind speed to define event matrix M
epsi=0; % min. probability threshold for discarding the AND of trippings
SOGLIA_LI = 5; % max nr of lines for each cluster: suggested a number included between 3 and 8.
pesoSil = 0.5; % weight for the slouette metrics importance in selecting the parameter patterns for clustering technique
min_corr_int_step1 = [0.3:0.1:0.7]; % min. correlation intra cluster at step 1 of clustering technique
max_dist_intercluster_step1 = [0.3:0.1:0.7]; % max distance between clusters at step 1
corr_residua_interclusters_step3 = 0.3; % residual correlation between clusters at step 3
nmax_it_step3 = 100; % max nr of iterations at step 3
avg_corr_min_interclustertopo_step3 = 0.3; % average min correlation between topological clusters in step 3
correlazio_max_interclustertopo_step3 = 0.6; % max correlation between topological clusters a step 3
peso_topologia = 0.5;
min_corr_media_cuts_topologici = 0.2; % min average correlation in toppological cutsets
kmax = SOGLIA_LI; % maximum order for N-k contingencies

%%%%%%%% RELIABILITY DATA FOR DC AND AC GRIDS
%%% for N-1 and DC grid ctgs and adequacy analysis
Llinee = mpc.branch(find(mpc.branch(:,9)==0),4)/0.02;
LeLinee = find(mpc.branch(:,9)==0);
iTrafi = find(mpc.branch(:,9)>0);
LlineeDc = mpc.branchdc(idexdc,3)/0.002;

% AC lines
MTTRbr_ac = 16; % hour for ac lines
MTTRbr_ac_resil = 4*30*24; % hour for ac lines for resiliencey analysis
MTTRtr_ac = 768; % hour for ac lines
MTTRbr_ohldc = 30;
MTTRcable_dc = 1500; % ore
MTTRgen_ac = 50;
MTTRconv_dc = 50;
MTTRbr_ohldc_resil = 4*30*24; % hour for ac lines for resiliencey analysis
MTTRcable_dc_resil = 1500; % hour for ac lines for resiliencey analysis

FOR_tr_ac = 0.02; % per yr
FOR_genac = 0.03;
FOR_convdc = 0.05;
%%%

MTTFbrs_ac(LeLinee) = 8760./(Llinee*0.52*1.6/100 + 0.22);
MTTFbrs_ac(iTrafi) = MTTRtr_ac*(1/FOR_tr_ac - 1);
MTTRbrs_ac(LeLinee) = MTTRbr_ac;
MTTRbrs_ac(iTrafi) = MTTRtr_ac;

%%%% GEN AC
MTTRsgen_ac = ones(1,size(mpc.gen,1)).*MTTRgen_ac;
MTTFsgen_ac = MTTRsgen_ac*(1/FOR_genac - 1);

%DC CONVERTER
MTTRsconv_dc = ones(1,size(mpc.convdc,1)).*MTTRconv_dc;
MTTFsconv_dc = MTTRsconv_dc.*(1/FOR_convdc - 1);

% DC BRANCH
MTTFbrs_dc(find(mpc.branchdc(idexdc,17)==0)) = 8760./(LlineeDc(find(mpc.branchdc(idexdc,17)==0))'.*0.07/100);
MTTRbrs_dc(find(mpc.branchdc(idexdc,17)==0)) = MTTRcable_dc; % suppose all are cables
MTTRbrs_dc(find(mpc.branchdc(idexdc,17)==1)) = MTTRbr_ohldc; % suppose all are cables
MTTFbrs_dc(find(mpc.branchdc(idexdc,17)==1)) = 8760./(LlineeDc(find(mpc.branchdc(idexdc,17)==1))'.*0.5/100 + 0.1);

corr_dcpoles_adequacy = 0.5; % correlation among the wires (P, N and R) of DC branches for adequacy assessment analyses on N-1 ctgs
DCdependent_adequacy = 0; % considering correlation among the wires of DC branches (1= considered, 0=not considered)

%%% FOR RESILIENCE ANALYSES
MTTRbrs_ac_resil = ones(1,length(MTTFbrs_ac))*MTTRbr_ac_resil;
MTTRbrs_dc_resil(find(mpc.branchdc(idexdc,17)==1)) = ones(1,length((find(mpc.branchdc(idexdc,17)==1))))*MTTRbr_ohldc_resil;
MTTRbrs_dc_resil(find(mpc.branchdc(idexdc,17)==0)) = ones(1,length(find(mpc.branchdc(idexdc,17)==0)))*MTTRcable_dc_resil;

corr_dcpoles_resilience = 0.7; % correlation among the wires (P, N and R) of DC branches for resilience assessment analyses (N-k ctgs)


%%%% END DATA INPUT SECTION
