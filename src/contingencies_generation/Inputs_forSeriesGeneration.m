%%% --- Main user inputs ---

mpc = grid_model;  % Matpower file name (to load the grid)
N = 5; % Number of time series to be sampled (for N-1 independent unavailabilities)
output_folder = ['availability_series']

%%% --- Reliability data ---

% Generator
MTTRgen_ac = 50; % hour
FOR_genac = 0.03; % unavailability rate (forced outage rate)

% Transformer
MTTRtr_ac = 768; % hour
FOR_tr_ac = 0.02; % unavailability rate

% Converter
MTTRconv_dc = 50; % hour
FOR_convdc = 0.05; % unavailability rate

% AC lines
MTTRbr_ac = 16; % hour (for normal events)
MTTRbr_ac_resil = 4*30*24; % hour (for extreme events)
failurerate_fixe_br_ac = 0.22; % per year
failurerate_perkm_br_ac = 0.52*1.6/100; % per year per km

% DC OHL
MTTRbr_ohldc = 30; % hour (for normal events)
MTTRbr_ohldc_resil = 4*30*24; % hour (for extreme events)
failurerate_fixe_ohldc = 0.1; % per year
failurerate_perkm_ohldc = 0.05/100; % per year per km

% DC cable
MTTRcable_dc = 1500; % hour (for normal events)
MTTRcable_dc_resil = 1500; % hour (for extreme events)
failurerate_perkm_cabledc = 0.07/100; % per year per km

% Correlation between pole failures
corr_dcpoles_adequacy = 0.5; % correlation among the wires (P, N and R) of DC branches for adequacy assessment analyses on N-1 ctgs
DCdependent_adequacy = 0; % considering correlation among the wires of DC branches (1= considered, 0=not considered)
corr_dcpoles_resilience = 0.7; % correlation among the wires (P, N and R) of DC branches for resilience assessment analyses (N-k ctgs)

%%% --- Other hypotheses ---

warning off

type_of_unavailabilities = 0; % 0 = only N-1, 1 = only N-k, 2 = both

%%%%%%%% adding data about type of DC branch (1 = overhead, 0 = cable)
mpc.branchdc(:,17) = 1; %
% adding percentage of cable type over the total length (not used in the % release)
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

%%% --- Computation of MTTR & MTTF for each component ---

%% Line lengths & transformer indexes
Llinee = mpc.branch(find(mpc.branch(:,9)==0),4)/0.02;
LeLinee = find(mpc.branch(:,9)==0);
iTrafi = find(mpc.branch(:,9)>0);
LlineeDc = mpc.branchdc(idexdc,3)/0.002;

% Generator
MTTRsgen_ac = ones(1,size(mpc.gen,1)).*MTTRgen_ac;
MTTFsgen_ac = MTTRsgen_ac*(1/FOR_genac - 1);

% Transformer
MTTFbrs_ac(iTrafi) = MTTRtr_ac*(1/FOR_tr_ac - 1);
MTTRbrs_ac(iTrafi) = MTTRtr_ac;

% Converter
MTTRsconv_dc = ones(1,size(mpc.convdc,1)).*MTTRconv_dc;
MTTFsconv_dc = MTTRsconv_dc.*(1/FOR_convdc - 1);

% AC lines
MTTRbrs_ac(LeLinee) = MTTRbr_ac;
MTTRbrs_ac_resil = ones(1,length(MTTFbrs_ac))*MTTRbr_ac_resil;
MTTFbrs_ac(LeLinee) = 8760./(Llinee*failurerate_perkm_br_ac + failurerate_fixe_br_ac);

% DC OHL
MTTRbrs_dc(find(mpc.branchdc(idexdc,17)==1)) = MTTRbr_ohldc; % suppose all are cables
MTTFbrs_dc(find(mpc.branchdc(idexdc,17)==1)) = 8760./(LlineeDc(find(mpc.branchdc(idexdc,17)==1))'.*failurerate_perkm_ohldc + failurerate_fixe_ohldc);
MTTRbrs_dc_resil(find(mpc.branchdc(idexdc,17)==1)) = ones(1,length((find(mpc.branchdc(idexdc,17)==1))))*MTTRbr_ohldc_resil;

% DC cable
MTTRbrs_dc(find(mpc.branchdc(idexdc,17)==0)) = MTTRcable_dc; % suppose all are cables
MTTFbrs_dc(find(mpc.branchdc(idexdc,17)==0)) = 8760./(LlineeDc(find(mpc.branchdc(idexdc,17)==0))'.*failurerate_perkm_cabledc);
MTTRbrs_dc_resil(find(mpc.branchdc(idexdc,17)==0)) = ones(1,length(find(mpc.branchdc(idexdc,17)==0)))*MTTRcable_dc_resil;
