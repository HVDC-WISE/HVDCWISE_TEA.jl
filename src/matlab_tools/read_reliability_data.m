function reliability_data = read_reliability_data(work_dir)
    % addpath(fileparts(mfilename('fullpath')));
    pkg load io

    assert(exist('work_dir', 'var') == 1, "work_dir is not defined")

    simulation_dir = [work_dir, '\simulation_interface'];
    macro_scenario = NaN;
    folders=dir(simulation_dir);
    folders(1:2)=[];
    for i=1:length(folders)
        fifo_name = folders(i).name;
        if fifo_name(end-1:end) == '.m'
            macro_scenario = fifo_name(1:end-2);
        endif
    endfor
    assert(isnan(macro_scenario) == 0, strcat("No .m file has been found in ", simulation_dir))

    %%% --- Main user inputs ---

    %% load grid model (.m file)
    addpath(simulation_dir);
    mpc = feval(macro_scenario);
    rmpath(simulation_dir);

    % mpc = grid_model;  % Matpower file name (to load the grid)
    output_folder = [simulation_dir, '\Input_series\Availability'];

    %%% --- Reliability data ---
    reliability_data_path = [work_dir, '\user_interface\inputs\reliability_data.xlsx'];
    assert(isfile(reliability_data_path) == 1, strcat(reliability_data_path, " does not exist"))

    [~, ~, sheet] = xlsread(reliability_data_path);

    column_titles = {"Parameter", "Value", "Unit", "Scope", "Meaning"};
    for j = 1:5
        assert(sheet{1,j} == column_titles{j}, strcat("First cell of column ", num2str(j), " should be ", column_titles{j}, ", not ", sheet{1,j}))
    end

    Parameter = sheet(2:end, 1);
    Values = cell2mat(sheet(2:end, 2)); % Convert column 1 into a matrix of floats

    parameters = {"N",
    "MTTR_genac", % hour
    "FOR_genac", % unavailability rate (forced outage rate)
    "MTTR_trac",
    "FOR_trac",
    "MTTR_convdc",
    "FOR_convdc",
    "MTTR_ohlac", % hour (for normal events)
    "MTTR_ohlac_resil", % hour (for extreme events)
    "failurerate_fixe_ohlac", % per year
    "failurerate_perkm_ohlac", % per year per km
    "MTTRbr_ohldc",
    "MTTRbr_ohldc_resil",
    "failurerate_fixe_ohldc",
    "failurerate_perkm_ohldc",
    "MTTR_cabledc",
    "MTTR_cabledc_resil",
    "failurerate_perkm_cabledc",
    };

    n_series = Values(1); % Number of time series to be sampled (for N-1 independent unavailabilities)
    [MTTR_genac, FOR_genac] = deal(num2cell(Values(2:3)){:}); % Generator
    [MTTR_trac, FOR_trac] = deal(num2cell(Values(4:5)){:}); % Transformer
    [MTTR_convdc, FOR_convdc] = deal(num2cell(Values(6:7)){:}); % Converter
    [MTTR_ohlac, MTTR_ohlac_resil, failurerate_fixe_ohlac, failurerate_perkm_ohlac] = deal(num2cell(Values(8:11)){:}); % AC lines
    [MTTRbr_ohldc, MTTRbr_ohldc_resil, failurerate_fixe_ohldc, failurerate_perkm_ohldc] = deal(num2cell(Values(12:15)){:}); % DC OHL
    [MTTR_cabledc, MTTR_cabledc_resil, failurerate_perkm_cabledc] = deal(num2cell(Values(16:18)){:}); % DC cable

    % Correlation between pole failures
    corr_dcpoles_adequacy = 0.5; % correlation among the wires (P, N and R) of DC branches for adequacy assessment analyses on N-1 ctgs
    DCdependent_adequacy = 0; % considering correlation among the wires of DC branches (1= considered, 0=not considered)
    corr_dcpoles_resilience = 0.7; % correlation among the wires (P, N and R) of DC branches for resilience assessment analyses (N-k ctgs)

    %%% --- Other hypotheses ---

    warning off

    %%%%%%%% adding data about type of DC branch (1 = overhead, 0 = cable)
    mpc.branchdc(:,17) = 0; %
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
    MTTRsgen_ac = ones(1,size(mpc.gen,1)).*MTTR_genac;
    MTTFsgen_ac = MTTRsgen_ac*(1/FOR_genac - 1);

    % Transformer
    MTTFbrs_ac(iTrafi) = MTTR_trac*(1/FOR_trac - 1);
    MTTRbrs_ac(iTrafi) = MTTR_trac;

    % Converter
    MTTRsconv_dc = ones(1,size(mpc.convdc,1)).*MTTR_convdc;
    MTTFsconv_dc = MTTRsconv_dc.*(1/FOR_convdc - 1);

    % AC lines
    MTTRbrs_ac(LeLinee) = MTTR_ohlac;
    MTTRbrs_ac_resil = ones(1,length(MTTFbrs_ac))*MTTR_ohlac_resil;
    MTTFbrs_ac(LeLinee) = 8760./(Llinee*failurerate_perkm_ohlac + failurerate_fixe_ohlac);

    % DC OHL
    MTTRbrs_dc(find(mpc.branchdc(idexdc,17)==1)) = MTTRbr_ohldc; % suppose all are cables
    MTTFbrs_dc(find(mpc.branchdc(idexdc,17)==1)) = 8760./(LlineeDc(find(mpc.branchdc(idexdc,17)==1))'.*failurerate_perkm_ohldc + failurerate_fixe_ohldc);
    MTTRbrs_dc_resil(find(mpc.branchdc(idexdc,17)==1)) = ones(1,length((find(mpc.branchdc(idexdc,17)==1))))*MTTRbr_ohldc_resil;

    % DC cable
    MTTRbrs_dc(find(mpc.branchdc(idexdc,17)==0)) = MTTR_cabledc; % suppose all are cables
    MTTFbrs_dc(find(mpc.branchdc(idexdc,17)==0)) = 8760./(LlineeDc(find(mpc.branchdc(idexdc,17)==0))'.*failurerate_perkm_cabledc);
    MTTRbrs_dc_resil(find(mpc.branchdc(idexdc,17)==0)) = ones(1,length(find(mpc.branchdc(idexdc,17)==0)))*MTTR_cabledc_resil;

    % Data to return
    reliability_data.n_series = n_series;
    reliability_data.mpc = mpc;
    # reliability_data.mpc0 = mpc0;
    reliability_data.output_folder = output_folder;

    # reliability_data.MTTR_genac = MTTR_genac;
    # reliability_data.FOR_genac = FOR_genac;
    # reliability_data.MTTR_trac = MTTR_trac;
    # reliability_data.FOR_trac = FOR_trac;
    # reliability_data.MTTR_convdc = MTTR_convdc;
    # reliability_data.FOR_convdc = FOR_convdc;
    # reliability_data.MTTR_ohlac = MTTR_ohlac;
    # reliability_data.MTTR_ohlac_resil = MTTR_ohlac_resil;
    # reliability_data.failurerate_fixe_ohlac = failurerate_fixe_ohlac;
    # reliability_data.failurerate_perkm_ohlac = failurerate_perkm_ohlac;
    # reliability_data.MTTRbr_ohldc = MTTRbr_ohldc;
    # reliability_data.MTTRbr_ohldc_resil = MTTRbr_ohldc_resil;
    # reliability_data.failurerate_fixe_ohldc = failurerate_fixe_ohldc;
    # reliability_data.failurerate_perkm_ohldc = failurerate_perkm_ohldc;
    # reliability_data.MTTR_cabledc = MTTR_cabledc;
    # reliability_data.MTTR_cabledc_resil = MTTR_cabledc_resil;
    # reliability_data.failurerate_perkm_cabledc = failurerate_perkm_cabledc;

    reliability_data.MTTRsgen_ac = MTTRsgen_ac;
    reliability_data.MTTFsgen_ac = MTTFsgen_ac;

    reliability_data.MTTFbrs_ac = MTTFbrs_ac;
    reliability_data.MTTRbrs_ac = MTTRbrs_ac;

    reliability_data.MTTRbrs_dc = MTTRbrs_dc;
    reliability_data.MTTFbrs_dc = MTTFbrs_dc;

    reliability_data.MTTRsconv_dc = MTTRsconv_dc;
    reliability_data.MTTFsconv_dc = MTTFsconv_dc;

    reliability_data.corr_dcpoles_adequacy = corr_dcpoles_adequacy;
    reliability_data.DCdependent_adequacy = DCdependent_adequacy;
    reliability_data.corr_dcpoles_resilience = corr_dcpoles_resilience;

    reliability_data.idexdc = idexdc;
    reliability_data.poli = poli;
end
