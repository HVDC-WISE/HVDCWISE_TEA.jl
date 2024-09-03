clear
clc

pkg load io

inputfile='grid_model'; # grid_model.m has to be in the same folder
simulation_results_dir = fileread('simulation_results_path.txt');
parent_names = strsplit(simulation_results_dir, '\')
macro_scenario = parent_names(end)
work_dir = strjoin(parent_names(1:end-2), '\')
user_results_dir = [work_dir, '\user_interface\results']

%% load simulation input/output data
mpc=eval(inputfile);

folders=dir(simulation_results_dir);
folders(1:2)=[];

for k=1:length(folders)
    folder_path = [folders(k).folder, '\', folders(k).name, '\'];

    KPI.scenario_id{k}=folders(k).name;

    # ps_storage{k,1}=readmatrix([folders(k).folder,'\',folders(k).name,'\storage\ps.csv'],'NumHeaderLines',1,'Delimiter',',')*mpc.baseMVA;
    ps_storage{k,1} = csvread([folder_path, 'storage\ps.csv'], 1, 0)*mpc.baseMVA;  # 1 and 0 means starting from row 2 and col 1

    # pflex_load{k,1}=readmatrix([folders(k).folder,'\',folders(k).name,'\load\pflex.csv'],'NumHeaderLines',1,'Delimiter',',')*mpc.baseMVA;
    pflex_load{k,1} = csvread([folder_path, 'load\pflex.csv'], 1, 0)*mpc.baseMVA;

    # pg_gen{k,1}=readmatrix([folders(k).folder,'\',folders(k).name,'\gen\pg.csv'],'NumHeaderLines',1,'Delimiter',',')*mpc.baseMVA;
    pg_gen{k,1} = csvread([folder_path, 'gen\pg.csv'], 1, 0)*mpc.baseMVA;

    # pg_curt{k,1}=readmatrix([folders(k).folder,'\',folders(k).name,'\gen\pgcurt.csv'],'NumHeaderLines',1,'Delimiter',',')*mpc.baseMVA;
    pg_curt{k,1} = csvread([folder_path, 'gen\pgcurt.csv'], 1, 0)*mpc.baseMVA;

    # pconv_dc{k,1}=readmatrix([folders(k).folder,'\',folders(k).name,'\convdc\pconv.csv'],'NumHeaderLines',1,'Delimiter',',')*mpc.baseMVA;
    % iconv_dc{k,1}=readmatrix([folders(k).folder,'\',folders(k).name,'\convdc\iconv_dc.csv'],'NumHeaderLines',1,'Delimiter',',')*mpc.baseMVA;
    pconv_dc{k,1} = csvread([folder_path, 'convdc\pconv.csv'], 1, 0)*mpc.baseMVA;


    # p_slack{k,1}=readmatrix([folders(k).folder,'\',folders(k).name,'\bus\p_slack_up.csv'],'NumHeaderLines',1,'Delimiter',',')*mpc.baseMVA-...
    #     readmatrix([folders(k).folder,'\',folders(k).name,'\bus\p_slack_down.csv'],'NumHeaderLines',1,'Delimiter',',')*mpc.baseMVA;
    p_slack{k,1} = csvread([folder_path, 'bus\p_slack_up.csv'], 1, 0)*mpc.baseMVA-...
        csvread([folder_path, 'bus\p_slack_down.csv'], 1, 0)*mpc.baseMVA;

    # i_dc{k,1}=readmatrix([folders(k).folder,'\',folders(k).name,'\branchdc\i_from.csv'],'NumHeaderLines',1,'Delimiter',',');
    i_dc{k,1} = csvread([folder_path, 'branchdc\i_from.csv'], 1, 0);

    # p_ac{k,1}=readmatrix([folders(k).folder,'\',folders(k).name,'\branch\pf.csv'],'NumHeaderLines',1,'Delimiter',',')*mpc.baseMVA;
    p_ac{k,1} = csvread([folder_path, 'branch\pf.csv'], 1, 0)*mpc.baseMVA;

    % validTimesteps{k}=find(~isnan(...
    %     sum(ps_storage{k},2)+...
    %     sum(pg_gen{k},2)+...
    %     sum(pconv_conv{k},2)+...
    %     sum(iconv_conv{k},2)+...
    %     sum(p_slack{k},2)+...
    %     sum(i_dc{k},2)+...
    %     sum(p_ac{k},2)));

    validTimesteps{k}=1:size(pg_gen{k},1);
    KPI.convergence(k,1)=true;
    if validTimesteps{k}(end)<8760
        KPI.convergence(k,1)=false;
        warning(['scenario ',KPI.scenario_id{k},' has only ',num2str(validTimesteps{k}(end)),' time steps'])
    endif


endfor

%% calculate KPIs for each scenario
for scenario=1:length(KPI.scenario_id)

%% SEW

    % calculation of the maximum price
    maxPrice=zeros(length(validTimesteps{scenario}),1);
    for t=validTimesteps{scenario}

        % find dispatched generators
        index=find(pg_gen{scenario}(t,:)>0.01); % tolerance of 0.01 MW

        % remove indexes referred to non-dispatchable generators
        index=index(index<=size(mpc.gencost,2));

        % maximum price
        maxPrice(t)=max(mpc.gencost(index,6));
    endfor

    % calculation of SEW
    SEW=0;
    for k=1:size(pg_gen{scenario},2)

        pg_gen{scenario}(isnan(pg_gen{scenario}(:,k)),k)=0;

        if k<=size(mpc.gencost,2)
            SEW=SEW+sum(pg_gen{scenario}(validTimesteps{scenario},k).*(maxPrice-mpc.gencost(k,6)));
        else
            % if non-dispatchable, price is set to zero
            SEW=SEW+sum(pg_gen{scenario}(validTimesteps{scenario},k).*maxPrice);
        endif

    endfor

    % extrapolation for 1-year time horizon
    SEW=SEW*8760/length(validTimesteps{scenario});

    KPI.SEW(scenario,1)=SEW;

    %% CO2 emissions

    % emission table (missing, hypotesized proportional to gencost)
    mpc.CO2_emission_factors=[mpc.gencost(:,6);zeros(size(pg_gen{scenario},2)-size(mpc.gencost,1),1)];


    CO2=0;
    for k=1:size(pg_gen{scenario},2)
        pg_gen{scenario}(isnan(pg_gen{scenario}(:,k)),k)=0;
        CO2=CO2+sum(pg_gen{scenario}(validTimesteps{scenario},k).*mpc.CO2_emission_factors(k));
    endfor

    % extrapolation for 1-year time horizon
    CO2=CO2*8760/length(validTimesteps{scenario});

    KPI.CO2emission(scenario,1)=CO2;


    %% RES integration
    % to be corrected

    if isfield(mpc, 'ndgen')
      % buses to which RES are connected
      busRES=mpc.ndgen(:,1);

      % generation curtailment is negative
      curtailmentRES=p_slack{scenario}(validTimesteps{scenario},busRES);
      curtailmentRES(curtailmentRES<0.01)=0; % tolerance of 0.01 MW

      # RES=sum(curtailmentRES,'all');
      RES=sum(sum(curtailmentRES));

      % extrapolation for 1-year time horizon
      RES=RES*8760/length(validTimesteps{scenario});

      KPI.REScurtailment(scenario,1)=RES;
    endif

    %% nonCO2 emissions

    % emission table (missing, hypotesized proportional to gencost)
    mpc.nonCO2_emission_factors=[mpc.gencost(:,6);zeros(size(pg_gen{scenario},2)-size(mpc.gencost,1),1)];


    nonCO2=0;
    for k=1:size(pg_gen{scenario},2)
        pg_gen{scenario}(isnan(pg_gen{scenario}(:,k)),k)=0;
        nonCO2=nonCO2+sum(pg_gen{scenario}(validTimesteps{scenario},k).*mpc.nonCO2_emission_factors(k));
    endfor

    % extrapolation for 1-year time horizon
    nonCO2=nonCO2*8760/length(validTimesteps{scenario});

    KPI.nonCO2emission(scenario,1)=nonCO2;

    %% grid losses

    % ac branch losses
    Elosses_acbranch=0;
    p_ac{scenario}(isnan(p_ac{scenario}))=0;
    for k=1:size(mpc.branch,1)

        % index=find(mpc.bus(:,1)==mpc.branch(k,2));
        % voltage=mpc.bus(index,10);

        cosphi=0.95;

        Rpu=mpc.branch(k,3);

        Elosses_acbranch=Elosses_acbranch+sum(p_ac{scenario}(validTimesteps{scenario},k).^2)*Rpu/mpc.baseMVA/cosphi^2;
    endfor

    % extrapolation for 1-year time horizon
    Elosses_acbranch=Elosses_acbranch*8760/length(validTimesteps{scenario});


    % dc branch losses
    Elosses_dcbranch=0;
    i_dc{scenario}(isnan(i_dc{scenario}))=0;
    for k=1:size(mpc.branchdc,1)

        R=mpc.branchdc(k,3);
        R0=mpc.branchdc(k,13);

        Elosses_dcbranch=Elosses_dcbranch +...
            sum(i_dc{scenario}(validTimesteps{scenario},k*3-2).^2)*R*mpc.baseMVA +...
            sum(i_dc{scenario}(validTimesteps{scenario},k*3-1).^2)*R*mpc.baseMVA +...
            sum(i_dc{scenario}(validTimesteps{scenario},k*3-0).^2)*R0*mpc.baseMVA;
    endfor

    % extrapolation for 1-year time horizon
    Elosses_dcbranch=Elosses_dcbranch*8760/length(validTimesteps{scenario});


    %
    % % storage losses (necessary?)
    % ps_storage{scenario}(isnan(ps_storage{scenario}))=0;
    % Elosses_storage=-sum(ps_storage{scenario},'all');
    %
    % % extrapolation for 1-year time horizon
    % Elosses_storage=Elosses_storage*8760/length(validTimesteps{scenario});


    % ac-dc converter losses
    pconv_dc{scenario}(isnan(pconv_dc{scenario}))=0;
    powerDC=pconv_dc{scenario}(validTimesteps{scenario},:);
    # Elosses_dcconv=sum(powerDC,'all');
    Elosses_dcconv=sum(sum(powerDC));

    % extrapolation for 1-year time horizon
    Elosses_dcconv=Elosses_dcconv*8760/length(validTimesteps{scenario});



    % transformer losses of ac-dc converter
    Elosses_dctrafo=0;
    for k=1:size(mpc.convdc,1)

        cosphi=0.95;

        Rpu=mpc.convdc(k,9);

        Elosses_dctrafo=Elosses_dctrafo+...
            sum(pconv_dc{scenario}(validTimesteps{scenario},k*2-1).^2)*Rpu/mpc.baseMVA/cosphi^2+...
            sum(pconv_dc{scenario}(validTimesteps{scenario},k*2-0).^2)*Rpu/mpc.baseMVA/cosphi^2;
    endfor

    % extrapolation for 1-year time horizon
    Elosses_dctrafo=Elosses_dctrafo*8760/length(validTimesteps{scenario});


    Elosses=Elosses_acbranch+Elosses_dcbranch+Elosses_dcconv+Elosses_dctrafo;

    KPI.Elosses(scenario,1)=Elosses;

    %% Adequacy

    % generation curtailment is negative
    curtailmentLoad=p_slack{scenario}(validTimesteps{scenario},:);
    curtailmentLoad(curtailmentLoad<0.01)=0; % tolerance of 0.01 MW

    # EENS=sum(curtailmentLoad,'all');
    EENS=sum(sum(curtailmentLoad));

    LOLE=sum(sum(curtailmentLoad,2)>0);

    % extrapolation for 1-year time horizon
    EENS=EENS*8760/length(validTimesteps{scenario});
    LOLE=LOLE*8760/length(validTimesteps{scenario});

    KPI.EENS_MWh(scenario,1)=EENS;
    KPI.LOLE(scenario,1)=LOLE;

    %% Power (Energy) transferred by the HVDC grids

    powerDC=pconv_dc{scenario}(validTimesteps{scenario},:);
    # lossesDC=sum(powerDC,'all');
    lossesDC=sum(sum(powerDC));

    energyDC=sum(powerDC(powerDC>0))-lossesDC;


    % extrapolation for 1-year time horizon
    energyDC=energyDC*8760/length(validTimesteps{scenario});

    KPI.energyDC(scenario,1)=energyDC;

endfor

%% create table of results and save them to excel

# KPI_table=struct2table(KPI); # Not implemented in Octave
# writetable(KPI_table,[user_results_dir, 'KPI_results.xlsx'],'sheet',macro_scenario)

KPI_names = fieldnames(KPI);
nCols = length(KPI_names);
nRows = length(KPI.scenario_id);

data = cell(nRows+1, nCols);
data(1, :) = KPI_names;

for i = 1:nRows
    for j = 1:nCols
      data_ij = KPI.(KPI_names{j})(i);
      if isa(data_ij, 'cell')
        data_ij = data_ij{1};
      endif
      data{i+1, j} = data_ij;
    endfor
endfor

xlswrite([user_results_dir, '\KPI_results.xlsx'], data);

disp('KPI computation completed')

