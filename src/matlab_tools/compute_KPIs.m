clear
clc

pkg load io

MWh_to_TWh = 10^-6;
t_to_Mt = 10^-6;
euro_to_Meuro = 10^-6;

inputfile='grid_model'; # grid_model.m has to be in the same folder
simulation_results_dir = fileread('simulation_results_path.txt');
parent_names = strsplit(simulation_results_dir, '\');
macro_scenario = parent_names(end);
work_dir = strjoin(parent_names(1:end-2), '\');
user_results_dir = [work_dir, '\user_interface\results'];

%% load simulation input/output data
mpc=eval(inputfile);

folders=dir(simulation_results_dir);
folders(1:2)=[];

for scenario_id=1:length(folders)
    folder_path = [folders(scenario_id).folder, '\', folders(scenario_id).name, '\'];
    scenario_name=folders(scenario_id).name;
    KPI.scenario_name{scenario_id} = scenario_name;
    if length(dir(folder_path)) > 2 % '.' and '..' and in an empty folder

        disp(strcat('KPI computation for scenario', num2str(scenario_id), ': ', scenario_name))

        assert(isfolder([folder_path, 'bus\']), strcat(folder_path, '\', 'bus', ' does not exist'))
        p_slack = csvread([folder_path, 'bus\p_slack_up.csv'], 1, 0)*mpc.baseMVA-...
            csvread([folder_path, 'bus\p_slack_down.csv'], 1, 0)*mpc.baseMVA;

        assert(isfolder([folder_path, 'gen\']), strcat(folder_path, '\', 'gen', 'does not exist'))
        pg_gen = csvread([folder_path, 'gen\pg.csv'], 1, 0)*mpc.baseMVA;
        pg_curt = csvread([folder_path, 'gen\pgcurt.csv'], 1, 0)*mpc.baseMVA;

        assert(isfolder([folder_path, 'load\']), strcat(folder_path, '\', 'load', 'does not exist'))
        pflex_load = csvread([folder_path, 'load\pflex.csv'], 1, 0)*mpc.baseMVA;

        if isfolder([folder_path, 'storage\'])
            ps_storage = csvread([folder_path, 'storage\ps.csv'], 1, 0)*mpc.baseMVA;  # 1 and 0 means starting from row 2 and col 1
        endif

        if isfolder([folder_path, 'convdc\'])
            pconv_dc = csvread([folder_path, 'convdc\pconv.csv'], 1, 0)*mpc.baseMVA;
        endif

        if isfolder([folder_path, 'branch\'])
            p_ac = csvread([folder_path, 'branch\pf.csv'], 1, 0)*mpc.baseMVA;
        endif

        if isfolder([folder_path, 'branchdc\'])
            i_dc = csvread([folder_path, 'branchdc\i_from.csv'], 1, 0);
        endif

        validTimesteps=1:size(pg_gen,1);
        KPI.convergence(scenario_id,1)=true;
        if validTimesteps(end)<8760
            KPI.convergence(scenario_id,1)=false;
            warning(['scenario ', scenario_name, ' has only ', num2str(validTimesteps(end)), ' time steps'])
        endif


        %% --- SEW ---

        number_of_dispatchable_gen = size(mpc.gencost,1);

        % calculation of the maximum price
        maxPrice=zeros(length(validTimesteps),1);
        for t=validTimesteps

            % find dispatched generators
            index=find(pg_gen(t,:)>0.01); % tolerance of 0.01 MW

            % maximum price
            index_gen=index(index <= number_of_dispatchable_gen);
            if length(index_gen) > 0
                maxPrice_gen = max(mpc.gencost(index_gen,5)); % €/MWh
            else
                maxPrice_gen = nan;
            endif

            index_ndgen=index(index > number_of_dispatchable_gen);
            if length(index_ndgen) > 0
                maxPrice_ndgen = max(mpc.ndgen(index_ndgen-number_of_dispatchable_gen,6)); % €/MWh
            else
                maxPrice_ndgen = nan;
            endif

            maxPrice(t) = max(maxPrice_gen, maxPrice_ndgen);
            assert(~isnan(maxPrice(t)))
        endfor

        % calculation of SEW
        SEW=0;

        % Dispatchable generators
        for gen = 1:min(size(pg_gen,1), number_of_dispatchable_gen)
            pg_gen(isnan(pg_gen(:,gen)),gen)=0;
            SEW = SEW + sum(pg_gen(validTimesteps,gen).*(maxPrice-mpc.gencost(gen,5)));
        endfor

        % Non-dispatchable generators
        for ndgen = number_of_dispatchable_gen+1:size(pg_gen,2)
            pg_gen(isnan(pg_gen(:,ndgen)),ndgen)=0;
            SEW = SEW + sum(pg_gen(validTimesteps,ndgen).*(maxPrice-mpc.ndgen(ndgen-number_of_dispatchable_gen,6)));
            SEW = SEW - sum(pg_curt(validTimesteps,ndgen).*mpc.ndgen(ndgen-number_of_dispatchable_gen,7));
        endfor

        SEW = SEW * 8760 / length(validTimesteps); % extrapolation for 1-year time horizon

        KPI.SEW(scenario_id,1) = SEW;

        %% --- CO2 emissions ---
        CO2=0;
        for k=1:size(pg_gen,2)
            pg_gen(isnan(pg_gen(:,k)),k)=0;
            CO2=CO2+sum(pg_gen(validTimesteps,k).*mpc.emission_factors(k,1));
        endfor

        % extrapolation for 1-year time horizon
        CO2=CO2*8760/length(validTimesteps);

        KPI.CO2emission(scenario_id,1)=CO2;


        %% --- RES integration ---

        if isfield(mpc, 'ndgen')
            % buses to which RES are connected
            busRES=mpc.ndgen(:,1);

            % replace NaN with 0   -- newly added
            pg_curt(isnan(pg_curt))=0;

            % generation curtailment
            curtailmentRES=pg_curt(validTimesteps,:);
            curtailmentRES(curtailmentRES<0.01)=0; % tolerance of 0.01 MW

            # RES=sum(curtailmentRES,'all');
            RES=sum(sum(curtailmentRES));

            % extrapolation for 1-year time horizon
            RES=RES*8760/length(validTimesteps);
        else
            RES = 0;
        endif
        KPI.REScurtailment(scenario_id,1)=RES;


        %% --- Grid losses ---

        % - AC branch losses -
        Elosses_acbranch=0;
        if exist('i_dc', 'var')
            p_ac(isnan(p_ac))=0;
            for k=1:size(mpc.branch,1)
                cosphi=0.95;

                Rpu=mpc.branch(k,3);

                Elosses_acbranch=Elosses_acbranch+sum(p_ac(validTimesteps,k).^2)*Rpu/mpc.baseMVA/cosphi^2;
            endfor
            Elosses_acbranch=Elosses_acbranch*8760/length(validTimesteps); % extrapolation for 1-year time horizon
        endif

        % - DC branch losses -
        Elosses_dcbranch=0;
        if exist('i_dc', 'var')
            i_dc(isnan(i_dc))=0;
            for k=1:size(mpc.branchdc,1)

                R=mpc.branchdc(k,3);
                R0=mpc.branchdc(k,13);

                Elosses_dcbranch=Elosses_dcbranch +...
                    sum(i_dc(validTimesteps,k*3-2).^2)*R*mpc.baseMVA +...
                    sum(i_dc(validTimesteps,k*3-1).^2)*R*mpc.baseMVA +...
                    sum(i_dc(validTimesteps,k*3-0).^2)*R0*mpc.baseMVA;
            endfor
            Elosses_dcbranch=Elosses_dcbranch*8760/length(validTimesteps); % extrapolation for 1-year time horizon
        endif

        % - Converter & Transformer losses -
        if exist('pconv_dc', 'var')
            pconv_dc(isnan(pconv_dc))=0;

            % ac-dc converter losses
            powerDC=pconv_dc(validTimesteps,:);
            Elosses_dcconv=sum(sum(powerDC));
            Elosses_dcconv=Elosses_dcconv*8760/length(validTimesteps); % extrapolation for 1-year time horizon

            % transformer losses of ac-dc converter
            Elosses_dctrafo=0;
            for k=1:size(mpc.convdc,1)
                cosphi=0.95;
                Rpu=mpc.convdc(k,9);
                Elosses_dctrafo=Elosses_dctrafo+...
                    sum(pconv_dc(validTimesteps,k*2-1).^2)*Rpu/mpc.baseMVA/cosphi^2+...
                    sum(pconv_dc(validTimesteps,k*2-0).^2)*Rpu/mpc.baseMVA/cosphi^2;
            endfor
            Elosses_dctrafo=Elosses_dctrafo*8760/length(validTimesteps); % extrapolation for 1-year time horizon
        else
            Elosses_dcconv = 0;
            Elosses_dctrafo = 0;
        endif

        Elosses=Elosses_acbranch+Elosses_dcbranch+Elosses_dcconv+Elosses_dctrafo;

        KPI.Elosses(scenario_id,1)=Elosses;


        %% --- Adequacy (ENS & LOLE) ---

        % Generation curtailment is negative
        curtailmentLoad=p_slack(validTimesteps,:);
        curtailmentLoad(curtailmentLoad<0.01)=0; % tolerance of 0.01 MW

        ENS=sum(sum(curtailmentLoad));
        LOLE=sum(sum(curtailmentLoad,2)>0);

        % extrapolation for 1-year time horizon
        ENS=ENS*8760/length(validTimesteps);
        LOLE=LOLE*8760/length(validTimesteps);

        KPI.ENS(scenario_id,1)=ENS;
        KPI.LOLE(scenario_id,1)=LOLE;


        %% --- Power (Energy) transferred by the HVDC grids ---
        if exist('pconv_dc', 'var')
            powerDC=pconv_dc(validTimesteps,:);
            lossesDC=sum(sum(powerDC));
            energyDC=sum(powerDC(powerDC>0))-lossesDC;
            energyDC=energyDC*8760/length(validTimesteps); % extrapolation for 1-year time horizon
        else
            energyDC = 0;
        endif
        KPI.energyDC(scenario_id,1)=energyDC;
    else
        KPI.SEW(scenario_id,1) = nan;
        KPI.REScurtailment(scenario_id,1) = nan;
        KPI.CO2emission(scenario_id,1) = nan;
        KPI.Elosses(scenario_id,1) = nan;
        KPI.ENS(scenario_id,1) = nan;
        KPI.LOLE(scenario_id,1) = nan;
        KPI.energyDC(scenario_id,1) = nan;
    endif
endfor


%% --- Save results ---

% - Convert KPI units -

KPI_names = fieldnames(KPI);
expected_KPI_names = {"scenario_name", "convergence", "SEW", "REScurtailment", "CO2emission", "Elosses", "ENS", "LOLE", "energyDC"};
missing_KPIs = setdiff(expected_KPI_names, KPI_names);
unexpected_KPIs = setdiff(KPI_names, expected_KPI_names);
if length(missing_KPIs) > 0
  display(missing_KPIs)
  error("Missing KPIs names")
elseif length(unexpected_KPIs) > 0
  display(unexpected_KPIs)
  error("Unexpected KPI names")
endif
units = {"", "", "M€/y", "TWh/y", "Mt/y", "TWh/y", "TWh/y", "h/y", "TWh/y"};

KPI.SEW = KPI.SEW * MWh_to_TWh;
KPI.REScurtailment = KPI.REScurtailment * MWh_to_TWh;
KPI.CO2emission = KPI.CO2emission * t_to_Mt;
KPI.Elosses = KPI.Elosses * MWh_to_TWh;
KPI.ENS = KPI.ENS * MWh_to_TWh;
KPI.ENS = KPI.LOLE * 1;
KPI.energyDC = KPI.energyDC * MWh_to_TWh;

% - Create table of results and save them to excel -

nCols = length(KPI_names);
nRows = length(KPI.scenario_name);

data = cell(nRows+1, nCols);
data(1, :) = KPI_names;
data(2, :) = units;

for i = 1:nRows
    for j = 1:nCols
      data_ij = KPI.(KPI_names{j})(i);
      if isa(data_ij, 'cell')
        data_ij = data_ij{1};
      endif
      data{i+2, j} = data_ij;
    endfor
endfor

xlswrite([user_results_dir, '\KPI_results.xlsx'], data);

disp('KPI computation completed')

