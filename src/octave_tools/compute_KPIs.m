% Code developed by Marco ROSSI (RSE) and Nicolas BARLA (SuperGrid Institute)
clc
vars_to_delete = setdiff(who, {'work_dir'});
if length(vars_to_delete) > 0
    clear(vars_to_delete{:});  % clear all variables except work_dir
endif

addpath(fileparts(mfilename('fullpath')));
pkg load io

MWh_to_TWh = 10^-6;
t_to_Mt = 10^-6;
euro_to_Meuro = 10^-6;

% work_dir is a global variable defined inside the run comand in the Julia code
% work_dir = strjoin(parent_names(1:end-2), '\');
assert(exist('work_dir', 'var') == 1, "work_dir is not defined")

simulation_dir = [work_dir, '\simulation_interface'];
macro_scenario = NaN;
folders=dir(simulation_dir);
folders(1:2)=[];
for i=1:length(folders)
    fifo_name = folders(i).name;
    if fifo_name(end-1:end) == '.m';
        macro_scenario = fifo_name(1:end-2);
    endif
endfor
assert(isnan(macro_scenario) == 0, strcat("No .m file has been found in ", simulation_dir))

simulation_results_dir = [simulation_dir, '\', macro_scenario];
user_results_dir = [work_dir, '\user_interface\results'];

%% load grid model (.m file)
addpath(simulation_dir);
mpc = feval(macro_scenario);
rmpath(simulation_dir);

% Matching load_id and bus_id
bus_to_load = [];  % 1 load id per bus. 0 means no load.
load_to_bus = [];  % bus id per load
load_to_load_extra = [];  % Row of the mpc.load_extra matrix where the flexibility costs (shifting, reduction, curtailment) are stored. 1 value per load. 0 means fixed load with default curtailment cost (3000 €/MWh).

load_id = 0;
for k=1:size(mpc.bus,1)
    if mpc.bus(k,3) != 0
        load_id = load_id + 1;
        bus_to_load = [bus_to_load; load_id];
        load_to_bus = [load_to_bus; k];
        if isfield(mpc, 'load_extra') && sum(mpc.load_extra(:,1)==load_id) > 0
            load_extra_id = find(mpc.load_extra(:,1)==load_id);
            assert(length(load_extra_id)==1, strcat("mpc.load_extra contains several times load_id==", load_id, "\nIt is not allowed"))  % Curtailed power is known per bus and curtailment cost is known per load.
            load_to_load_extra = [load_to_load_extra; load_extra_id];
        else
            load_to_load_extra = [load_to_load_extra; 0];  % Fixed load with default curtailment cost
        endif
    else
        bus_to_load = [bus_to_load; 0];  % Not a load
    endif
endfor


folders=dir(simulation_results_dir);
folders(1:2)=[];

for scenario_id=1:length(folders)
    folder_path = [folders(scenario_id).folder, '\', folders(scenario_id).name, '\'];
    scenario_name=folders(scenario_id).name;
    KPI.scenario_name{scenario_id} = scenario_name;
    if length(dir(folder_path)) > 2 % '.' and '..' and in an empty folder

        disp(strcat('KPI computation for scenario', num2str(scenario_id), ': ', scenario_name))

        assert(isfolder([folder_path, 'bus\']) == 1, strcat(folder_path, '\', 'bus', ' does not exist'))
        p_slack = csvread([folder_path, 'bus\p_slack_up.csv'], 1, 0)*mpc.baseMVA-...
            csvread([folder_path, 'bus\p_slack_down.csv'], 1, 0)*mpc.baseMVA;

        assert(isfolder([folder_path, 'gen\']) == 1, strcat(folder_path, '\', 'gen', 'does not exist'))
        pg_gen = csvread([folder_path, 'gen\pg.csv'], 1, 0)*mpc.baseMVA;
        pg_curt = csvread([folder_path, 'gen\pgcurt.csv'], 1, 0)*mpc.baseMVA;

        assert(isfolder([folder_path, 'load\']) == 1, strcat(folder_path, '\', 'load', 'does not exist'))
        pflex = csvread([folder_path, 'load\pflex.csv'], 1, 0)*mpc.baseMVA;
        pred = csvread([folder_path, 'load\pred.csv'], 1, 0)*mpc.baseMVA;
        pshift_up = csvread([folder_path, 'load\pshift_up.csv'], 1, 0)*mpc.baseMVA;
        pshift_down = csvread([folder_path, 'load\pshift_down.csv'], 1, 0)*mpc.baseMVA;


        if isfolder([folder_path, 'storage\'])
            ps_storage = csvread([folder_path, 'storage\ps.csv'], 1, 0)*mpc.baseMVA;  % 1 and 0 means starting from row 2 and col 1
        endif

        if isfolder([folder_path, 'convdc\'])
            pconv_dc_with_headers = csvread([folder_path, 'convdc\pconv.csv'], 0, 0);
            pconv_dc_headers = pconv_dc_with_headers(1,:);
            nb_rows = size(pconv_dc_with_headers,1);
            pconv_dc = pconv_dc_with_headers(2:nb_rows,:)*mpc.baseMVA;
            % Remove columns with only NaN
            pconv_dc_headers = pconv_dc_headers(:,~all(isnan(pconv_dc)));
            pconv_dc = pconv_dc(:,~all(isnan(pconv_dc)));
        endif

        if isfolder([folder_path, 'branch\'])
            p_ac = csvread([folder_path, 'branch\pf.csv'], 1, 0)*mpc.baseMVA;
        endif

        if isfolder([folder_path, 'branchdc\'])
            i_dc_with_headers = csvread([folder_path, 'branchdc\i_from.csv'], 0, 0);
            i_dc_headers = i_dc_with_headers(1,:);
            nb_rows = size(i_dc_with_headers,1);
            i_dc = i_dc_with_headers(2:nb_rows,:);
            % Remove columns with only NaN
            i_dc_headers = i_dc_headers(:,~all(isnan(i_dc)));
            i_dc = i_dc(:,~all(isnan(i_dc)));
        endif

        validTimesteps=1:size(pg_gen,1);
        KPI.convergence(scenario_id,1)=true;
        if validTimesteps(end)<8760
            KPI.convergence(scenario_id,1)=false;
            warning(['scenario ', scenario_name, ' has only ', num2str(validTimesteps(end)), ' time steps'])
        endif


        %% --- Economic indicators ---

        number_of_dispatchable_gen = size(mpc.gencost,1);

        % calculation of the maximum price
        maxPrice=zeros(length(validTimesteps),1);  % Market price each hour, considering 1 market zone
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

        cost_generation = 0;  % Generation cost (€) of the generators
        cost_curtailment_ndgen = 0;  % Cost (€) of RES power curtailment
        cost_curtailment_load = 0;  % Cost (€) of load shedding
        cost_flexibility = 0;  % Cost (€) of load shifting & load voluntary reduction
        welfare_producers = 0;  % Benefits/Surplus (€) for producers (production * (market price - generation cost)) & storages ((production - consumption) * market price)
        welfare_consumers = 0;  % Benefits/Surplus (€) for consumers (consumption * (curtailment cost - market price))

        % Dispatchable generators
        for k_gen = 1:min(size(pg_gen,2), number_of_dispatchable_gen)
            pg_gen(isnan(pg_gen(:,k_gen)),k_gen)=0;  % MW, vector(time)
            cost_generation = cost_generation + sum(pg_gen(validTimesteps,k_gen).*mpc.gencost(k_gen,5));  % €
            welfare_producers = welfare_producers + sum(pg_gen(validTimesteps,k_gen).*(maxPrice-mpc.gencost(k_gen,5)));  % €
        endfor

        % Non-dispatchable generators
        if isfield(mpc, 'ndgen')
            for k_ndgen = number_of_dispatchable_gen+1:size(pg_gen,2)
                pg_gen(isnan(pg_gen(:,k_ndgen)),k_ndgen)=0;  % MW, vector(time)
                cost_generation = cost_generation + sum(pg_gen(validTimesteps,k_ndgen).*mpc.ndgen(k_ndgen-number_of_dispatchable_gen,6));  % €
                cost_curtailment_ndgen = cost_curtailment_ndgen + sum(pg_curt(validTimesteps,k_ndgen).*mpc.ndgen(k_ndgen-number_of_dispatchable_gen,7));  % €
                welfare_producers = welfare_producers + sum(pg_gen(validTimesteps,k_ndgen).*(maxPrice-mpc.ndgen(k_ndgen-number_of_dispatchable_gen,6)));  % €
                % Should we do welfare_producers = welfare_producers - cost_curtailment_ndgen ?
            endfor
        endif

        % Storages
        if isfield(mpc, 'storage')
            for k = 1:size(ps_storage,2)
                ps_storage(isnan(ps_storage(:,k)),k)=0;  % MW, vector(time)
                welfare_producers = welfare_producers - sum(ps_storage(validTimesteps,k).*maxPrice);
            endfor
        endif

        % Loads
        for k = 1:size(pflex,2)
            pflex(isnan(pflex(:,k)),k)=0;  % MW, vector(time)
            pred(isnan(pred(:,k)),k)=0;  % MW, vector(time)
            pshift_up(isnan(pshift_up(:,k)),k)=0;  % MW, vector(time)
            pshift_down(isnan(pshift_down(:,k)),k)=0;  % MW, vector(time)

            bus_id = load_to_bus(k);
            load_extra_id = load_to_load_extra(k);
            if load_extra_id == 0  % Fixed load with default curtailment cost
                cost_curt = 3000;  % Default curtailment cost for fixed loads (€/MWh)
            else
                cost_curt = mpc.load_extra(load_extra_id, 12);  % €/MWh
                cost_shift = mpc.load_extra(load_extra_id, 10);  % €/MWh
                cost_red = mpc.load_extra(load_extra_id, 11);  % €/MWh
                cost_flexibility = cost_flexibility + sum(pred(validTimesteps,k).*cost_red) + sum(pshift_up(validTimesteps,k).*cost_shift);  % €.   % Shifted down and shifted up energy are supposed to be equal
            endif
            curtailed_power = p_slack(validTimesteps,bus_id);  % MW, vector(time)
            consumed_power = pflex(validTimesteps,k) - curtailed_power + pshift_up(validTimesteps,k) - pshift_down(validTimesteps,k);  % MW, vector(time)
            cost_curtailment_load = cost_curtailment_load + sum(curtailed_power.*cost_curt);  % €
            welfare_consumers = welfare_consumers + sum(consumed_power.*(cost_curt - maxPrice));  % €
        endfor

        % extrapolation for 1-year time horizon
        KPI.Generation_cost(scenario_id,1) = cost_generation * 8760 / length(validTimesteps);  % €/y
        KPI.RES_curtailment_cost(scenario_id,1) = cost_curtailment_ndgen * 8760 / length(validTimesteps);  % €/y
        KPI.Load_shedding_cost(scenario_id,1) = cost_curtailment_load * 8760 / length(validTimesteps);  % €/y
        KPI.Load_flexibility_cost(scenario_id,1) = cost_flexibility * 8760 / length(validTimesteps);  % €/y
        KPI.Producers_surplus(scenario_id,1) = welfare_producers * 8760 / length(validTimesteps);  % €/y
        KPI.Consumers_surplus(scenario_id,1) = welfare_consumers * 8760 / length(validTimesteps);  % €/y


        %% --- CO2 emissions ---
        CO2=0;
        for k=1:number_of_dispatchable_gen  % Use k=1:size(pg_gen,2) to consider also ndgen emissions
            pg_gen(isnan(pg_gen(:,k)),k)=0;  % MW, vector(time)
            CO2=CO2+sum(pg_gen(validTimesteps,k).*mpc.emission_factors(k,1));  % tCO2
        endfor

        % extrapolation for 1-year time horizon
        CO2=CO2*8760/length(validTimesteps);  % tCO2/y

        KPI.CO2emission(scenario_id,1)=CO2;


        %% --- RES integration ---

        if isfield(mpc, 'ndgen')
            % buses to which RES are connected
            busRES=mpc.ndgen(:,1);

            % replace NaN with 0   -- newly added
            pg_curt(isnan(pg_curt))=0;  % MW, vector(time)

            % generation curtailment
            curtailmentRES=pg_curt(validTimesteps,:);  % MW, vector(time)
            curtailmentRES(curtailmentRES<0.01)=0; % tolerance of 0.01 MW

            % RES=sum(curtailmentRES,'all');
            RES=sum(sum(curtailmentRES));  % MWh

            % extrapolation for 1-year time horizon
            RES=RES*8760/length(validTimesteps);  % MWh/y
        else
            RES = 0;
        endif
        KPI.REScurtailment(scenario_id,1)=RES;


        %% --- Grid losses ---

        % - AC branch losses -
        Elosses_acbranch=0;
        if exist('i_dc', 'var')
            p_ac(isnan(p_ac))=0;  % MW, vector(time)
            for k=1:size(mpc.branch,1)
                cosphi=0.95;

                Rpu=mpc.branch(k,3);  % Rpu = R / R_base, with R_base = V_base^2/mpc.baseMVA

                Elosses_acbranch=Elosses_acbranch+sum(p_ac(validTimesteps,k).^2)*Rpu/mpc.baseMVA/cosphi^2;  % MWh
            endfor
            Elosses_acbranch=Elosses_acbranch*8760/length(validTimesteps); % extrapolation for 1-year time horizon (MWh/y)
        endif

        % - DC branch losses -
        Elosses_dcbranch=0;
        if exist('i_dc', 'var')
            i_dc(isnan(i_dc))=0;

            column = 1;
            for k=1:size(mpc.branchdc ,1)

                R=mpc.branchdc(k,3);  % pu, with R_base = V_base^2/mpc.baseMVA
                R0=mpc.branchdc(k,12);  % pu, with R_base = V_base^2/mpc.baseMVA

                assert(i_dc_headers(1, column)==k);
                assert(i_dc_headers(1, column + 1)==k);

                if column + 2 <= size(i_dc_headers, 2) && i_dc_headers(1, column + 2)==k  % Bipolar branch has 3 columns
                    Elosses_dcbranch=Elosses_dcbranch +...
                        sum(i_dc(validTimesteps, column).^2)*R*mpc.baseMVA +...
                        sum(i_dc(validTimesteps, column + 1).^2)*R*mpc.baseMVA +...
                        sum(i_dc(validTimesteps, column + 2).^2)*R0*mpc.baseMVA;  % MWh
                    column = column + 3;
                else
                    Elosses_dcbranch=Elosses_dcbranch +...
                        sum(i_dc(validTimesteps, column).^2)*R*mpc.baseMVA +...
                        sum(i_dc(validTimesteps, column + 1).^2)*R*mpc.baseMVA;  % MWh
                    column = column + 2;
                endif
            endfor
            Elosses_dcbranch=Elosses_dcbranch*8760/length(validTimesteps); % extrapolation for 1-year time horizon (MWh/y)
        endif

        % - Converter & Transformer losses -
        if exist('pconv_dc', 'var')
            pconv_dc(isnan(pconv_dc))=0;

            % ac-dc converter losses
            powerDC=pconv_dc(validTimesteps,:);  % MW, vector(time)
            Elosses_dcconv=sum(sum(powerDC));  % MWh
            Elosses_dcconv=Elosses_dcconv*8760/length(validTimesteps); % extrapolation for 1-year time horizon (MWh/y)

            % transformer losses of ac-dc converter
            Elosses_dctrafo=0;
            column=1;
            for k=1:size(mpc.convdc,1)
                cosphi=0.95;
                Rpu=mpc.convdc(k,9); % column rtf in mpc.convdc. pu, with R_base = V_base^2/mpc.baseMVA

                assert(pconv_dc_headers(1, column)==k);
                if column + 1 <= size(pconv_dc_headers, 2) && pconv_dc_headers(1, column + 1)==k  % Bipolar has two columns
                    Elosses_dctrafo=Elosses_dctrafo+...
                        sum(pconv_dc(validTimesteps, column).^2)*Rpu/mpc.baseMVA/cosphi^2+...
                        sum(pconv_dc(validTimesteps, column + 1).^2)*Rpu/mpc.baseMVA/cosphi^2;  % MWh
                    column = column + 2;
                else
                    Elosses_dctrafo=Elosses_dctrafo+...
                        sum(pconv_dc(validTimesteps, column).^2)*Rpu/mpc.baseMVA/cosphi^2;  % MWh
                    column = column + 1;
                endif
            endfor
            Elosses_dctrafo=Elosses_dctrafo*8760/length(validTimesteps); % extrapolation for 1-year time horizon (MWh/y)
        else
            Elosses_dcconv = 0;
            Elosses_dctrafo = 0;
        endif

        Elosses=Elosses_acbranch+Elosses_dcbranch+Elosses_dcconv+Elosses_dctrafo;  % MWh/y

        KPI.Energy_losses(scenario_id,1)=Elosses;


        %% --- Adequacy (ENS & LOLE) ---

        % Generation curtailment is negative
        curtailmentLoad=p_slack(validTimesteps,:);  % MW, vector(time)
        curtailmentLoad(curtailmentLoad<0.01)=0; % tolerance of 0.01 MW

        ENS=sum(sum(curtailmentLoad));  % MWh
        LOLE=sum(sum(curtailmentLoad,2)>0);  % h

        % extrapolation for 1-year time horizon
        ENS=ENS*8760/length(validTimesteps);  % MWh/y
        LOLE=LOLE*8760/length(validTimesteps);  % h/y

        KPI.ENS(scenario_id,1)=ENS;
        KPI.LOLE(scenario_id,1)=LOLE;


        %% --- Power (Energy) transferred by the HVDC grids ---
        if exist('pconv_dc', 'var')
            powerDC=pconv_dc(validTimesteps,:);  % MW, vector(time)
            lossesDC=sum(sum(powerDC));  % MWh
            energyDC=sum(powerDC(powerDC>0))-lossesDC;  % MWh
            energyDC=energyDC*8760/length(validTimesteps); % extrapolation for 1-year time horizon (MWh/y)
        else
            energyDC = 0;
        endif
        KPI.energyDC(scenario_id,1)=energyDC;
    else
        KPI.Generation_cost(scenario_id,1) = nan;
        KPI.RES_curtailment_cost(scenario_id,1) = nan;
        KPI.Load_shedding_cost(scenario_id,1) = nan;
        KPI.Load_flexibility_cost(scenario_id,1) = nan;
        KPI.Producers_surplus(scenario_id,1) = nan;
        KPI.Consumers_surplus(scenario_id,1) = nan;
        KPI.CO2emission(scenario_id,1) = nan;
        KPI.REScurtailment(scenario_id,1) = nan;
        KPI.Energy_losses(scenario_id,1) = nan;
        KPI.energyDC(scenario_id,1) = nan;
        KPI.ENS(scenario_id,1) = nan;
        KPI.LOLE(scenario_id,1) = nan;
    endif
endfor


%% --- Save results ---

% - Convert KPI units -

KPI_names = fieldnames(KPI);
expected_KPI_names = {"scenario_name", "convergence", "Generation_cost", "RES_curtailment_cost", "Load_shedding_cost", "Load_flexibility_cost", "Producers_surplus", "Consumers_surplus", "CO2emission", "REScurtailment", "Energy_losses", "energyDC", "ENS", "LOLE"};
units = {"", "", "M€/y", "M€/y", "M€/y", "M€/y", "M€/y", "M€/y", "Mt/y", "TWh/y", "TWh/y", "TWh/y", "TWh/y", "h/y"};
missing_KPIs = setdiff(expected_KPI_names, KPI_names);
unexpected_KPIs = setdiff(KPI_names, expected_KPI_names);
if length(missing_KPIs) > 0
    display(missing_KPIs)
    error("Missing KPIs names")
elseif length(unexpected_KPIs) > 0
    display(unexpected_KPIs)
    error("Unexpected KPI names")
endif

for j = 3:8
    KPI.(expected_KPI_names{j}) = KPI.(expected_KPI_names{j}) * euro_to_Meuro;  % M€/y
endfor
KPI.CO2emission = KPI.CO2emission * t_to_Mt;  % MtCO2/y
for j = 10:13
    KPI.(expected_KPI_names{j}) = KPI.(expected_KPI_names{j}) * MWh_to_TWh;  % TWh/y
endfor
KPI.LOLE = KPI.LOLE * 1;  % h/y

% - Create table of results and save them to excel -

nCols = length(expected_KPI_names);
nRows = length(KPI.scenario_name);

data = cell(nRows+1, nCols);
data(1, :) = expected_KPI_names;
data(2, :) = units;

for i = 1:nRows
    for j = 1:nCols
        data_ij = KPI.(expected_KPI_names{j})(i);  % expected_KPI_names is used instead of KPI_names to ensure coherence with units
        if isa(data_ij, 'cell')
        data_ij = data_ij{1};
        endif
        data{i+2, j} = data_ij;
    endfor
endfor

xlswrite([user_results_dir, '\KPI_results.xlsx'], data);

disp('KPI computation completed')
