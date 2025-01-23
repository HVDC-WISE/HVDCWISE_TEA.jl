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

    line_titles = {"Attribute", "Unit", "Description"};
    for i = 1:3
        assert(sheet{i,1} == line_titles{i}, strcat("Cell A", num2str(i), " in ", reliability_data_path, " should be ", line_titles{i}, ", not ", sheet{i,1}))
    endfor

    attributes = {"Attribute", "MTTR", "FOR", "Failure rate 0", "Failure rate per km", "Correlation between poles"};
    units = {"Unit", "h", "%", "1/y", "1/y/km", "Number between 0 & 1"};
    for j = 2:6
        assert(strcmp(sheet{1,j}, attributes{j}), strcat("First cell of column ", num2str(j), " in ", reliability_data_path, " should be ", attributes{j}, ", not ", sheet{1,j}))
        assert(strcmp(sheet{2,j}, units{j}), strcat("Unit for attribute ", attributes{j}, " in ", reliability_data_path, " should be ", units{j}, ", not ", sheet{2,j}))
    end

    % Input data
    input_data = struct();

    % Components without a length
    components_localised = {"Generator", "Transformer", "Converter ACDC"};
    for i = 4:6
        comp_name = components_localised{i-3};
        assert(sheet{i,1} == comp_name, strcat("Cell A", num2str(i), " in ", reliability_data_path, " should be ", comp_name, ", not ", sheet{i,1}))
        for j = 2:3
            assert(isempty(sheet{i,j}) == 0, strcat("No cell in B", num2str(i), ":C", num2str(i), " in ", reliability_data_path, " should be empty. ", mat2str(sheet{i,2:3})))
        endfor
        comp_name = strrep(comp_name, " ", "_");
        input_data.(comp_name) = struct();
        input_data.(comp_name).MTTR = sheet{i,2};
        input_data.(comp_name).FOR = sheet{i,3};
    endfor

    % Components with a length
    components_linear = {"AC OHL", "DC OHL", "DC cable"};
    for i = 7:9
        comp_name = components_linear{i-6};
        assert(strcmp(sheet{i,1}, comp_name), strcat("Cell A", num2str(i), " in ", reliability_data_path, " should be ", comp_name, ", not ", sheet{i,1}))
        assert(isempty(sheet{i,2}) == 0, strcat("Cell B", num2str(i), " in ", reliability_data_path, " should not be empty. ", mat2str(sheet{i,2})))
        for j = 4:5
            assert(isempty(sheet{i,j}) == 0, strcat("No cell in D", num2str(i), ":E", num2str(i), " in ", reliability_data_path, " should be empty. ", mat2str(sheet{i,4:5})))
        endfor
        comp_name = strrep (comp_name, " ", "_");
        input_data.(comp_name) = struct();
        input_data.(comp_name).MTTR = sheet{i,2};
        input_data.(comp_name).failurerate_0 = sheet{i,4};
        input_data.(comp_name).failurerate_perkm = sheet{i,5};
    endfor

    % DC components (with several poles)
    dc_components = {"Converter ACDC", "DC OHL", "DC cable"};
    dc_rows = [6, 8, 9];
    for i = 1:3
        row = dc_rows(i);
        comp_name = strrep(dc_components{i}, " ", "_");
        assert(isempty(sheet{row,6}) == 0, strcat("Cell F", num2str(row), " in ", reliability_data_path, " should not be empty."))
        correlation_poles.(comp_name) = sheet{row,6};
    endfor

    warning off

    % Length & type per component
    macro_data_path = [work_dir, '\user_interface\inputs\', macro_scenario, '_model.xlsx'];
    assert(isfile(macro_data_path) == 1, strcat(macro_data_path, " does not exist"))

    [~, ~, branch_sheet] = xlsread(macro_data_path, 'branch');
    assert(branch_sheet{1,5} == "type", strcat("Cell E1 in sheet branch in ", macro_data_path, " should be type, not ", branch_sheet{1,5}))
    assert(branch_sheet{1,6} == "length", strcat("Cell F1 in sheet branch in ", macro_data_path, " should be type, not ", branch_sheet{1,6}))

    [~, ~, branchdc_sheet] = xlsread(macro_data_path, 'branchdc');
    assert(branchdc_sheet{1,5} == "type", strcat("Cell E1 in sheet branch in ", macro_data_path, " should be type, not ", branchdc_sheet{1,5}))
    assert(branchdc_sheet{1,6} == "length", strcat("Cell F1 in sheet branch in ", macro_data_path, " should be type, not ", branchdc_sheet{1,6}))
    branchdc_lengths = branchdc_sheet{4:3+size(mpc.branchdc,1),6};  % km

    % Giving names to buses AC and DC
    for j = 1:size(mpc.bus,1)
        mpc.bus_name{j} = ['AC-' num2str(mpc.bus(j,1))];
    end
    for j = 1:size(mpc.busdc,1)
        mpc.busdc_name{j} = ['DC-' num2str(mpc.busdc(j,1))];
    end

    mpc0 = mpc;
    % Correspondence between wires N,R,P and DC branches id & type & length
    branchdc_id_per_pole = [];  % branchdc indexes (1 value per pole P/N/R)
    pole_ids=[];  % 1 value per pole per branchdc: 0=R, 1=P, 2=N
    branchdc_lengths_per_pole = [];
    branchdc_types = {};  % 'DC cable' or 'DC OHL'
    branchdc_type_per_pole = [];  % 1 for DC cable, 2 for DC OHL
    dc_type_cable = 1;
    dc_type_ohl = 2;
    for i = 1:size(mpc0.branchdc,1)
        dc_length = branchdc_sheet{3+i,6};
        branchdc_type = branchdc_sheet{3+i,5};
        branchdc_types{i} = branchdc_type;
        if strcmp(branchdc_type, "DC cable")
            dc_type = dc_type_cable;
        else
            assert(strcmp(branchdc_type, "DC OHL"), strcat("Type for DC branch ", i, " should be 'DC cable' or 'DC OHL', not ", branchdc_type))
            dc_type = dc_type_ohl;
        endif
        if mpc0.branchdc(i,10)==2  % Column 10 is line_confi (1=monopolar, 2=bipolar)
            for pole_id = 1:3
                branchdc_lengths_per_pole = [branchdc_lengths_per_pole, dc_length];
                branchdc_type_per_pole = [branchdc_type_per_pole, dc_type];
                branchdc_id_per_pole = [branchdc_id_per_pole i];
                pole_ids=[pole_ids pole_id-1];
            end
        else  % Monopolar branchdc
            used_poles = find(mpc.branchdc(i,14:16)==1);  % 14=status_p, 15=status_n, 16=status_r
            for pole_id = 1:length(used_poles)
                branchdc_lengths_per_pole = [branchdc_lengths_per_pole, dc_length];
                branchdc_type_per_pole = [branchdc_type_per_pole, dc_type];
                branchdc_id_per_pole = [branchdc_id_per_pole i];
                switch used_poles(pole_id)
                    case 1
                        pole_ids=[pole_ids 1];
                    case 2
                        pole_ids=[pole_ids 2];
                    case 3
                        pole_ids=[pole_ids 0];
                end
            end
        end
    end

    dc_cable_id_per_pole = find(branchdc_type_per_pole == dc_type_cable);
    dc_ohl_id_per_pole = find(branchdc_type_per_pole == dc_type_ohl);

    % Indexes per type of AC branch
    branch_types = [];  % 1 for AC OHL, 2 for Transformer
    ac_type_ohl = 1;
    ac_type_transformer = 2;
    ac_ohl_lengths = [];  % Length (in km) of each AC OHL

    for i=1:size(mpc.branch,1)
        branch_type = branch_sheet{3+i,5};
        if strcmp(branch_type, "AC OHL")
            ac_type = ac_type_ohl;
            ac_ohl_lengths = [ac_ohl_lengths, branch_sheet{3+i,6}];  % km
        else
            assert(strcmp(branch_type, "Transformer"), strcat("Type for AC branch ", i, " should be 'AC OHL' or 'Transformer', not ", branch_type))
            ac_type = ac_type_transformer;
        endif
        branch_types = [branch_types, ac_type];
    endfor
    ac_ohl_ids = find(branch_types == ac_type_ohl);
    transfo_ids = find(branch_types == ac_type_transformer);

    %%% --- Computation of MTTR & MTTF for each component ---

    % Generator
    MTTRsgen_ac = ones(1,size(mpc.gen,1)).*input_data.Generator.MTTR;
    MTTFsgen_ac = MTTRsgen_ac*(1/input_data.Generator.FOR - 1);

    % Transformer
    MTTRbrs_ac(transfo_ids) = input_data.Transformer.MTTR;
    MTTFbrs_ac(transfo_ids) = input_data.Transformer.MTTR*(1/input_data.Transformer.FOR - 1);

    % Converter
    MTTRsconv_dc = ones(1,size(mpc.convdc,1)).*input_data.Converter_ACDC.MTTR;
    MTTFsconv_dc = MTTRsconv_dc.*(1/input_data.Converter_ACDC.FOR - 1);

    % AC OHL
    MTTRbrs_ac(ac_ohl_ids) = input_data.AC_OHL.MTTR;
    MTTFbrs_ac(ac_ohl_ids) = 8760./(ac_ohl_lengths(ac_ohl_ids).*input_data.AC_OHL.failurerate_perkm + input_data.AC_OHL.failurerate_0);

    % DC OHL
    MTTRbrs_dc(dc_ohl_id_per_pole) = input_data.DC_OHL.MTTR;
    MTTFbrs_dc(dc_ohl_id_per_pole) = 8760./(branchdc_lengths_per_pole(dc_ohl_id_per_pole).*input_data.DC_OHL.failurerate_perkm + input_data.DC_OHL.failurerate_0);

    % DC cable
    MTTRbrs_dc(dc_cable_id_per_pole) = input_data.DC_cable.MTTR;
    MTTFbrs_dc(dc_cable_id_per_pole) = 8760./(branchdc_lengths_per_pole(dc_cable_id_per_pole).*input_data.DC_cable.failurerate_perkm + input_data.DC_cable.failurerate_0);

    % Data to return
    reliability_data.mpc = mpc;
    reliability_data.output_folder = output_folder;

    reliability_data.MTTRsgen_ac = MTTRsgen_ac;
    reliability_data.MTTFsgen_ac = MTTFsgen_ac;

    reliability_data.MTTFbrs_ac = MTTFbrs_ac;
    reliability_data.MTTRbrs_ac = MTTRbrs_ac;

    reliability_data.MTTRbrs_dc = MTTRbrs_dc;
    reliability_data.MTTFbrs_dc = MTTFbrs_dc;

    reliability_data.MTTRsconv_dc = MTTRsconv_dc;
    reliability_data.MTTFsconv_dc = MTTFsconv_dc;

    reliability_data.branchdc_types = branchdc_types;
    reliability_data.correlation_poles = correlation_poles;

    reliability_data.branchdc_id_per_pole = branchdc_id_per_pole;
    reliability_data.pole_ids = pole_ids;
end
