% Code developed by Andrea PITTO (RSE) and Nicolas BARLA (SuperGrid Institute)"
% Script for the generation of unavailabilities seriesfor N-1 outages including HVDC assets

% Uncomment (and update) the following 4 lines to run this script directly with Octave
% n_series = 2;
% n_hours = 3;
% work_dir = "c:\\Users\\n.barla\\Documents\\Local_codes\\HVDCWISE_TEA.jl\\studies\\study_2";
% previous_work_dir = "c:\\Users\\n.barla\\Documents\\Local_codes\\HVDCWISE_TEA.jl\\studies\\study_1";


clc
vars_to_delete = setdiff(who, {'n_hours', 'n_series', 'previous_work_dir', 'work_dir'});
if length(vars_to_delete) > 0
    clear(vars_to_delete{:});  % clear all variables except the expected ones
endif

pkg load statistics
format long

addpath(fileparts(mfilename('fullpath')));

%%% Retrieve input data
assert(exist('work_dir', 'var') == 1, "work_dir is not defined")
assert(exist('previous_work_dir', 'var') == 1, "previous_work_dir is not defined")
assert(exist('n_series', 'var') == 1, "n_series is not defined")
assert(exist('n_hours', 'var') == 1, "n_hours is not defined")
if ~isnumeric(n_series)
    n_series = str2num(n_series);
endif
if ~isnumeric(n_hours)
    n_hours = str2num(n_hours);
endif

reliability_data = read_reliability_data(work_dir);

mpc = reliability_data.mpc;
has_dc = reliability_data.has_dc;

folderout = reliability_data.output_folder;
confirm_recursive_rmdir(0);
if isfolder(folderout)
    rmdir(folderout, 's');  % 's' to also delete all subfolders and files
endif
mkdir(folderout);

MTTRsgen_ac = reliability_data.MTTRsgen_ac;
MTTFsgen_ac = reliability_data.MTTFsgen_ac;

MTTRbrs_ac = reliability_data.MTTRbrs_ac;
MTTFbrs_ac = reliability_data.MTTFbrs_ac;
branch_n_parallel = reliability_data.branch_n_parallel;  % Number of real lines per model (aggregated) line

MTTRbrs_dc = reliability_data.MTTRbrs_dc;
MTTFbrs_dc = reliability_data.MTTFbrs_dc;

MTTRsconv_dc = reliability_data.MTTRsconv_dc;
MTTFsconv_dc = reliability_data.MTTFsconv_dc;

correlation_poles_conv = reliability_data.correlation_poles.Converter_ACDC;
correlation_poles_dc_ohl = reliability_data.correlation_poles.DC_OHL;
correlation_poles_dc_cable = reliability_data.correlation_poles.DC_cable;
branchdc_types = reliability_data.branchdc_types;

branchdc_id_per_pole = reliability_data.branchdc_id_per_pole;  % branchdc indexes (1 value per pole P/N/R)
pole_ids = reliability_data.pole_ids;  % 1 value per pole per branchdc: 0=R, 1=P, 2=N

% Initialization
old_n_series = 0;
old_n_hours = 0;
old_n_branch = 0;
old_n_branchdc = 0;
old_n_gen = 0;
old_n_conv = 0;

if isfolder(previous_work_dir)
    % Read number of series and components per type

    old_availability_series_dir = [previous_work_dir, '\simulation_interface\Input_series\Availability'];

    first_folder = [old_availability_series_dir filesep 'YR_1'];
    assert(isfolder(first_folder), strcat('Folder YR_1 is missing in ', old_availability_series_dir))

    old_n_hours = size(csvread([first_folder, '\branch\br_status.csv']), 1) - 1;  % First row is for component ids

    old_n_branch = size(csvread([first_folder, '\branch\br_status.csv']), 2);
    old_n_branchdc = size(csvread([first_folder, '\branchdc\status_p.csv']), 2);
    old_n_gen = size(csvread([first_folder, '\gen\gen_status.csv']), 2);
    old_n_conv = size(csvread([first_folder, '\convdc\status_p.csv']), 2);

    % folder_list = ls(old_availability_series_dir);
    old_folders = dir(old_availability_series_dir);
    old_folders(1:2) = [];
    old_n_series = length(old_folders);

    % old_n_series = size(folder_list, 1);
    for s=1:old_n_series
        expected_folder_name = strcat("YR_", num2str(s));
        assert(
            any(strcmp({old_folders.name}, expected_folder_name)),
            strcat('There are ', old_n_series, ' folders in ', old_availability_series_dir, '. Folder ', expected_folder_name, ' is missing.', old_folders.name)
            )
    endfor

else
    assert(length(previous_work_dir) == 0, strcat('The provided previous_work_dir does not exist ', previous_work_dir))
    disp('No previous work dir provided')
endif

n_h = min(n_hours, old_n_hours);
n_s = min(n_series, old_n_series);

%%% START - GENERATION of unavailabilities (considering normal events)


% AC  branches
clear br_status br_rating

n_branch = size(mpc.branch,1);
br_status = ones(n_hours*n_series, n_branch);
br_rating = ones(n_hours*n_series, n_branch);
br_resistance = ones(n_hours*n_series, n_branch);
for branch_id = 1:n_branch  % branch id
    n_parallel = branch_n_parallel(branch_id);
    detailed_status = ones(n_hours*n_series, n_parallel);
    for i = 1:n_parallel
        t_ref = 1;
        while t_ref <= n_hours*n_series
            r = rand;  % Random number between 0 & 1
            if detailed_status(t_ref, i) == 1  % Available
                TTF = ceil(-MTTFbrs_ac(branch_id)*log(r));  % Time to failure. Probability density function (TTF=t) = exp(-t/MTTF)
                t1 = min(t_ref + TTF - 1, n_hours*n_series);  % Availabiliy ends
                detailed_status(t_ref:t1, i) = 1;
            else
                TTR = ceil(-MTTRbrs_ac(branch_id)*log(r));  % Time to repair
                t1 = min(t_ref + TTR - 1, n_hours*n_series);  % Unavailabiliy ends
                detailed_status(t_ref:t1, i) = 0;
            end
            if t1 < n_hours*n_series
                detailed_status(t1 + 1, i) = 1 - detailed_status(t_ref, i);
            endif
            t_ref = t1 + 1;  % (Un)Availabiliy starts
        end
    endfor
    br_status(1:n_hours*n_series,branch_id) = (sum(detailed_status, 2) > 0);
    br_rating(1:n_hours*n_series,branch_id) = mpc.branch(branch_id,6) * sum(detailed_status, 2) / n_parallel;
    br_resistance(1:n_hours*n_series,branch_id) = mpc.branch(branch_id,3) * n_parallel ./ max(1, sum(detailed_status, 2));
    br_reactance(1:n_hours*n_series,branch_id) = mpc.branch(branch_id,4) * n_parallel ./ max(1, sum(detailed_status, 2));
end
for s = 1:n_series
    branch_status{s} = br_status(1+(s-1)*n_hours:n_hours*s, :);
    branch_rate_a{s} = br_rating(1+(s-1)*n_hours:n_hours*s, :);
    branch_resistance{s} = br_resistance(1+(s-1)*n_hours:n_hours*s, :);
    branch_reactance{s} = br_reactance(1+(s-1)*n_hours:n_hours*s, :);
end
n_comp = min(n_branch, old_n_branch);
if n_s * n_h * n_comp > 0
    for s = 1:n_s
        old_folder = [old_availability_series_dir filesep strcat('YR_', num2str(s))];
        % TODO assert coherent n_parallel in previous_work_dir

        old_status_file = [old_folder, '\branch\br_status.csv'];
        old_status_data = csvread(old_status_file);
        assert(old_status_data(1, 1:n_comp) == 1:n_comp, strcat('Component ids in ', old_status_file, ' should be ', num2str(1:n_comp), ' not ', num2str(old_status_data(1, 1:n_comp))));
        branch_status{s}(1:n_h, 1:n_comp) = old_status_data(2:n_h+1, 1:n_comp);  % First row is for component ids

        old_rating_file = [old_folder, '\branch\rate_a.csv'];
        old_rating_data = csvread(old_rating_file);
        assert(old_rating_data(1, 1:n_comp) == 1:n_comp, strcat('Component ids in ', old_rating_file, ' should be ', num2str(1:n_comp), ' not ', num2str(old_rating_data(1, 1:n_comp))));
        branch_rate_a{s}(1:n_h, 1:n_comp) = old_rating_data(2:n_h+1, 1:n_comp);  % First row is for component ids

        old_resistance_file = [old_folder, '\branch\br_r.csv'];
        old_resistance_data = csvread(old_resistance_file);
        assert(old_resistance_data(1, 1:n_comp) == 1:n_comp, strcat('Component ids in ', old_resistance_file, ' should be ', num2str(1:n_comp), ' not ', num2str(old_resistance_data(1, 1:n_comp))));
        branch_resistance{s}(1:n_h, 1:n_comp) = old_resistance_data(2:n_h+1, 1:n_comp);  % First row is for component ids

        old_reactance_file = [old_folder, '\branch\br_x.csv'];
        old_reactance_data = csvread(old_reactance_file);
        assert(old_reactance_data(1, 1:n_comp) == 1:n_comp, strcat('Component ids in ', old_reactance_file, ' should be ', num2str(1:n_comp), ' not ', num2str(old_reactance_data(1, 1:n_comp))));
        branch_reactance{s}(1:n_h, 1:n_comp) = old_reactance_data(2:n_h+1, 1:n_comp);  % First row is for component ids
    endfor
endif


% AC Gens
clear Uacgen

Uacgen = ones(n_hours*n_series, size(mpc.gen,1));
n_gen = size(mpc.gen,1);
for gen_id = 1:n_gen
    t_ref = 1;
    while t_ref <= n_hours*n_series
        r = rand;  % Random number between 0 & 1
        if Uacgen(t_ref, gen_id) == 1  % Available
            TTF = ceil(-MTTFsgen_ac(gen_id)*log(r));  % Time to failure. Probability density function (TTF=t) = exp(-t/MTTF)
            t1 = min(t_ref + TTF - 1, n_hours*n_series);  % Availabiliy ends
            Uacgen(t_ref:t1, gen_id) = 1;
        else
            TTR = ceil(-MTTRsgen_ac(gen_id)*log(r));  % Time to repair
            t1 = min(t_ref + TTR - 1, n_hours*n_series);  % Unavailabiliy ends
            Uacgen(t_ref:t1, gen_id) = 0;
        end
        if t1 < n_hours*n_series
            Uacgen(t1 + 1, gen_id) = 1 - Uacgen(t_ref, gen_id);
        endif
        t_ref = t1 + 1;  % (Un)Availabiliy starts
    end
end
for i = 1:n_series
    Uacgens{i} = Uacgen(1+(i-1)*n_hours:n_hours*i, :);
end
n_comp = min(n_gen, old_n_gen);
if n_s * n_h * n_comp > 0
    for s = 1:n_s
        old_folder = [old_availability_series_dir filesep strcat('YR_', num2str(s))];
        old_file = [old_folder, '\gen\gen_status.csv'];
        old_data = csvread(old_file);
        assert(old_data(1, 1:n_comp) == 1:n_comp, strcat('Component ids in', old_file, ' should be ', num2str(1:n_comp), ' not ', num2str(old_data(1, 1:n_comp))));
        Uacgens{s}(1:n_h, 1:n_comp) = old_data(2:n_h+1, 1:n_comp);
    endfor
endif


if has_dc
    A=[]; BV=[];

    % DC branches

    n_branchdc = size(mpc.branchdc,1);
    Udcline_r = ones(n_hours*n_series, n_branchdc);
    Udcline_p = ones(n_hours*n_series, n_branchdc);
    Udcline_n = ones(n_hours*n_series, n_branchdc);

    for branchdc_id = 1:n_branchdc
        branchdc_type = branchdc_types{branchdc_id};
        if strcmp(branchdc_type, 'DC cable')
            correlation = correlation_poles_dc_cable;
        else
            assert(strcmp(branchdc_type, 'DC OHL'), strcat('branchdc_type should be "DC cable" or "DC OHL", not ', branchdc_type))
            correlation = correlation_poles_dc_ohl;
        endif
        idxs = find(branchdc_id_per_pole == branchdc_id);
        used_poles = pole_ids(idxs);
        t_ref = 1;
        while t_ref <= n_hours*n_series
            r = rand;

            for ipol = 1:length(used_poles)
                if used_poles(ipol) == 0  % Metallic return
                    Status_ini(ipol) = Udcline_r(t_ref, branchdc_id);  % Initial status
                elseif used_poles(ipol) == 1  % Pole P
                    Status_ini(ipol) = Udcline_p(t_ref, branchdc_id);
                else  % Pole N
                    Status_ini(ipol) = Udcline_n(t_ref, branchdc_id);
                endif

                if  Status_ini(ipol)  % Pole is available
                    transition_rate(ipol) = 1/MTTFbrs_dc(idxs(ipol));  % Failure rate (per hour)
                    trans_tipo(ipol)=0;  % Transition type (0 = to failure)
                else  % Pole is unavailable
                    transition_rate(ipol) = 1/MTTRbrs_dc(idxs(ipol));  % Repair rate (per hour)
                    trans_tipo(ipol) = 1;  % Transition type (1 = to repair)
                end
            endfor

            if correlation > 0  % Correlation between pole failures

                transitio = [];
                Vs = [];
                for k = 0:length(used_poles)
                    v = nchoosek([1:length(used_poles)], k);
                    for iv = 1:size(v,1)
                        Vs = [Vs; Status_ini];
                        Vs(end,v(iv,:))=not((Status_ini(v(iv,:))));

                        a = [transition_rate];
                        b = [ones(1,length(transition_rate))];
                        AB = [a;b];
                        dummy = 0;
                        quali_compl = (Vs(end,:)-Status_ini==0)*2+(abs(Vs(end,:)-Status_ini)>0)*1;

                        correlation_matrix = eye(length(used_poles));  % Correlation between poles
                        idxf = find(Status_ini==1);
                        for i1 = 1:length(idxf)-1
                            for i2 = i1+1:length(idxf)
                                correlation_matrix(idxf(i1),idxf(i2)) = correlation;
                                correlation_matrix(idxf(i2),idxf(i1)) = correlation;
                            end
                        end

                        for in1 = 1:quali_compl(1)
                            for in2 = 1:quali_compl(2)
                                if length(used_poles) == 3
                                    for in3 = 1:quali_compl(3)
                                        Vx = [in1 in2 in3];
                                        Vss = (abs(Vs(end,:)-Status_ini)==0)*2+(abs(Vs(end,:)-Status_ini)>0)*1;
                                        dummy = dummy + ((-1).^(sum(abs(Vss-Vx)))).*copulacdf('gaussian',[AB(in1,1) AB(in2,2) AB(in3,3)], correlation_matrix);
                                    end
                                else
                                    Vx = [in1 in2 ];
                                    Vss = (abs(Vs(end,:)-Status_ini)==0)*2+(abs(Vs(end,:)-Status_ini)>0)*1;
                                    dummy = dummy + ((-1).^(sum(abs(Vss-Vx)))).*copulacdf('gaussian',[AB(in1,1) AB(in2,2)], correlation_matrix);
                                end
                            end
                        end
                        transitio = [transitio; dummy];
                    end
                end
                iqualeini = find(ismember(Vs,Status_ini,'rows'));
                transitio(iqualeini) = [];
                Vs(iqualeini,:)=[];

                TT = ceil(-1/sum(transitio)*log(r));  % Transition time
                rd = rand;

                transition_ = transitio./sum(transitio);

                idx = find(cumsum(transition_) > rd);

                qualeTran = idx(1);

                dquali = Vs(qualeTran,:) - Status_ini;

                dummytrans0 = find(not(dquali == 0));

                for idum = 1:length(dummytrans0)
                    ipol = dummytrans0(idum);
                    t1 = min(t_ref + TT - 1, n_hours*n_series);  % (Un)Availabiliy ends
                    switch used_poles(ipol)
                        case 0
                            if trans_tipo(ipol) == 0  % Transition type (0 = to failure)
                                Udcline_r(t_ref:t1, branchdc_id) = 1;
                            else
                                Udcline_r(t_ref:t1, branchdc_id) = 0;
                            end
                            if t1 < n_hours*n_series
                                Udcline_r(t1 + 1, branchdc_id) = 1 - Udcline_r(t_ref, branchdc_id);
                            endif
                        case 1
                            if trans_tipo(ipol) == 0
                                Udcline_p(t_ref:t1, branchdc_id) = 1;
                            else
                                Udcline_p(t_ref:t1, branchdc_id) = 0;
                            end
                            if t1 < n_hours*n_series
                                Udcline_p(t1 + 1, branchdc_id) = 1 - Udcline_p(t_ref, branchdc_id);
                            endif
                        case 2
                            if trans_tipo(ipol) == 0
                                Udcline_n(t_ref:t1, branchdc_id) = 1;
                            else
                                Udcline_n(t_ref:t1, branchdc_id) = 0;
                            end
                            if t1 < n_hours*n_series
                                Udcline_n(t1 + 1, branchdc_id) = 1 - Udcline_n(t_ref, branchdc_id);
                            endif
                    end
                    t_ref = t1 + 1;  % (Un)Availability starts
                end
                clear transition_rate Status_ini Vs transitio

            else  % No correlation between pole failures

                TT = ceil(-1/sum(transition_rate)*log(r));
                rd = rand;

                transition_rate = transition_rate./sum(transition_rate);

                idx = find(cumsum(transition_rate) > rd);

                ipol = idx(1);

                t1 = min(t_ref + TT - 1, n_hours*n_series);  % Availabiliy ends
                switch used_poles(ipol)
                    case 0
                        if trans_tipo(ipol) == 0  % Transition type (0 = to failure)
                            Udcline_r(t_ref:t1, branchdc_id) = 1;
                        else
                            Udcline_r(t_ref:t1, branchdc_id) = 0;
                        end
                        if t1 < n_hours*n_series
                            Udcline_r(t1 + 1, branchdc_id) = 1 - Udcline_r(t_ref, branchdc_id);
                        endif
                    case 1
                        if trans_tipo(ipol) == 0
                            Udcline_p(t_ref:t1, branchdc_id) = 1;
                        else
                            Udcline_p(t_ref:t1, branchdc_id) = 0;
                        end
                        if t1 < n_hours*n_series
                            Udcline_p(t1 + 1, branchdc_id) = 1 - Udcline_p(t_ref, branchdc_id);
                        endif
                    case 2
                        if trans_tipo(ipol) == 0
                            Udcline_n(t_ref:t1, branchdc_id) = 1;
                        else
                            Udcline_n(t_ref:t1, branchdc_id) = 0;
                        end
                        if t1 < n_hours*n_series
                            Udcline_n(t1 + 1, branchdc_id) = 1 - Udcline_n(t_ref, branchdc_id);
                        endif
                end
                t_ref = t1 + 1;  % Unavailability starts

                clear transition_rate Status_ini
            endif
        end
    endfor

    for s = 1:n_series
        Udclines_r{s} = Udcline_r(1+(s-1)*n_hours:n_hours*s, :);
        Udclines_p{s} = Udcline_p(1+(s-1)*n_hours:n_hours*s, :);
        Udclines_n{s} = Udcline_n(1+(s-1)*n_hours:n_hours*s, :);
    end
    n_comp = min(n_branchdc, old_n_branchdc);
    if n_s * n_h * n_comp > 0
        for s = 1:n_s
            old_folder = [old_availability_series_dir filesep strcat('YR_', num2str(s))];
            old_file_r = [old_folder, '\branchdc\status_r.csv'];
            old_file_p = [old_folder, '\branchdc\status_p.csv'];
            old_file_n = [old_folder, '\branchdc\status_n.csv'];
            old_data_r = csvread(old_file_r);
            old_data_p = csvread(old_file_p);
            old_data_n = csvread(old_file_n);
            assert(old_data_r(1, 1:n_comp) == 1:n_comp, strcat('Component ids in', old_file_r, ' should be ', num2str(1:n_comp), ' not ', num2str(old_data_r(1, 1:n_comp))));
            assert(old_data_p(1, 1:n_comp) == 1:n_comp, strcat('Component ids in', old_file_p, ' should be ', num2str(1:n_comp), ' not ', num2str(old_data_p(1, 1:n_comp))));
            assert(old_data_n(1, 1:n_comp) == 1:n_comp, strcat('Component ids in', old_file_n, ' should be ', num2str(1:n_comp), ' not ', num2str(old_data_n(1, 1:n_comp))));
            Udclines_r{s}(1:n_h, 1:n_comp) = old_data_r(2:n_h+1, 1:n_comp);
            Udclines_p{s}(1:n_h, 1:n_comp) = old_data_p(2:n_h+1, 1:n_comp);
            Udclines_n{s}(1:n_h, 1:n_comp) = old_data_n(2:n_h+1, 1:n_comp);
        endfor
    endif


    % DC converters
    clear Udcconv_p Udcconv_n

    n_conv = size(mpc.convdc,1);
    Udcconv_p = ones(n_hours*n_series, n_conv);
    Udcconv_n = ones(n_hours*n_series, n_conv);
    if correlation_poles_conv == 0  % No correlation between pole failures
        % Positive poles
        for conv_id = 1:n_conv
            t_ref = 1;
            while t_ref <= n_hours*n_series
                r = rand;
                if Udcconv_p(t_ref, conv_id) == 1
                    TTF = ceil(-MTTFsconv_dc(conv_id)*log(r));  % Time to failure. Probability density function (TTF=t) = exp(-t/MTTF)
                    t1 = min(t_ref + TTF - 1, n_hours*n_series);  % Availabiliy ends
                    Udcconv_p(t_ref:t1, conv_id) = 1;
                else
                    TTR = ceil(-MTTRsconv_dc(conv_id)*log(r));  % Time to repair
                    t1 = min(t_ref + TTR - 1, n_hours*n_series);  % Unavailabiliy ends
                    Udcconv_p(t_ref:t1, conv_id) = 0;
                end
                if t1 < n_hours*n_series
                    Udcconv_p(t1 + 1, conv_id) = 1 - Udcconv_p(t_ref, conv_id);
                endif
                t_ref = t1 + 1;  % (Un)Availabiliy starts
            end
        end
        % Negative poles
        for conv_id = 1:n_conv
            t_ref = 1;
            while t_ref <= n_hours*n_series
                r = rand;
                if Udcconv_n(t_ref, conv_id) == 1
                    TTF = ceil(-MTTFsconv_dc(conv_id)*log(r));  % Time to failure. Probability density function (TTF=t) = exp(-t/MTTF)
                    t1 = min(t_ref + TTF - 1, n_hours*n_series);  % Availabiliy ends
                    Udcconv_n(t_ref:t1, conv_id) = 1;
                else
                    TTR = ceil(-MTTRsconv_dc(conv_id)*log(r));  % Time to repair
                    t1 = min(t_ref + TTR - 1, n_hours*n_series);  % Unavailabiliy ends
                    Udcconv_n(t_ref:t1, conv_id) = 0;
                end
                if t1 < n_hours*n_series
                    Udcconv_n(t1 + 1, conv_id) = 1 - Udcconv_n(t_ref, conv_id);
                endif
                t_ref = t1 + 1;  % (Un)Availabiliy starts
            end
        end

    else  % Correlation between pole failures
        for conv_id = 1:n_conv
            t_ref = 1;
            while t_ref <= n_hours*n_series
                used_poles = [1 2];
                r = rand;

                for pole_id = 1:length(used_poles)
                    if used_poles(pole_id) == 1  % Pole P
                        StatusC_ini(pole_id) = Udcconv_p(t_ref, conv_id);
                    else  % Pole N
                        StatusC_ini(pole_id) = Udcconv_n(t_ref, conv_id);
                    endif

                    if  StatusC_ini(pole_id)  % Pole is available
                        transition_rate(pole_id) = 1/MTTFsconv_dc(conv_id);  % Failure rate (per hour)
                        trans_tipo(pole_id)=0;  % Transition type (0 = to failure)
                    else  % Pole is unavailable
                        transition_rate(pole_id) = 1/MTTRsconv_dc(conv_id);  % Repair rate (per hour)
                        trans_tipo(pole_id) = 1;  % Transition type (1 = to repair)
                    end
                endfor

                transitio = [];
                Vs = [];

                for k = 0:length(used_poles)
                    v = nchoosek([1:length(used_poles)], k);

                    for iv = 1:size(v,1)
                        Vs = [Vs; StatusC_ini];
                        Vs(end,v(iv,:))=not((StatusC_ini(v(iv,:))));

                        a = [transition_rate];
                        b = [ones(1,length(transition_rate))];
                        AB = [a;b];
                        dummy = 0;
                        quali_compl = (Vs(end,:)-StatusC_ini==0)*2+(abs(Vs(end,:)-StatusC_ini)>0)*1;

                        correlation_matrix = eye(length(used_poles) );
                        idxf = find(StatusC_ini==1);
                        for i1 = 1:length(idxf)-1
                            for i2 = i1+1:length(idxf)
                                correlation_matrix(idxf(i1),idxf(i2)) = correlation_poles_conv;
                                correlation_matrix(idxf(i2),idxf(i1)) = correlation_poles_conv;
                            end
                        end
                        for in1 = 1:quali_compl(1)
                            for in2 = 1:quali_compl(2)
                                Vx = [in1 in2 ];
                                Vss = (abs(Vs(end,:)-StatusC_ini)==0)*2+(abs(Vs(end,:)-StatusC_ini)>0)*1;
                                dummy = dummy + ((-1).^(sum(abs(Vss-Vx)))).*copulacdf('gaussian',[AB(in1,1) AB(in2,2)],correlation_matrix);
                            end
                        end
                        transitio = [transitio; dummy];
                    end
                end
                iqualeini = find(ismember(Vs,StatusC_ini,'rows'));
                transitio(iqualeini) = [];
                Vs(iqualeini,:)=[];

                TT = ceil(-1/sum(transitio)*log(r));  % Transition time
                rd = rand;

                transition_ = transitio./sum(transitio);

                idx = find(cumsum(transition_) > rd);

                qualeTran = idx(1);

                dquali = Vs(qualeTran,:) - StatusC_ini;

                dummytrans0 = find(not(dquali == 0));

                for idum = 1:length(dummytrans0)
                    pole_id = dummytrans0(idum);
                    t1 = min(t_ref + TT - 1, n_hours*n_series);  % (Un)Availabiliy ends
                    switch used_poles(pole_id)
                        case 1
                            if trans_tipo(pole_id) == 0
                                Udcconv_p(t_ref:t1, conv_id) = 1;
                            else
                                Udcconv_p(t_ref:t1, conv_id) = 0;
                            end
                            if t1 < n_hours*n_series
                                Udcconv_p(t1 + 1, conv_id) = 1 - Udcconv_p(t_ref, conv_id);
                            endif
                        case 2
                            if trans_tipo(pole_id) == 0
                                Udcconv_n(t_ref:t1, conv_id) = 1;
                            else
                                Udcconv_n(t_ref:t1, conv_id) = 0;
                            end
                            if t1 < n_hours*n_series
                                Udcconv_n(t1 + 1, conv_id) = 1 - Udcconv_n(t_ref, conv_id);
                            endif
                    end
                    if t1 < n_hours*n_series
                        Udcconv_p(t1 + 1, conv_id) = 1 - Udcconv_p(t_ref, conv_id);
                        Udcconv_n(t1 + 1, conv_id) = 1 - Udcconv_n(t_ref, conv_id);
                    endif
                    t_ref = t1 + 1;  % (Un)Availability starts
                end
                clear transition_rate StatusC_ini Vs transitio
            end
        end
    end
    for s = 1:n_series
        Udcconvs_n{s} = Udcconv_n(1+(s-1)*n_hours:n_hours*s, :);
        Udcconvs_p{s} = Udcconv_p(1+(s-1)*n_hours:n_hours*s, :);
    end
    n_comp = min(n_conv, old_n_conv);
    if n_s * n_h * n_comp > 0
        for s = 1:n_s
            old_folder = [old_availability_series_dir filesep strcat('YR_', num2str(s))];
            old_file_p = [old_folder, '\convdc\status_p.csv'];
            old_file_n = [old_folder, '\convdc\status_n.csv'];
            old_data_p = csvread(old_file_p);
            old_data_n = csvread(old_file_n);
            assert(old_data_p(1, 1:n_comp) == 1:n_comp, strcat('Component ids in', old_file_p, ' should be ', num2str(1:n_comp), ' not ', num2str(old_data_p(1, 1:n_comp))));
            assert(old_data_n(1, 1:n_comp) == 1:n_comp, strcat('Component ids in', old_file_n, ' should be ', num2str(1:n_comp), ' not ', num2str(old_data_n(1, 1:n_comp))));
            Udclines_p{s}(1:n_h, 1:n_comp) = old_data_p(2:n_h+1, 1:n_comp);
            Udclines_n{s}(1:n_h, 1:n_comp) = old_data_n(2:n_h+1, 1:n_comp);
        endfor
    endif
endif

%%% PRINT CSV FILES FOR CSV TREE STRUCTURE

for i = 1:n_series
    folderctgout = [folderout filesep 'YR_' num2str(i)];

    mkdir(folderctgout);
    mkdir([folderctgout filesep 'branch'])
    csvwrite([folderctgout filesep 'branch' filesep 'br_status.csv'],[[1:size(mpc.branch,1)];branch_status{i}]);
    csvwrite([folderctgout filesep 'branch' filesep 'rate_a.csv'],[[1:size(mpc.branch,1)];branch_rate_a{i}]);
    csvwrite([folderctgout filesep 'branch' filesep 'br_r.csv'],[[1:size(mpc.branch,1)];branch_resistance{i}]);
    csvwrite([folderctgout filesep 'branch' filesep 'br_x.csv'],[[1:size(mpc.branch,1)];branch_reactance{i}]);
    mkdir([folderctgout filesep 'gen'])
    csvwrite([folderctgout filesep 'gen' filesep 'gen_status.csv'],[[1:size(mpc.gen,1)];Uacgens{i}]);
    if has_dc
        mkdir([folderctgout filesep 'branchdc'])
        Udcbr_r{i}=Udclines_r{i}(1:n_hours, branchdc_id_per_pole(find(pole_ids==0)));
        csvwrite([folderctgout filesep 'branchdc' filesep 'status_r.csv'],[[branchdc_id_per_pole(find(pole_ids==0))];Udcbr_r{i}])
        Udcbr_p{i}=Udclines_p{i}(1:n_hours, branchdc_id_per_pole(find(pole_ids==1)));
        csvwrite([folderctgout filesep 'branchdc' filesep 'status_p.csv'],[[branchdc_id_per_pole(find(pole_ids==1))];Udcbr_p{i}])
        Udcbr_n{i}=Udclines_n{i}(1:n_hours, branchdc_id_per_pole(find(pole_ids==2)));
        csvwrite([folderctgout filesep 'branchdc' filesep 'status_n.csv'],[[branchdc_id_per_pole(find(pole_ids==2))];Udcbr_n{i}])
        mkdir([folderctgout filesep 'convdc'])
        csvwrite([folderctgout filesep 'convdc' filesep 'status_p.csv'],[[1:size(mpc.convdc,1)];Udcconvs_p{i}]);
        csvwrite([folderctgout filesep 'convdc' filesep 'status_n.csv'],[[1:size(mpc.convdc,1)];Udcconvs_n{i}]);
    endif
endfor

warning on
disp('Contingencies time series generation completed')
