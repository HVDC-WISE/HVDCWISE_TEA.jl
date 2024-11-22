% script for the generation of unavailabilities seriesfor N-1 outages including HVDC assets
clc
close all
if length(who) > 2
    clear(setdiff(who, {'n_hours', 'work_dir'}){:});  % clear all variables except work_dir and n_hours
endif

format long

addpath(fileparts(mfilename('fullpath')));

%%% Retrieve input data
assert(exist('work_dir', 'var') == 1, "work_dir is not defined")
assert(exist('n_hours', 'var') == 1, "n_hours is not defined")
if ~isnumeric(n_hours)
    n_hours = str2num(n_hours);
endif

# read_reliability_data;
reliability_data = read_reliability_data(work_dir);

mpc = reliability_data.mpc;
n_series = reliability_data.n_series;

folderout = reliability_data.output_folder;
if isfolder(folderout)
    rmdir(folderout, 's');  % 's' to also delete all subfolders and files
endif
mkdir(folderout);

MTTRsgen_ac = reliability_data.MTTRsgen_ac;
MTTFsgen_ac = reliability_data.MTTFsgen_ac;
MTTFbrs_ac = reliability_data.MTTFbrs_ac;
MTTRbrs_ac = reliability_data.MTTRbrs_ac;
MTTRbrs_dc = reliability_data.MTTRbrs_dc;
MTTFbrs_dc = reliability_data.MTTFbrs_dc;
MTTRsconv_dc = reliability_data.MTTRsconv_dc;
MTTFsconv_dc = reliability_data.MTTFsconv_dc;

corr_dcpoles_adequacy = reliability_data.corr_dcpoles_adequacy;
DCdependent_adequacy = reliability_data.DCdependent_adequacy;
corr_dcpoles_resilience = reliability_data.corr_dcpoles_resilience;

idexdc = reliability_data.idexdc;
poli = reliability_data.poli;

%%% START - GENERATION of N-1 and HVDC unavailabilities

A=[];BV=[];
clear Uacline Uacgen  Udcconv_p Udcconv_n

Udcline_r = ones(size(mpc.branchdc,1),n_hours*n_series);
Udcline_p = ones(size(mpc.branchdc,1),n_hours*n_series);
Udcline_n = ones(size(mpc.branchdc,1),n_hours*n_series);
Uacline = ones(size(mpc.branch,1),n_hours*n_series);
Uacgen = ones(size(mpc.gen,1),n_hours*n_series);
Udcconv_p = ones(size(mpc.convdc,1),n_hours*n_series);
Udcconv_n = ones(size(mpc.convdc,1),n_hours*n_series);

% AC  branches
for ibr = 1:size(mpc.branch,1)
    iyr=1;dt=1;in=1;
    while in < n_series
        r=rand;
        if Uacline(ibr,iyr)==1
            dt=ceil(-MTTFbrs_ac(ibr)*log(r));
            t1r = iyr+dt-1;
            t2r = iyr+dt;
            % iyr = t2r;
            if t2r< n_hours*n_series
                Uacline(ibr,iyr:t1r)=1;
                Uacline(ibr,t2r)=0;
                in = floor(t2r/n_hours);
                iyr=t2r;
            else
                Uacline(ibr,iyr:n_hours*n_series)=1;
                in=n_series;
            end
        else
            dt=ceil(-MTTRbrs_ac(ibr)*log(r));
            t1r = iyr+dt-1;
            t2r = iyr+dt;
            % iyr = t2r;
            if t2r< n_hours*n_series
                Uacline(ibr,iyr:t1r)=0;
                Uacline(ibr,t2r)=1;
                in =floor(t2r/n_hours);
                iyr=t2r;
            else
                Uacline(ibr,iyr:n_hours*n_series)=0;
                in=n_series;
            end
        end
    end
end
for i = 1:n_series
    Uaclines{i} = Uacline(:,1+(i-1)*n_hours:n_hours*i);
end
% DC branches
for dcbr = 1:size(mpc.branchdc,1)
    quali_poli = poli(find(idexdc==dcbr));
    idxs = find(idexdc==dcbr);
    iyr=1;dt=1;in=1;
    while in < n_series
        if DCdependent_adequacy == 1
            r = rand;
            for ipol = 1:length(quali_poli)
                if quali_poli(ipol)==0
                    try
                    Status_ini(ipol) = Udcline_r(dcbr,iyr) ;
                    catch err
                        keyboard
                    end
                    if  Status_ini(ipol)
                        transition_rate(ipol) = 1/MTTFbrs_dc(idxs(ipol));
                        trans_tipo(ipol)=0;
                    else
                        transition_rate(ipol) = 1/MTTRbrs_dc(idxs(ipol));
                        trans_tipo(ipol)=1;
                    end

                elseif quali_poli(ipol)==1
                    try
                        Status_ini(ipol) = Udcline_p(dcbr,iyr) ;
                    catch err
                        keyboard
                    end
                    if Status_ini(ipol)
                        transition_rate(ipol) = 1/MTTFbrs_dc(idxs(ipol));
                        trans_tipo(ipol)=0;
                    else
                        transition_rate(ipol) = 1/MTTRbrs_dc(idxs(ipol));
                        trans_tipo(ipol)=1;
                    end
                else
                    Status_ini(ipol) = Udcline_n(dcbr,iyr) ;
                    if Status_ini(ipol)
                        transition_rate(ipol) = 1/MTTFbrs_dc(idxs(ipol));
                        trans_tipo(ipol)=0;
                    else
                        transition_rate(ipol) = 1/MTTRbrs_dc(idxs(ipol));
                        trans_tipo(ipol)=1;
                    end
                end
            end
            transitio = [];Vs = [];
            for k = 0:length(quali_poli)
                v = nchoosek([1:length(quali_poli)], k);
                for iv = 1:size(v,1)
                    Vs = [Vs; Status_ini];
                    Vs(end,v(iv,:))=not((Status_ini(v(iv,:))));

                    a = [transition_rate];
                    b = [ones(1,length(transition_rate))];
                    AB = [a;b];
                    dummy = 0;
                    quali_compl = (Vs(end,:)-Status_ini==0)*2+(abs(Vs(end,:)-Status_ini)>0)*1;
                    if length(quali_poli) ==3
                        CORRE = eye(length(quali_poli) );
                        idxf = find(Status_ini==1);
                        for i1 = 1:length(idxf)-1
                            for i2 = i1+1:length(idxf)
                                CORRE(idxf(i1),idxf(i2)) = corr_dcpoles_adequacy;
                                CORRE(idxf(i2),idxf(i1)) = corr_dcpoles_adequacy;
                            end
                        end

                        for in1 = 1:quali_compl(1)
                            for in2 = 1:quali_compl(2)
                                for in3 = 1:quali_compl(3)

                                    Vx = [in1 in2 in3];
                                    Vss = (abs(Vs(end,:)-Status_ini)==0)*2+(abs(Vs(end,:)-Status_ini)>0)*1;

                                    dummy = dummy + ((-1).^(sum(abs(Vss-Vx)))).*copulacdf('gaussian',[AB(in1,1) AB(in2,2) AB(in3,3)],CORRE);
                                end
                            end
                        end
                    else
                        CORRE = eye(length(quali_poli) );
                        CORRE = eye(length(quali_poli) );
                        idxf = find(Status_ini==1);
                        for i1 = 1:length(idxf)-1
                            for i2 = i1+1:length(idxf)
                                CORRE(idxf(i1),idxf(i2)) = corr_dcpoles_adequacy;
                                CORRE(idxf(i2),idxf(i1)) = corr_dcpoles_adequacy;
                            end
                        end
                        for in1 = 1:quali_compl(1)
                            for in2 = 1:quali_compl(2)
                                Vx = [in1 in2 ];
                                Vss = (abs(Vs(end,:)-Status_ini)==0)*2+(abs(Vs(end,:)-Status_ini)>0)*1;
                                dummy = dummy + ((-1).^(sum(abs(Vss-Vx)))).*copulacdf('gaussian',[AB(in1,1) AB(in2,2)],CORRE);
                            end
                        end
                    end
                    transitio = [transitio; dummy];
                end
            end
            iqualeini = find(ismember(Vs,Status_ini,'rows'));
            transitio(iqualeini) = [];
            Vs(iqualeini,:)=[];

            dt=ceil(-1/sum(transitio)*log(1-r));
            rd = rand;

            transition_ = transitio./sum(transitio);

            idx = find(cumsum(transition_) > rd);

            qualeTran = idx(1);

            dquali = Vs(qualeTran,:) - Status_ini;

            t1r = iyr+dt-1;
            t2r = t1r+1;

            dummytrans0 = find(not(dquali == 0));

            for idum = 1:length(dummytrans0)
                ipol = dummytrans0(idum);
                switch quali_poli(ipol)
                    case 0
                        if trans_tipo(ipol)==0

                            if t2r< n_hours*n_series
                                Udcline_r(dcbr,iyr:t1r)=1;
                                Udcline_r(dcbr,t2r)=0;
                                in =floor(t2r/n_hours);
                                iyr=t2r;
                            else
                                Udcline_r(dcbr,iyr:n_hours*n_series)=1;
                                in=n_series;
                            end
                        else
                            if t2r< n_hours*n_series
                                Udcline_r(dcbr,iyr:t1r)=0;
                                Udcline_r(dcbr,t2r)=1;
                                in =floor(t2r/n_hours);
                                iyr=t2r;
                            else
                                Udcline_r(dcbr,iyr:n_hours*n_series)=0;
                                in=n_series;
                            end
                        end
                    case 1
                        if trans_tipo(ipol)==0

                            if t2r< n_hours*n_series
                                Udcline_p(dcbr,iyr:t1r)=1;
                                Udcline_p(dcbr,t2r)=0;
                                in =floor(t2r/n_hours);
                                iyr=t2r;
                            else
                                Udcline_p(dcbr,iyr:n_hours*n_series)=1;
                                in=n_series;
                            end
                        else
                            if t2r< n_hours*n_series
                                Udcline_p(dcbr,iyr:t1r)=0;
                                Udcline_p(dcbr,t2r)=1;
                                in =floor(t2r/n_hours);
                                iyr=t2r;
                            else
                                Udcline_p(dcbr,iyr:n_hours*n_series)=0;
                                in=n_series;
                            end
                        end
                    case 2
                        if trans_tipo(ipol)==0

                            if t2r< n_hours*n_series
                                Udcline_n(dcbr,iyr:t1r)=1;
                                Udcline_n(dcbr,t2r)=0;
                                in =floor(t2r/n_hours);
                                iyr=t2r;
                            else
                                Udcline_n(dcbr,iyr:n_hours*n_series)=1;
                                in=n_series;
                            end
                        else
                            if t2r< n_hours*n_series
                                Udcline_n(dcbr,iyr:t1r)=0;
                                Udcline_n(dcbr,t2r)=1;
                                in =floor(t2r/n_hours);
                                iyr=t2r;
                            else
                                Udcline_n(dcbr,iyr:n_hours*n_series)=0;
                                in=n_series;
                            end
                        end
                end
            end
            clear transition_rate Status_ini Vs transitio

        else
            r=rand;
            for ipol = 1:length(quali_poli)
                if quali_poli(ipol)==0
                    Status_ini(ipol) = Udcline_r(dcbr,iyr) ;
                    if Status_ini(ipol)
                        transition_rate(ipol) = 1/MTTFbrs_dc(idxs(ipol));
                        trans_tipo(ipol)=0;
                    else
                        transition_rate(ipol) = 1/MTTRbrs_dc(idxs(ipol));
                        trans_tipo(ipol)=1;
                    end

                elseif quali_poli(ipol)==1
                    Status_ini(ipol) = Udcline_p(dcbr,iyr) ;
                    if Status_ini(ipol)
                        transition_rate(ipol) = 1/MTTFbrs_dc(idxs(ipol));
                        trans_tipo(ipol)=0;
                    else
                        transition_rate(ipol) = 1/MTTRbrs_dc(idxs(ipol));
                        trans_tipo(ipol)=1;
                    end
                else
                    Status_ini(ipol) = Udcline_n(dcbr,iyr) ;
                    if Status_ini(ipol)
                        transition_rate(ipol) = 1/MTTFbrs_dc(idxs(ipol));
                        trans_tipo(ipol)=0;
                    else
                        transition_rate(ipol) = 1/MTTRbrs_dc(idxs(ipol));
                        trans_tipo(ipol)=1;
                    end
                end
            end
            dt=ceil(-1/sum(transition_rate)*log(1-r));
            rd = rand;

            transition_rate = transition_rate./sum(transition_rate);

            idx = find(cumsum(transition_rate) > rd);

            t1r = iyr+dt-1;
            t2r = t1r+1;

            ipol = idx(1);
            switch quali_poli(ipol)
                case 0
                    if trans_tipo(ipol)==0

                        if t2r< n_hours*n_series
                            Udcline_r(dcbr,iyr:t1r)=1;
                            Udcline_r(dcbr,t2r)=0;
                            in =floor(t2r/n_hours);
                            iyr=t2r;
                        else
                            Udcline_r(dcbr,iyr:n_hours*n_series)=1;
                            in=n_series;
                        end
                    else
                        if t2r< n_hours*n_series
                            Udcline_r(dcbr,iyr:t1r)=0;
                            Udcline_r(dcbr,t2r)=1;
                            in =floor(t2r/n_hours);
                            iyr=t2r;
                        else
                            Udcline_r(dcbr,iyr:n_hours*n_series)=0;
                            in=n_series;
                        end
                    end
                case 1
                    if trans_tipo(ipol)==0
                        if t2r< n_hours*n_series
                            Udcline_p(dcbr,iyr:t1r)=1;
                            Udcline_p(dcbr,t2r)=0;
                            in =floor(t2r/n_hours);
                            iyr=t2r;
                        else
                            Udcline_p(dcbr,iyr:n_hours*n_series)=1;
                            in=n_series;
                        end
                    else
                        if t2r< n_hours*n_series
                            Udcline_p(dcbr,iyr:t1r)=0;
                            Udcline_p(dcbr,t2r)=1;
                            in =floor(t2r/n_hours);
                            iyr=t2r;
                        else
                            Udcline_p(dcbr,iyr:n_hours*n_series)=0;
                            in=n_series;
                        end
                    end
                case 2
                    if trans_tipo(ipol)==0
                        if t2r< n_hours*n_series
                            Udcline_n(dcbr,iyr:t1r)=1;
                            Udcline_n(dcbr,t2r)=0;
                            in =floor(t2r/n_hours);
                            iyr=t2r;
                        else
                            Udcline_n(dcbr,iyr:n_hours*n_series)=1;
                            in=n_series;
                        end
                    else
                        if t2r< n_hours*n_series
                            Udcline_n(dcbr,iyr:t1r)=0;
                            Udcline_n(dcbr,t2r)=1;
                            in =floor(t2r/n_hours);
                            iyr=t2r;
                        else
                            Udcline_n(dcbr,iyr:n_hours*n_series)=0;
                            in=n_series;
                        end
                    end
            end
            clear transition_rate Status_ini
        end
    end
end

for i = 1:n_series
    Udclines_n{i} = Udcline_n(:,1+(i-1)*n_hours:n_hours*i);
    Udclines_r{i} = Udcline_r(:,1+(i-1)*n_hours:n_hours*i);
    Udclines_p{i} = Udcline_p(:,1+(i-1)*n_hours:n_hours*i);
end
% AC Gens
for ibr = 1:size(mpc.gen,1)
    iyr=1;dt=1;in=1;
    while in < n_series
        r=rand;
        if Uacgen(ibr,iyr)==1
            dt=ceil(-MTTFsgen_ac(ibr)*log(r));
            t1r = iyr+dt-1;
            t2r = iyr+dt;
            if t2r< n_hours*n_series
                Uacgen(ibr,iyr:t1r)=1;
                Uacgen(ibr,t2r)=0;
                in =floor(t2r/n_hours);
                iyr=t2r;
            else
                Uacgen(ibr,iyr:n_hours*n_series)=1;
                in=n_series;
            end
        else
            dt = ceil(-MTTRsgen_ac(ibr)*log(r));
            t1r = iyr+dt-1;
            t2r = iyr+dt;
            if t2r< n_hours*n_series
                Uacgen(ibr,iyr:t1r)=0;
                Uacgen(ibr,t2r)=1;
                in = floor(t2r/n_hours);
                iyr = t2r;
            else
                Uacgen(ibr,iyr:n_hours*n_series)=0;
                in = n_series;
            end
        end
    end
end
for i = 1:n_series
    Uacgens{i} = Uacgen(:,1+(i-1)*n_hours:n_hours*i);
end
% DC converters
if DCdependent_adequacy == 0
    for ibr = 1:size(Udcconv_p,1)
        iyr=1;dt=1;in =1;
        while in < n_series
            r = rand;
            if Udcconv_p(ibr,iyr)==1
                dt=ceil(-MTTFsconv_dc(ibr)*log(r));
                t1r = iyr+dt-1;
                t2r = iyr+dt;

                if t2r < n_hours*n_series
                    Udcconv_p(ibr,iyr:t1r)=1;
                    Udcconv_p(ibr,t2r)=0;
                    in = floor(t2r/n_hours);
                    iyr = t2r;
                else
                    Udcconv_p(ibr,iyr:n_hours*n_series)=1;
                    in = n_series;
                end
            else
                dt = ceil(-MTTRsconv_dc(ibr)*log(r));
                t1r = iyr+dt-1;
                t2r = iyr+dt;
                if t2r < n_hours*n_series
                    Udcconv_p(ibr,iyr:t1r)=0;
                    Udcconv_p(ibr,t2r)=1;
                    in = floor(t2r/n_hours);
                    iyr = t2r;
                else
                    Udcconv_p(ibr,iyr:n_hours*n_series) = 0;
                    in = n_series;
                end
            end
        end
    end

    for ibr = 1:size(Udcconv_n,1)
        iyr=1;dt=1;in =1;
        while in < n_series
            r = rand;
            if Udcconv_n(ibr,iyr) == 1
                dt = ceil(-MTTFsconv_dc(ibr)*log(r));
                t1r = iyr+dt-1;
                t2r = iyr+dt;

                if t2r< n_hours*n_series
                    Udcconv_n(ibr,iyr:t1r)=1;
                    Udcconv_n(ibr,t2r)=0;
                    in = floor(t2r/n_hours);
                    iyr = t2r;
                else
                    Udcconv_n(ibr,iyr:n_hours*n_series)=1;
                    in=n_series;
                end

            else
                dt = ceil(-MTTRsconv_dc(ibr)*log(r));
                t1r = iyr+dt-1;
                t2r = iyr+dt;
                if t2r < n_hours*n_series
                    Udcconv_n(ibr,iyr:t1r)=0;
                    Udcconv_n(ibr,t2r)=1;
                    in = floor(t2r/n_hours);
                    iyr = t2r;
                else
                    Udcconv_n(ibr,iyr:n_hours*n_series) = 0;
                    in = n_series;
                end
            end
        end
    end
else
    for dcbr = 1:size(mpc.convdc,1)
        iyr=1;dt=1;in=1;
        while in < n_series
            quali_conv = [1 2 ];
            r = rand;
            for ipol = 1:length(quali_conv)
                if quali_poli(ipol)==1
                    StatusC_ini(ipol) = Udcconv_p(dcbr,iyr) ;
                    if StatusC_ini(ipol)
                        transition_rate(ipol) = 1/MTTFsconv_dc(dcbr);
                        trans_tipo(ipol)=0;
                    else
                        transition_rate(ipol) = 1/MTTRsconv_dc(dcbr);
                        trans_tipo(ipol)=1;
                    end
                else
                    StatusC_ini(ipol) = Udcconv_n(dcbr,iyr) ;
                    if StatusC_ini(ipol)
                        transition_rate(ipol) = 1/MTTFsconv_dc(dcbr);
                        trans_tipo(ipol)=0;
                    else
                        transition_rate(ipol) = 1/MTTRsconv_dc((dcbr));
                        trans_tipo(ipol)=1;
                    end
                end
            end
            transitio = [];Vs = [];

            for k = 0:length(quali_conv)
                v = nchoosek([1:length(quali_conv)], k);

                for iv = 1:size(v,1)
                    Vs = [Vs; Status_ini];
                    Vs(end,v(iv,:))=not((StatusC_ini(v(iv,:))));

                    a = [transition_rate];
                    b = [ones(1,length(transition_rate))];
                    AB = [a;b];
                    dummy = 0;
                    quali_compl = (Vs(end,:)-StatusC_ini==0)*2+(abs(Vs(end,:)-StatusC_ini)>0)*1;

                    CORRE = eye(length(quali_conv) );
                    CORRE = eye(length(quali_poli) );
                    idxf = find(StatusC_ini==1);
                    for i1 = 1:length(idxf)-1
                        for i2 = i1+1:length(idxf)
                            CORRE(idxf(i1),idxf(i2)) = corr_dcpoles_adequacy;
                            CORRE(idxf(i2),idxf(i1)) = corr_dcpoles_adequacy;
                        end
                    end
                    for in1 = 1:quali_compl(1)
                        for in2 = 1:quali_compl(2)
                            Vx = [in1 in2 ];
                            Vss = (abs(Vs(end,:)-StatusC_ini)==0)*2+(abs(Vs(end,:)-StatusC_ini)>0)*1;
                            dummy = dummy + ((-1).^(sum(abs(Vss-Vx)))).*copulacdf('gaussian',[AB(in1,1) AB(in2,2)],CORRE);
                        end
                    end
                    transitio = [transitio; dummy];
                end
            end
            iqualeini = find(ismember(Vs,StatusC_ini,'rows'));
            transitio(iqualeini) = [];
            Vs(iqualeini,:)=[];

            dt=ceil(-1/sum(transitio)*log(1-r));
            rd = rand;

            transition_ = transitio./sum(transitio);

            idx = find(cumsum(transition_) > rd);

            qualeTran = idx(1);

            dquali = Vs(qualeTran,:) - StatusC_ini;
            t1r = iyr+dt-1;
            t2r = t1r+1;

            dummytrans0 = find(not(dquali == 0));

            for idum = 1:length(dummytrans0)
                ipol = dummytrans0(idum);
                switch quali_poli(ipol)
                    case 1
                        if trans_tipo(ipol)==0
                            if t2r< n_hours*n_series
                                Udcconv_p(dcbr,iyr:t1r)=1;
                                Udcconv_p(dcbr,t2r)=0;
                                in =floor(t2r/n_hours);
                                iyr=t2r;
                            else
                                Udcconv_p(dcbr,iyr:n_hours*n_series)=1;
                                in=n_series;
                            end
                        else
                            if t2r< n_hours*n_series
                                Udcconv_p(dcbr,iyr:t1r)=0;
                                Udcconv_p(dcbr,t2r)=1;
                                in =floor(t2r/n_hours);
                                iyr=t2r;
                            else
                                Udcconv_p(dcbr,iyr:n_hours*n_series)=0;
                                in=n_series;
                            end
                        end
                    case 2
                        if trans_tipo(ipol)==0
                            if t2r< n_hours*n_series
                                Udcconv_n(dcbr,iyr:t1r)=1;
                                Udcconv_n(dcbr,t2r)=0;
                                in =floor(t2r/n_hours);
                                iyr=t2r;
                            else
                                Udcconv_n(dcbr,iyr:n_hours*n_series)=1;
                                in=n_series;
                            end
                        else
                            if t2r< n_hours*n_series
                                Udcconv_n(dcbr,iyr:t1r)=0;
                                Udcconv_n(dcbr,t2r)=1;
                                in =floor(t2r/n_hours);
                                iyr=t2r;
                            else
                                Udcconv_n(dcbr,iyr:n_hours*n_series)=0;
                                in=n_series;
                            end
                        end
                end
            end
            clear transition_rate StatusC_ini Vs transitio
        end
    end
end
for i = 1:n_series
    Udcconvs_n{i} = Udcconv_n(:,1+(i-1)*n_hours:n_hours*i);
    Udcconvs_p{i} = Udcconv_p(:,1+(i-1)*n_hours:n_hours*i);
end


%%% PRINT CSV FILES FOR CSV TREE STRUCTURE

for in = 1:n_series
    % folderctgout = [root_path '\' folderout filesep 'YR_' num2str(in)];
    folderctgout = [folderout filesep 'YR_' num2str(in)];

    mkdir(folderctgout);
    mkdir([folderctgout filesep 'branch'])
    csvwrite([folderctgout filesep 'branch' filesep 'br_status.csv'],[[1:size(mpc.branch,1)];Uaclines{in}']);
    mkdir([folderctgout filesep 'gen'])
    csvwrite([folderctgout filesep 'gen' filesep 'gen_status.csv'],[[1:size(mpc.gen,1)];Uacgens{in}']);
    mkdir([folderctgout filesep 'branchdc'])
    Udcbr_r{in}=Udclines_r{in}(idexdc(find(poli==0)),1:n_hours);
    csvwrite([folderctgout filesep 'branchdc' filesep 'status_r.csv'],[[idexdc(find(poli==0))];Udcbr_r{in}'])
    Udcbr_p{in}=Udclines_p{in}(idexdc(find(poli==1)),1:n_hours);
    csvwrite([folderctgout filesep 'branchdc' filesep 'status_p.csv'],[[idexdc(find(poli==1))];Udcbr_p{in}'])
    Udcbr_n{in}=Udclines_n{in}(idexdc(find(poli==2)),1:n_hours);
    csvwrite([folderctgout filesep 'branchdc' filesep 'status_n.csv'],[[idexdc(find(poli==2))];Udcbr_n{in}'])
    mkdir([folderctgout filesep 'convdc'])
    csvwrite([folderctgout filesep 'convdc' filesep 'status_p.csv'],[[1:size(mpc.convdc,1)];Udcconvs_p{in}']);
    csvwrite([folderctgout filesep 'convdc' filesep 'status_n.csv'],[[1:size(mpc.convdc,1)];Udcconvs_n{in}']);
end

warning on
disp('Contingencies time series generation completed')
