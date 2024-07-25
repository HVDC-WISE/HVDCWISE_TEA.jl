% script for the generation of unavailabilities seriesfor N-1 and N-k
% outages including HVDC assets
clear all
clc
close all
format long

% fiel for inputting data useful for simulations
Inputs_forSeriesGeneration
%
switch type_of_unavailabilities
    case 0
        disp('unavailabilities generated only for N-1 outages ')
    otherwise
        if type_of_unavailabilities == 1
            disp('unavailabilities generated only for N-k outages ')
        else
           disp('unavailabilities generated only for both N-1 and N-k outages ')
        end
%% START generating GEV using the Block maxima technique
for ilo = 1:length(LONGI_)
    for ilat = 1: length(LATI_)
        y0 = squeeze(Wind(LONGI_(ilo),LATI_(ilat),:));

        % convert to 3s gust
        y = y0*(1.277+0.296*tanh(0.9*log10(45/3)));

        if any(isnan(y)) == 0
            if metho == 1
                n = length(y);
                % method POT (not used in the sequel)
                Thr = [0.1:0.1:max(y)*0.8];

                for jt = 1:length(Thr)
                    try
                    [parametri]=gpfit(max(eps,y-Thr(jt)));
                    ParamPOT(jt,:)=parametri;
                    fu(jt) = length(find(y > Thr(jt)))/n;
                    catch err
                        break
                    end
                end

                ParamsPOTgeo{LONGI_(ilo),LATI_(ilat)} = ParamPOT;
                ThrsGeo{LONGI_(ilo),LATI_(ilat)} = Thr(1:jt-1);
                FuPOTgeo{LONGI_(ilo),LATI_(ilat)} = fu;

                clear  ParamPOT fu
            else
                n = length(y);
                M = [7 30 365]; % weekly, monthly, yearly
                for j = 1:length(M)
                    Qj = reshape(y(1:floor(n/M(j))*M(j)),floor(n/M(j)),M(j));
                    QMj = max(Qj,[],2);

                    [parametri]=mle(QMj,'Distribution','Generalized Extreme Value');
                    ParamBM(j,:)=parametri;

                    clear Qj QMj
                end
                ParamsBMgeo{LONGI_(ilo),LATI_(ilat)} = ParamBM;
                M_Geo{LONGI_(ilo),LATI_(ilat)} = M;
            end
        end
    end
end
%%%%%%%% END generating GEV using the Block maxima technique

for i = 1:size(mpc0.branch,1)
    [intercepto]=bresenham_mi(LATIc_,LONGIc_,[mpc0.bus_coord2(mpc0.branch(i,1),1),mpc0.bus_coord2(mpc0.branch(i,2),1)],[mpc0.bus_coord2(mpc0.branch(i,1),2),mpc0.bus_coord2(mpc0.branch(i,2),2)]);
    C{i} = [LATI_(intercepto(:,1))' LONGI_(intercepto(:,2))'];
    if mpc0.branch(i,9)==0
        V{i} = ones(size(intercepto,1),1)*[19 0.2];%[;19.9 0.2]; % parametri media,sigma di lognormale
    else
        V{i} = ones(size(intercepto,1),1)*[1e5 0.2];%
    end
end
nn = size(mpc.branch,1);

ndc = 0;
for i = 1:size(mpc0.branchdc,1)
    if mpc0.branchdc(i,10)==2
        for ipol = 1:3

            ndc = ndc+1;
            [intercepto]=bresenham_mi(LATIc_,LONGIc_,[mpc0.busdc_coord2(mpc0.branchdc(i,1),1),mpc0.busdc_coord2(mpc0.branchdc(i,2),1)],[mpc0.busdc_coord2(mpc0.branchdc(i,1),2),mpc0.busdc_coord2(mpc0.branchdc(i,2),2)]);
            C{nn+ndc} = [LATI_(intercepto(:,1))' LONGI_(intercepto(:,2))'];
            if mpc0.branchdc(i,17)
                V{nn+ndc} = ones(size(intercepto,1),1)*[19 0.2];
            else
                V{nn+ndc} = ones(size(intercepto,1),1)*[1e5 0.2];%
            end
        end
    else
        quali_poli = find(mpc.branchdc(i,14:16)==1); %1 =p,2=n, 3=r
        for ipol = 1:length(quali_poli)

            ndc = ndc+1;
            [intercepto]=bresenham_mi(LATIc_,LONGIc_,[mpc0.busdc_coord2(mpc0.branchdc(i,1),1),mpc0.busdc_coord2(mpc0.branchdc(i,2),1)],[mpc0.busdc_coord2(mpc0.branchdc(i,1),2),mpc0.busdc_coord2(mpc0.branchdc(i,2),2)]);
            C{nn+ndc} = [LATI_(intercepto(:,1))' LONGI_(intercepto(:,2))'];
            if mpc0.branchdc(i,17)
                V{nn+ndc} = ones(size(intercepto,1),1)*[meanVu sigmaVu];%[;19.9 0.2]; % parametri media,sigma di lognormale
            else
                V{nn+ndc} = ones(size(intercepto,1),1)*[1e5 0.2];%
            end
        end
    end
end

Indexes_branch = [[1:size(mpc0.branch,1)]';idexdc'];

% start definition of event matrix M

M = zeros(length(tempo),length(C));

for i = 1:length(C)
   for j = 1:size(C{i})
        if any(squeeze(Wind(C{i}(j,1),C{i}(j,2),:)) > SOGLIA_M)
            eventi = find(squeeze(Wind(C{i}(j,1),C{i}(j,2),:)) > SOGLIA_M);
            M(eventi,i)=1;
        end
    end
end
% END definition of event matrix M

% start definition of correlation matrix R
CORREL = eye(length(C));
for j = 1:length(C)-1
    for i = j+1:length(C)
        nia = length(find(M(:,i)==1));
        nja = length(find(M(:,j)==1));
        n_nia = length(find(M(:,i)==0));
        n_nja = length(find(M(:,j)==0));
        ni_j = length(intersect(find(M(:,i)==1),find(M(:,j)==1)));
        n_ni_j = length(intersect(find(M(:,i)==0),find(M(:,j)==1)));
        n_i_nj = length(intersect(find(M(:,i)==1),find(M(:,j)==0)));
        n_ni_nj = length(intersect(find(M(:,i)==0),find(M(:,j)==0)));

        CORREL(i,j) = (ni_j*n_ni_nj - n_ni_j*n_i_nj)/max(eps,sqrt(nia*nja*n_nia*n_nja));
        CORREL(j,i)=CORREL(i,j);
    end
end
% start definition of correlation matrix R

W = [0.1:0.1:60]; % m/s
% START asset TR calculation
for i = 1:length(C)
    for elms=1:size(V{i},1)

        miso = V{i}(elms,1);
        viso = ( V{i}(elms,2))^2;
        mu_ = log((miso^2)/sqrt(viso+miso^2));
        sigma_ = sqrt(log(viso/(miso^2)+1));

        PV{i}(elms,:) = logncdf(W,mu_,sigma_);
        if isempty(ParamsBMgeo{C{i}(elms,2),C{i}(elms,1)})==0
            pdfx{i}(elms,:) = gevpdf(W,ParamsBMgeo{C{i}(elms,2),C{i}(elms,1)}(3,1),ParamsBMgeo{C{i}(elms,2),C{i}(elms,1)}(3,2),ParamsBMgeo{C{i}(elms,2),C{i}(elms,1)}(3,3));
        else
            pdfx{i}(elms,:) = zeros(1,length(W));
        end
        TRr{i}(elms) = 1/max(eps,trapz(W,pdfx{i}(elms,:).*PV{i}(elms,:)));
    end
    TR(i) = min(TRr{i});
end
TR(TR>1e6)=Inf;
% END asset TR calculation

%% line clustering
Thorizon=8760;
P = 1 - exp(-1./TR);
P0=P;
P = P./ Thorizon;

RHOCORR0=eye(size(M,2));

for j = 1:size(M,2)-1
    for k = j+1:size(M,2)
        if Indexes_branch(i) == Indexes_branch(k) % si tratta di due cavi dello stesso DC branch
            RHOCORR0(j,k) = corr_dcpoles_resilience;
            RHOCORR0(k,j) = corr_dcpoles_resilience;
        else
        pjk = length(find(ismember(M(:,[j k]),[1 1],'rows')))/size(M,1);
        pj = length(find(M(:,[j ])==1))/size(M,1);
        pk = length(find(M(:,[k ])==1))/size(M,1);
        RHOCORR0(j,k) = (pjk - pj*pk)/sqrt(pj*(1-pj)*pk*(1-pk));
        RHOCORR0(k,j) = (pjk - pj*pk)/sqrt(pj*(1-pj)*pk*(1-pk));
        end
    end
end
RHOCORR0(isnan(RHOCORR0))=0;

if det(RHOCORR0)< 1e-10
    RHOCORR = closest_corr(RHOCORR0);
else
    RHOCORR = RHOCORR0;
end

if length(P)> SOGLIA_LI
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% START
    nitmax = 99;

    clear GRUPOS dimensione
    RHOCORR = RHOCORR - eye(size(RHOCORR));
    D = 1 - abs(RHOCORR);
    D0 = D;ng=Inf;NG0=[1:size(D,1)];

    for i = 1:size(RHOCORR,1)
        GRUPOS{i} = i;
        dimensione(i) = 1;
        valo_corr_int(i)=1;
    end
    %%%% STEP 1: clustering based on correlation matrix R
    %
    save step1_data.mat
    D2 = D;nt = 0;
    dimensione2 = dimensione;
    valo2 = valo_corr_int;
    for par1 = min_corr_int_step1
        for par2 = max_dist_intercluster_step1
            nt = nt+1;
            coppia(nt,:)=[par1 par2];
            GRUPOS1 =GRUPOS;
            D = D2;
            dimensione = dimensione2;
            valo_corr_int = valo2;

            while (min(valo_corr_int) > par1 && min(min(D)) < par2)  || max(dimensione) < SOGLIA_LI
                ng = length(GRUPOS1);
                [x y] = find(D == min(min(D)));min(min(D));
                quali = unique(union(x,y));
                try
                    GRUPOS1{ng+1} = unique([GRUPOS1{quali}]);
                    valo_corr_int(ng+1)= sum(sum(abs(RHOCORR(GRUPOS1{ng+1},GRUPOS1{ng+1}))))/prod(length(GRUPOS1{ng+1})*length(GRUPOS1{ng+1}) - length(GRUPOS1{ng+1}));%sum(sum(abs(RHOCORR(GRUPOS{ng+1},GRUPOS{ng+1}))))/prod(length(GRUPOS{ng+1})*length(GRUPOS{ng+1}) - length(GRUPOS{ng+1}));
                    dimensione(ng+1)=length(GRUPOS1{ng+1});
                    gruppi_residui = setdiff([1:size(D,1)],quali);
                    for i = gruppi_residui
                        if i == ng+1
                            D(ng+1,i)=1;
                        else
                            D(ng+1,i) = max(max(D0(GRUPOS1{ng+1},GRUPOS1{i})));%sum(sum(D0(GRUPOS1{ng+1},GRUPOS1{i})))/prod(length(GRUPOS1{i})*length(GRUPOS1{ng+1}));%max(max(D0(GRUPOS1{ng+1},GRUPOS1{i})));%sum(sum(D0(GRUPOS{ng+1},GRUPOS{i})))/prod(length(GRUPOS{i})*length(GRUPOS{ng+1}));%%max(max(D0(GRUPOS{ng+1},GRUPOS{i})));%sum(sum(D0(GRUPOS{ng+1},GRUPOS{i})))/prod(length(GRUPOS{i})*length(GRUPOS{ng+1}));%%
                        end
                        D(i,ng+1) = D(ng+1,i);
                    end

                    D(ng+1,ng+1)=1;
                    GRUPOS1(quali)=[];
                    D(quali,:)=[];
                    D(:,quali)=[];
                    dimensione(quali) = [];
                    valo_corr_int(quali)=[];

                catch err
                    keyboard
                end
            end
            silou = define_silhouette(GRUPOS1,D0);
            Dimens{nt}=dimensione;
            Valo_correl{nt}=valo_corr_int;
            Silouets{nt}=silou;
            QuantilSil(nt,:) = quantile(silou,[0.05 0.25 0.5 0.75 0.95]);
            QuantilCorrInt(nt,:) = quantile(valo_corr_int(dimensione>1),[0.05 0.25 0.5 0.75 0.95]);
            GRUPPING(nt).groups = GRUPOS1;
        end
    end
    % selection of the pair of parameters bringing to best clustering
    % performance
    Indice = (1-pesoSil)*QuantilCorrInt + (pesoSil)*QuantilSil;
    idq = find(Indice(:,1) == max(Indice(:,1)));
    GRUPOS = GRUPPING(idq(1)).groups;
    dimensione = Dimens{idq(1)};
    valo_corr_int = Valo_correl{idq(1)};
    silouette = Silouets{idq(1)};

    try
        if figu_plot == 1
            figure
            plot(Silouets{idq(1)}),title('silouette metrics for lines'),xlabel('LINE ID')
            RHOCORR = RHOCORR +eye(size(RHOCORR));
            %
            figure,
            subplot(2,2,1)
            plot(dimensione),title('cluster cardinality')
            subplot(2,2,3)
            plot(valo_corr_int),title('intebnal correlation value')
            subplot(2,2,[2 , 4])
            bar3(RHOCORR([GRUPOS{find(dimensione>1)}],[GRUPOS{find(dimensione==1)}])),title({'residual correlation between clustered and not clustered lines ','nr of card-1 clustersi: ' num2str(length(find(dimensione==1)))})

            gruppi_gr = find(dimensione > 1);
            for i1 = 1:length(gruppi_gr)-1
                for i2 = i1+1:length(gruppi_gr)
                    RHOGRAF(i1,i2) =  mean(mean(abs(RHOCORR(GRUPOS{gruppi_gr(i1)},GRUPOS{gruppi_gr(i2)}))));
                    RHOGRAF(i2,i1) =RHOGRAF(i1,i2);
                end
            end

            figure
            b = bar3(RHOGRAF)
            for k = 1:length (b)
                zdata = b(k).ZData;
                b(k).CData = zdata;
                b(k).FaceColor = 'interp';
            end
            colorbar
        end

        catch err
    end
    % STEp 2:separation of big clusters( dim > SOGLIA_LI) based on topological
    % information

    stazioni = mpc.bus_name;

    linee_ = [1:size(mpc.branch,1)];

    for il = 1:length(linee_)
        linee{il} = ['ACBR_' num2str(linee_(il))];
        estremi10{il} = stazioni{mpc.branch(il,1)};
        estremi20{il} = stazioni{mpc.branch(il,2)};
        nomi{1,il} = linee{il} ;
        nomi{2,il} = ['Line_' linee{il}] ;
        tiplinee{il} = 1;
    end

    nn = length(linee_);

    stazionidc = mpc.busdc_name;

    lineedc_ = idexdc;

    for il = 1:length(lineedc_)
        if poli(il) == 0
            linee{il+nn} = ['DCBR_' num2str(idexdc(il)) '_R'];
        elseif poli(il) == 1
            linee{il+nn} = ['DCBR_' num2str(idexdc(il)) '_P'];
        elseif   poli(il) == 2
            linee{il+nn} = ['DCBR_' num2str(idexdc(il)) '_N'];
        end
        estremi10{il+nn} = stazionidc{mpc.branchdc(idexdc(il),1)};
        estremi20{il+nn} = stazionidc{mpc.branchdc(idexdc(il),2)};
        nomi{1,il+nn} = linee{il+nn} ;
        nomi{2,il+nn} = ['Line_' linee{il+nn}] ;
        tiplinee{il+nn} = 2;
    end

    ns = length(stazioni);
    %%
    nl=0;
    MATR = [];
    for ili = 1:length(linee)
        if tiplinee{ili} ==1
            try
                if  all(isnan(linee{ili})==0) && not(isnumeric(linee{ili}))
                    idx1 = find(ismember(stazioni,estremi10{ili}));
                    idx2 = find(ismember(stazioni,estremi20{ili}));
                    if isempty(strfind(linee{ili},'VT_ERROR')) && isempty(idx1)==0 &&  isempty(idx2)==0 && not(idx1 == idx2)
                        nl = nl+1;

                        MATR = [MATR; zeros(1,length(stazioni)+length(stazionidc))];
                        MATR(end,idx1(1)) = 1;
                        MATR(end,idx2(1)) = -1;

                        nome_linee{nl}=linee{ili};
                        estremi1{nl}=estremi10{ili};
                        estremi2{nl}=estremi20{ili};
                    end
                end
            catch err
                keyboard
            end
        else
            if  all(isnan(linee{ili})==0) && not(isnumeric(linee{ili}))
                idx1 = find(ismember(stazionidc,estremi10{ili}));
                idx2 = find(ismember(stazionidc,estremi20{ili}));
            if isempty(strfind(linee{ili},'VT_ERROR')) && isempty(idx1)==0 &&  isempty(idx2)==0 && not(idx1 == idx2)
                nl = nl+1;

                MATR = [MATR; zeros(1,length(stazioni)+length(stazionidc))];
                MATR(end,ns+idx1(1)) = 1;
                MATR(end,ns+idx2(1)) = -1;

                nome_linee{nl}=linee{ili};
                estremi1{nl}=estremi10{ili};
                estremi2{nl}=estremi20{ili};
            end
        end
    end
end

seleziona=[];quali = [];quali_ma=[];
for isap = 1:size(nomi,2)
   if isempty(find(ismember(nome_linee,nomi{1,isap})))==0
        seleziona = [seleziona find(ismember(nome_linee,nomi{1,isap}))];
        quali = [quali find(sum((MATR(find(ismember(nome_linee,nomi{1,isap})),:)),1)>0)];
        CONNEX(find(ismember(nome_linee,nomi{1,isap})))=isap;
    else
        quali_ma = [quali_ma isap];
        disp([nomi{1,isap} '--' nomi{2,isap}])
    end
end
MATR_1 = zeros(size(nomi,2),length(stazioni)+length(stazionidc));
elimina=[];selecio=[];
for inol = 1:size(nomi,2)
    if isempty(find(CONNEX == inol))==0
        try
            if length(find(CONNEX == inol))==1
                MATR_1(inol,:) = (MATR(find(CONNEX == inol),:));
            elseif length(find(CONNEX == inol))==2
                dummy=find(CONNEX == inol);
                MATR_1(inol,:) = (MATR(dummy(1),:));
            else
                dummys = find(CONNEX == inol);
                qualis = find(sum(abs(MATR(dummys,:)),1)>0);
                qualis = setdiff(qualis,qualis(sum(abs(MATR(dummys,qualis)),1) == max(sum(abs(MATR(dummys,qualis)),1))))
                if length(qualis)==3
                    MATR_1(inol,qualis(1)) = 2;
                    MATR_1(inol,qualis(2)) = -1;
                    MATR_1(inol,qualis(3)) = -1;
                else
                    MATR_1(inol,qualis(1)) = 1;
                    MATR_1(inol,qualis(2)) = -1;
               end
            end
        catch err
            keyboard
        end
        selecio = [selecio inol];
    else
        elimina = [elimina inol];
    end
end

MATR_2 = MATR_1;
MATR_2(elimina,:)=[];
quali_staz_elim = find(all(MATR_2==0,1));
MATR_2(:,quali_staz_elim)=[];
stazioni_ = [stazioni stazionidc];
stazioni_(quali_staz_elim)=[];
Stopo = MATR_2;

for icu = 1:length(stazioni_)
   x = ones(length(stazioni_),1);
   x(icu)=-1;
   y = Stopo*x;
   idx_quali = find(y ~=0);
   cutsets{icu} = selecio(idx_quali);
   dimensione_clu(icu)=length(idx_quali);
end
Ncutsets = length(cutsets);

nsg=0;
disp('analysis of clusters so far generated')
for ig = 1:length(GRUPOS)
    disp(['CLUSTER ' num2str(ig)])
        MATRICE_AG=[];
        sottogruppo=[];
        dim_sottogruppo =[];
        valo_corr_sottogruppo =[];
    if length(GRUPOS{ig})> 1 && length(GRUPOS{ig}) <= SOGLIA_LI


    elseif length(GRUPOS{ig}) > SOGLIA_LI

        quali_lin = find(ismember(selecio,GRUPOS{ig}));
        %

        S_red = Stopo(quali_lin,:);
        quali_st = find(sum(abs(S_red),1) > 0);
        stazioni_2 = stazioni_(quali_st);
        S_red = S_red(:,quali_st);

        restanti = (GRUPOS{ig});
        nsg = 0;
        while (length(restanti)) > SOGLIA_LI
            [a b] = sort(sum(abs(S_red),1),'descend');
            nsg = nsg + 1;
            try
                sottogruppo{nsg} = (selecio(quali_lin(find(abs(S_red(:,b(1)))>0))));
            catch err
                keyboard
            end
            dim_sottogruppo(nsg) = length(sottogruppo{nsg});
            valo_corr_sottogruppo(nsg) = sum(sum(RHOCORR(sottogruppo{nsg},sottogruppo{nsg})-eye((dim_sottogruppo(nsg)))))/(dim_sottogruppo(nsg)^2 - dim_sottogruppo(nsg));
            restanti = setdiff(restanti,sottogruppo{nsg});
            quali = find(ismember(selecio(quali_lin),sottogruppo{nsg}));
            S_red((quali),:)=[];
            quali_lin(quali)=[];
        end

        nsg = nsg + 1;
        sottogruppo{nsg} = restanti;
        dim_sottogruppo(nsg) = length(restanti);
        valo_corr_sottogruppo(nsg) = sum(sum(RHOCORR(sottogruppo{nsg},sottogruppo{nsg})-eye((dim_sottogruppo(nsg)))))/(dim_sottogruppo(nsg)^2 - dim_sottogruppo(nsg));

        SS = eye(nsg,nsg);
        elim=[];nsg1 = nsg;fatti=[];
        for iinsg1 = 1:nsg-1
            for iinsg2 = iinsg1+1:nsg
                ncut = 0;sottogruppo0=cell(0);dim_sottogruppo0=[];valo_corr_sottogruppo0=[];
                for icut = 1:length(cutsets)
                    if (all(ismember(cutsets{icut},union(sottogruppo{iinsg1},sottogruppo{iinsg2}))) && any(ismember(cutsets{icut},sottogruppo{iinsg2})) && any(ismember(cutsets{icut},sottogruppo{iinsg1}))) && dim_sottogruppo(iinsg1)+dim_sottogruppo(iinsg2) <= SOGLIA_LI && ismember(iinsg1,fatti)  ==0 && ismember(iinsg2,fatti)  ==0
                        ncut = ncut +1;

                        sottogruppo0{ncut} = [sottogruppo{iinsg1} sottogruppo{iinsg2}];
                        AGGREGA(ncut,:) = [iinsg1 iinsg2];
                        dim_sottogruppo0(ncut) = length(sottogruppo0{ncut});
                        valo_corr_sottogruppo0(ncut) = sum(sum(RHOCORR(sottogruppo0{ncut},sottogruppo0{ncut})-eye((dim_sottogruppo0(ncut)))))/(dim_sottogruppo0(ncut)^2 - dim_sottogruppo0(ncut));
                    end
                end
                if ncut > 0
                    quale_cut = find(valo_corr_sottogruppo0==max(valo_corr_sottogruppo0));
                    nsg1 = nsg1+1;
                        sottogruppo{nsg1} = [sottogruppo0{quale_cut(1)}];
                        fatti = [fatti AGGREGA(quale_cut(1),:)];
                    dim_sottogruppo(nsg1) = length(sottogruppo{nsg1});
                    valo_corr_sottogruppo(nsg1) = sum(sum(RHOCORR(sottogruppo{nsg1},sottogruppo{nsg1})-eye((dim_sottogruppo(nsg1)))))/(dim_sottogruppo(nsg1)^2 - dim_sottogruppo(nsg1));
                end
                clear AGGREGA
            end
        end
        sottogruppo(fatti)=[];
        dim_sottogruppo(fatti) =[];
        valo_corr_sottogruppo(fatti) =[];

        disp('BIG CLUSTER detected: let''s use topological information')
   end
   AGGREGAZ{ig}=MATRICE_AG;
   DETTAGLI_GRUPPO(ig).sottogruppo=sottogruppo;
   DETTAGLI_GRUPPO(ig).dim_sottogruppo=dim_sottogruppo;
   DETTAGLI_GRUPPO(ig).valo_corr_sottogruppo=valo_corr_sottogruppo;
   DETTAGLI_GRUPPO(ig).cutset =[];
    for icut = 1:length(cutsets)
        if any(ismember(cutsets{icut},GRUPOS{ig}))
            DETTAGLI_GRUPPO(ig).cutset = [DETTAGLI_GRUPPO(ig).cutset  icut];
        end
    end
    DETTAGLI_GRUPPO(ig).disaggregato=0;
    clear S  nomis dim_sottogruppo sottogruppo  valo_corr_sottogruppo
end
NG = length(GRUPOS);elim=[];
for ig = 1:length(GRUPOS)
    DETTAGLI_GRUPPO(NG).disaggregato=0;
    if isempty(DETTAGLI_GRUPPO(ig).sottogruppo)==0
        for ik = 1:length(DETTAGLI_GRUPPO(ig).sottogruppo)
            NG = NG+1;
            GRUPOS{NG}=DETTAGLI_GRUPPO(ig).sottogruppo{ik};
            dimensione(NG)=DETTAGLI_GRUPPO(ig).dim_sottogruppo(ik);
            valo_corr_int(NG)=DETTAGLI_GRUPPO(ig).valo_corr_sottogruppo(ik);
            DETTAGLI_GRUPPO(NG)=DETTAGLI_GRUPPO(ig);
            DETTAGLI_GRUPPO(NG).disaggregato=1;
        end
        elim = [elim ig];
    end
end
GRUPOS(elim)=[];
dimensione(elim)=[];
valo_corr_int(elim)=[];
DETTAGLI_GRUPPO(elim)=[];
nit=1;
try
    if figu_plot == 1
        figure,
        subplot(2,2,1)
        plot(dimensione),title('cluster cardinality')
        subplot(2,2,3)
        plot(valo_corr_int),title('internal correlation ')
        subplot(2,2,[2,4])
        bar3(RHOCORR([GRUPOS{find(dimensione>1)}],[GRUPOS{find(dimensione==1)}]))
        title({'residual correlation between clustered lines and not-clustered lines ','number of clusters with cardinality 1: ' num2str(length(find(dimensione==1)))})
    end
catch err
end
silouette_step2 = define_silhouette(GRUPOS,D0);

save test3_cluster.mat


RHOGRU0 = ones(length(GRUPOS),length(GRUPOS));
fatti = 1;
%
% STEP 3: definition of aggregations between clusters with reduction of single clusters and residual correlation between clusters %%%
if length(find(dimensione == 1))>0
    while max(max(RHOCORR([GRUPOS{find(dimensione>1)}],[GRUPOS{find(dimensione==1)}]))) > corr_residua_interclusters_step3 && nit < nmax_it_step3 && isempty(fatti)==0
        nit = nit+1;
        max(median(RHOGRU0(intersect(find([DETTAGLI_GRUPPO.disaggregato]==0),find(dimensione>1)),intersect(find([DETTAGLI_GRUPPO.disaggregato]==0),find(dimensione>1)))));
        length(GRUPOS);
        NGS = length(GRUPOS);
        SS = eye(NGS,NGS);
        elim=[];nsg1 = NGS;MATRICE_AGG=[];fatti=[];
        for ig1 = 1:NGS-1
            %     ig1
            for ig2 = ig1+1:NGS

                correlazio = mean(mean(abs(RHOCORR(GRUPOS{ig1},GRUPOS{ig2}))));
                correlazio_max = max(max(abs(RHOCORR(GRUPOS{ig1},GRUPOS{ig2}))));

                RHOGRU(ig1,ig2) = correlazio;
                COND(ig1,ig2)=(dimensione(ig1)+dimensione(ig2) <= SOGLIA_LI);
                quali_cutsets = intersect(DETTAGLI_GRUPPO(ig1).cutset,DETTAGLI_GRUPPO(ig2).cutset);
                RHOTOPO(ig1,ig2) =0;
                for icut = 1:length(quali_cutsets)

                    RHOTOPO(ig1,ig2) = ((all(ismember(cutsets{quali_cutsets(icut)},union(GRUPOS{ig1},GRUPOS{ig2}))) && any(ismember(cutsets{quali_cutsets(icut)},GRUPOS{ig2})) && any(ismember(cutsets{quali_cutsets(icut)},GRUPOS{ig1})))  && (correlazio > avg_corr_min_interclustertopo_step3 || correlazio_max > correlazio_max_interclustertopo_step3));

                    if RHOTOPO(ig1,ig2) == 1
                        break
                    end
                end
                %
                RHOGRU(ig2,ig1) =RHOGRU(ig1,ig2);
                RHOTOPO(ig2,ig1) =RHOTOPO(ig1,ig2);
                COND(ig2,ig1)=COND(ig1,ig2);
            end
        end

        PERFORMANCE_IDX = (RHOGRU.*(1-peso_topologia.*(RHOTOPO > 0))+RHOTOPO.*(peso_topologia)).*COND;


        [x y] = find(PERFORMANCE_IDX == max(max(PERFORMANCE_IDX)));


        nsg1 = nsg1 +1;

        GRUPOS{nsg1} = [GRUPOS{x(1)} GRUPOS{y(1)}];

        dimensione(nsg1) = length(GRUPOS{nsg1});
        valo_corr_int(nsg1) = sum(sum(abs(RHOCORR(GRUPOS{nsg1},GRUPOS{nsg1}))-eye((dimensione(nsg1)))))/(dimensione(nsg1)^2 - dimensione(nsg1));
        DETTAGLI_GRUPPO(nsg1).cutset = union(DETTAGLI_GRUPPO(x(1)).cutset,DETTAGLI_GRUPPO(y(1)).cutset);
        DETTAGLI_GRUPPO(nsg1).disaggregato = DETTAGLI_GRUPPO(x(1)).disaggregato+DETTAGLI_GRUPPO(y(1)).disaggregato;
        fatti = [ x(1) y(1)];
        MATRICE_AGG = [MATRICE_AGG; x(1) y(1) nsg1];


        GRUPOS(fatti)=[];
        dimensione(fatti)=[];
        DETTAGLI_GRUPPO(fatti)=[];
        valo_corr_int(fatti)=[];
        RHOGRU0=RHOGRU;
        clear RHOGRU RHOTOPO RHOTOT COND
        if length(find(dimensione==1))==0
            break
        end
    end
end

else
     GRUPOS{1}=[1:length(P)];
end

PRAND=[];V=[];
RHOCORR = RHOCORR + eye(size(RHOCORR));
elimin = find(dimensione==0);
GRUPOS(elimin)=[];
DETTAGLI_GRUPPO(elimin)=[];
valo_corr_int(elimin)=[];
%%%%% END CLUSTERING METHOD

%%% START IDENTIFICATION OF THE CONTINGENCIES
save testctg.mat

for g = 1:length(GRUPOS)
    VA=zeros(1,length(GRUPOS{g}));
    if det(RHOCORR(GRUPOS{g},GRUPOS{g})) < 1e-2
        COP = closest_corr(RHOCORR(GRUPOS{g},GRUPOS{g}));
    else
        COP = RHOCORR(GRUPOS{g},GRUPOS{g});
    end

    for i = 1:length(GRUPOS{g})
        dummy=zeros(1,length(P));
        dummy(GRUPOS{g}(i))=1;
        %         VA = [VA; dummy];
        if P(GRUPOS{g}(i)) > epsi
            V = [V; dummy];
            PRAND = [PRAND; P(GRUPOS{g}(i))];
        end
    end

    v2 = [1:length(GRUPOS{g})]';
    E=[];EPRAND=[];

    for i = 2:min(kmax,length(GRUPOS{g}))
        SS = [];v=[];
        for jj = 1:size(v2,1)
            SS = [SS jj];
            dummys = setdiff([1:length(GRUPOS{g})],intersect(find([1:length(GRUPOS{g})] <= max(v2(jj,:))),unique(reshape(v2(SS,:),size(v2(SS,:),2)*size(v2(SS,:),1),1))));
            if length(dummys)>0
                v = [v;repmat(v2(jj,:),size(nchoosek(dummys,1),1),1) nchoosek(dummys,1)];
            end
        end
        pran=[];vectore=ones(size(v,1),length(GRUPOS{g}));
        for vv = 1:size(v,1)
            dummy=zeros(1,length(P));
            dummy(GRUPOS{g}(v(vv,:)))=1;
            V = [V; dummy];
            vectore(vv,v(vv,:))=P(GRUPOS{g}(v(vv,:)));
        end
        vectore(vectore==0)=eps;

            [dummys dum] = my_copulacdf_fast('Gaussian',vectore,COP);

        PRAND=[PRAND; dummys];
        pran=[pran; dummys];

        filtered = find(pran >= epsi);
        filtered2n = find(PRAND < epsi);

        v2 = v(filtered,:);
        E = [E; V(filtered2n,:)];
        EPRAND = [EPRAND; PRAND(filtered2n,:)];
        [EE iE ]= unique(E,'rows');
        EPRAND = EPRAND(iE);
        E = EE;

    end
    idxF = find(PRAND >= epsi);
    V = V(idxF,:);
    PRAND=PRAND(idxF);

    for jj = 1:length(GRUPOS{g})
        idx = nchoosek([1:length(GRUPOS{g})],jj);
        for jjj = 1:size(idx,1)
            dummy=zeros(1,length(GRUPOS{g}));
            dummy(idx(jjj,:))=1;
            VA = [VA;dummy];

        end
    end

    VAA{g}=VA;
    VAP = 1-VA;
    DV = ones(1,length(GRUPOS{g})) - VAP;

    idx_ = find(all(DV>=0,2));
    S_qua = length(find(ones(1,length(GRUPOS{g}))>0)) - sum( ones(1,length(GRUPOS{g})) - DV(idx_,:),2);
    for jj = 1:length(idx_)
        vectore(jj,find(VAA{g}(idx_(jj),:)==1))=P(GRUPOS{g}(find(VAA{g}(idx_(jj),:)==1)));
        vectore(jj,find(VAA{g}(idx_(jj),:)==0))=1;
    end
    [vectore ia ic] = unique(vectore,'rows');
    vectore(vectore==0)=eps;

    if size(vectore,2)>1
        S_qua = S_qua(ia);
        [Ps errS] = my_copulacdf_fast('Gaussian',vectore,COP);
        NULLCTG(g) = max(0,sum(((-1).^S_qua).*Ps));
         Parallelepipedo{g,1} = Ps;
    Parallelepipedo{g,2} = vectore;
        Parallelepipedo{g,3} = errS;
    else
        NULLCTG(g) = (1-P(GRUPOS{g}));
        Ps = [];
    end
    clear VA
end

VB1=[];
J = 3-2*max(min(1,(max(TR)-500)/(1500 - 500)),0);
%%% START IDENTIFICATION OF THE CONTINGENCIES
% START - CALCULATION OF CTG PROBABILITIES
parfor h = 1:size(V,1)
    VB1(h,:) =V(h,:);
    dimensione_ctg(h) = sum(V(h,:));

        [ A B C D E F G] = multiple_prob_2( GRUPOS,V,RHOCORR,P,h,VAA,NULLCTG,J,Parallelepipedo);

    PROB1s(h,:) = [A B C max(0,min([A B C])) min(1,max([A B C])) 0 D E max(0,min([A D E])) max(0,min([A D E])) 0 F G max(0,min([A F G])) min(1,max([A F G]))];

end

PROB1 = PROB1s(:,1);

PROB1yr = PROB1*8760;
% END - CALCULATION OF CTG PROBABILITIES
if figu_plot == 1
    figure
    bar(PROB1)
    title('probabilities of analysed combinations')
    xlabel('CTG ID'),ylabel('probability')

    for i = 1:size(V,2)
        idxi1 = find(VB1(:,i)==1);
        somma1(i) = sum(PROB1(idxi1));
        TEST{i} = PROB1(idxi1);
        INDEX{i} = (idxi1);
    end
    figure
    bar([P0'./8760 somma1'])
    title('verification of toal probabilities')
    xlabel('line ID')
    ylabel('p1-hour probability of contingencies')
    legend({' P original' 'P from copula application'})
    grid on
end

nm_k = intersect(find(sum(VB1,2)>1),find(PROB1yr>0));
nm_1 = intersect(find(sum(VB1,2)==1),find(PROB1yr>0));
PROBnm_k=PROB1(nm_k);
PROBnm_kyr=PROB1yr(nm_k);
PROBnm_1=PROB1(nm_1);
VBnm_1 = VB1(nm_1,:);
VBnm_k = VB1(nm_k,:);

%%%
%folderout = ['outputsNmk' strrep(char(datetime),':','p') ]; % FIXME commented
folderout = output_folder; % FIXME added
mkdir(folderout);
% print file csv for ctg pobabilities
csvwrite([folderout filesep 'ctg_probs.csv'],[[1:length(nm_k)]' PROBnm_kyr]);
%
parfor h = 1:length(PROBnm_k)
    folderctgout = [folderout filesep 'CTG_' num2str(h)];
    mkdir(folderctgout);
    U = ones(size(mpc.branch,1)+length(idexdc),8760);
    fh = find(VBnm_k(h,:)>0);
    for hh = 1:length(fh)
        if fh(hh) > size(mpc.branch,1)
            U(fh(hh),1:MTTRbrs_dc_resil(fh(hh)-size(mpc.branch,1)))=0;
        else
            U(fh(hh),1:MTTRbrs_ac_resil(fh(hh)))=0;
        end
    end

    %%%PRINT FILE FOR CSV TREE STRUCTURE
    nr_aclines = size(mpc.branch,1);
    mkdir([folderctgout filesep 'branch'])
    csvwrite([folderctgout filesep 'branch' filesep 'br_status.csv'],[[1:nr_aclines];U(1:nr_aclines,:)']);
    mkdir([folderctgout filesep 'branchdc'])
    nr_dcbr_r = nr_aclines+find(poli == 0);
    csvwrite([folderctgout filesep 'branchdc' filesep 'status_r.csv'],[idexdc(find(poli == 0)); U(nr_dcbr_r,:)']);
    nr_dcbr_p = nr_aclines+find(poli == 1);
    csvwrite([folderctgout filesep 'branchdc' filesep 'status_p.csv'],[idexdc(find(poli == 1)); U(nr_dcbr_p,:)']);
    nr_dcbr_n = nr_aclines+find(poli == 2);
    csvwrite([folderctgout filesep 'branchdc' filesep 'status_n.csv'],[idexdc(find(poli == 2)); U(nr_dcbr_n,:)']);

end
end

%%% START - GENERATION of N-1 and HVDC unavailabilities
if type_of_unavailabilities == 0 % FIXME 2 instead of 0
    %folderout = ['outputsNm1_' strrep(char(datetime) ,':','p')]; % FIXME
    folderout = output_folder; % FIXME
    mkdir(folderout);

    A=[];BV=[];
    clear Uacline Uacgen  Udcconv_p Udcconv_n

    Udcline_r = ones(size(mpc.branchdc,1),8760*N);
    Udcline_p= ones(size(mpc.branchdc,1),8760*N);
    Udcline_n = ones(size(mpc.branchdc,1),8760*N);
    Uacline= ones(size(mpc.branch,1),8760*N);
    Uacgen= ones(size(mpc.gen,1),8760*N);
    Udcconv_p= ones(size(mpc.convdc,1),8760*N);
    Udcconv_n= ones(size(mpc.convdc,1),8760*N);

    % AC  branches
    for ibr = 1:size(mpc.branch,1)
        iyr=1;dt=1;in=1;
        while in < N
            r=rand;
            if Uacline(ibr,iyr)==1
                dt=ceil(-MTTFbrs_ac(ibr)*log(r));
                t1r = iyr+dt-1;
                t2r = iyr+dt;
                % iyr = t2r;
                if t2r< 8760*N
                    Uacline(ibr,iyr:t1r)=1;
                    Uacline(ibr,t2r)=0;
                    in = floor(t2r/8760);
                    iyr=t2r;
                else
                    Uacline(ibr,iyr:8760*N)=1;
                    in=N;
                end
            else
                dt=ceil(-MTTRbrs_ac(ibr)*log(r));
                t1r = iyr+dt-1;
                t2r = iyr+dt;
                % iyr = t2r;
                if t2r< 8760*N
                    Uacline(ibr,iyr:t1r)=0;
                    Uacline(ibr,t2r)=1;
                    in =floor(t2r/8760);
                    iyr=t2r;
                else
                    Uacline(ibr,iyr:8760*N)=0;
                    in=N;
                end
            end
        end
    end
    for i = 1:N
        Uaclines{i} = Uacline(:,1+(i-1)*8760:8760*i);
    end
    % DC branches
    for dcbr = 1:size(mpc.branchdc,1)
        quali_poli = poli(find(idexdc==dcbr));
        idxs = find(idexdc==dcbr);
        iyr=1;dt=1;in=1;
        while in <N
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

                                    if t2r< 8760*N
                                        Udcline_r(dcbr,iyr:t1r)=1;
                                        Udcline_r(dcbr,t2r)=0;
                                        in =floor(t2r/8760);
                                        iyr=t2r;
                                    else
                                        Udcline_r(dcbr,iyr:8760*N)=1;
                                        in=N;
                                    end
                                else
                                    if t2r< 8760*N
                                        Udcline_r(dcbr,iyr:t1r)=0;
                                        Udcline_r(dcbr,t2r)=1;
                                        in =floor(t2r/8760);
                                        iyr=t2r;
                                    else
                                        Udcline_r(dcbr,iyr:8760*N)=0;
                                        in=N;
                                    end
                                end
                            case 1
                                if trans_tipo(ipol)==0

                                    if t2r< 8760*N
                                        Udcline_p(dcbr,iyr:t1r)=1;
                                        Udcline_p(dcbr,t2r)=0;
                                        in =floor(t2r/8760);
                                        iyr=t2r;
                                    else
                                        Udcline_p(dcbr,iyr:8760*N)=1;
                                        in=N;
                                    end
                                else
                                    if t2r< 8760*N
                                        Udcline_p(dcbr,iyr:t1r)=0;
                                        Udcline_p(dcbr,t2r)=1;
                                        in =floor(t2r/8760);
                                        iyr=t2r;
                                    else
                                        Udcline_p(dcbr,iyr:8760*N)=0;
                                        in=N;
                                    end
                                end
                            case 2
                                if trans_tipo(ipol)==0

                                    if t2r< 8760*N
                                        Udcline_n(dcbr,iyr:t1r)=1;
                                        Udcline_n(dcbr,t2r)=0;
                                        in =floor(t2r/8760);
                                        iyr=t2r;
                                    else
                                        Udcline_n(dcbr,iyr:8760*N)=1;
                                        in=N;
                                    end
                                else
                                    if t2r< 8760*N
                                        Udcline_n(dcbr,iyr:t1r)=0;
                                        Udcline_n(dcbr,t2r)=1;
                                        in =floor(t2r/8760);
                                        iyr=t2r;
                                    else
                                        Udcline_n(dcbr,iyr:8760*N)=0;
                                        in=N;
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

                                        if t2r< 8760*N
                                            Udcline_r(dcbr,iyr:t1r)=1;
                                            Udcline_r(dcbr,t2r)=0;
                                            in =floor(t2r/8760);
                                            iyr=t2r;
                                        else
                                            Udcline_r(dcbr,iyr:8760*N)=1;
                                            in=N;
                                        end
                                    else
                                        if t2r< 8760*N
                                            Udcline_r(dcbr,iyr:t1r)=0;
                                            Udcline_r(dcbr,t2r)=1;
                                            in =floor(t2r/8760);
                                            iyr=t2r;
                                        else
                                            Udcline_r(dcbr,iyr:8760*N)=0;
                                            in=N;
                                        end
                                    end
                                case 1
                                    if trans_tipo(ipol)==0
                                        if t2r< 8760*N
                                            Udcline_p(dcbr,iyr:t1r)=1;
                                            Udcline_p(dcbr,t2r)=0;
                                            in =floor(t2r/8760);
                                            iyr=t2r;
                                        else
                                            Udcline_p(dcbr,iyr:8760*N)=1;
                                            in=N;
                                        end
                                    else
                                        if t2r< 8760*N
                                            Udcline_p(dcbr,iyr:t1r)=0;
                                            Udcline_p(dcbr,t2r)=1;
                                            in =floor(t2r/8760);
                                            iyr=t2r;
                                        else
                                            Udcline_p(dcbr,iyr:8760*N)=0;
                                            in=N;
                                        end
                                    end
                                case 2
                                    if trans_tipo(ipol)==0
                                        if t2r< 8760*N
                                            Udcline_n(dcbr,iyr:t1r)=1;
                                            Udcline_n(dcbr,t2r)=0;
                                            in =floor(t2r/8760);
                                            iyr=t2r;
                                        else
                                            Udcline_n(dcbr,iyr:8760*N)=1;
                                            in=N;
                                        end
                                    else
                                        if t2r< 8760*N
                                            Udcline_n(dcbr,iyr:t1r)=0;
                                            Udcline_n(dcbr,t2r)=1;
                                            in =floor(t2r/8760);
                                            iyr=t2r;
                                        else
                                            Udcline_n(dcbr,iyr:8760*N)=0;
                                            in=N;
                                        end
                                    end
                                end
                                clear transition_rate Status_ini
                            end
                        end
                    end

                    for i = 1:N
                        Udclines_n{i} = Udcline_n(:,1+(i-1)*8760:8760*i);
                        Udclines_r{i} = Udcline_r(:,1+(i-1)*8760:8760*i);
                        Udclines_p{i} = Udcline_p(:,1+(i-1)*8760:8760*i);
                    end
                    % AC Gens
                    for ibr = 1:size(mpc.gen,1)
                        iyr=1;dt=1;in=1;
                        while in < N
                            r=rand;
                            if Uacgen(ibr,iyr)==1
                                dt=ceil(-MTTFsgen_ac(ibr)*log(r));
                                t1r = iyr+dt-1;
                                t2r = iyr+dt;
                                if t2r< 8760*N
                                    Uacgen(ibr,iyr:t1r)=1;
                                    Uacgen(ibr,t2r)=0;
                                    in =floor(t2r/8760);
                                    iyr=t2r;
                                else
                                    Uacgen(ibr,iyr:8760*N)=1;
                                    in=N;
                                end
                            else
                                dt = ceil(-MTTRsgen_ac(ibr)*log(r));
                                t1r = iyr+dt-1;
                                t2r = iyr+dt;
                                if t2r< 8760*N
                                    Uacgen(ibr,iyr:t1r)=0;
                                    Uacgen(ibr,t2r)=1;
                                    in = floor(t2r/8760);
                                    iyr = t2r;
                                else
                                    Uacgen(ibr,iyr:8760*N)=0;
                                    in = N;
                                end
                            end
                        end
                    end
                    for i = 1:N
                        Uacgens{i} = Uacgen(:,1+(i-1)*8760:8760*i);
                    end
                    % DC converters
                    if DCdependent_adequacy == 0
                        for ibr = 1:size(Udcconv_p,1)
                            iyr=1;dt=1;in =1;
                            while in < N
                                r = rand;
                                if Udcconv_p(ibr,iyr)==1
                                    dt=ceil(-MTTFsconv_dc(ibr)*log(r));
                                    t1r = iyr+dt-1;
                                    t2r = iyr+dt;

                                    if t2r < 8760*N
                                        Udcconv_p(ibr,iyr:t1r)=1;
                                        Udcconv_p(ibr,t2r)=0;
                                        in = floor(t2r/8760);
                                        iyr = t2r;
                                    else
                                        Udcconv_p(ibr,iyr:8760*N)=1;
                                        in = N;
                                    end
                                else
                                    dt = ceil(-MTTRsconv_dc(ibr)*log(r));
                                    t1r = iyr+dt-1;
                                    t2r = iyr+dt;
                                    if t2r < 8760*N
                                    Udcconv_p(ibr,iyr:t1r)=0;
                                    Udcconv_p(ibr,t2r)=1;
                                    in = floor(t2r/8760);
                                    iyr = t2r;
                                    else
                                    Udcconv_p(ibr,iyr:8760*N) = 0;
                                    in = N;
                                    end
                                end
                            end
                        end

                        for ibr = 1:size(Udcconv_n,1)
                            iyr=1;dt=1;in =1;
                            while in < N
                                r = rand;
                                if Udcconv_n(ibr,iyr) == 1
                                    dt = ceil(-MTTFsconv_dc(ibr)*log(r));
                                    t1r = iyr+dt-1;
                                    t2r = iyr+dt;

                                    if t2r< 8760*N
                                        Udcconv_n(ibr,iyr:t1r)=1;
                                        Udcconv_n(ibr,t2r)=0;
                                        in = floor(t2r/8760);
                                        iyr = t2r;
                                    else
                                        Udcconv_n(ibr,iyr:8760*N)=1;
                                        in=N;
                                    end

                                else
                                    dt = ceil(-MTTRsconv_dc(ibr)*log(r));
                                    t1r = iyr+dt-1;
                                    t2r = iyr+dt;
                                    if t2r < 8760*N
                                        Udcconv_n(ibr,iyr:t1r)=0;
                                        Udcconv_n(ibr,t2r)=1;
                                        in = floor(t2r/8760);
                                        iyr = t2r;
                                    else
                                        Udcconv_n(ibr,iyr:8760*N) = 0;
                                        in = N;
                                    end
                                end
                            end
                        end
                        else
                            for dcbr = 1:size(mpc.convdc,1)
                                iyr=1;dt=1;in=1;
                                while in < N
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
                                                    if t2r< 8760*N
                                                        Udcconv_p(dcbr,iyr:t1r)=1;
                                                        Udcconv_p(dcbr,t2r)=0;
                                                        in =floor(t2r/8760);
                                                        iyr=t2r;
                                                    else
                                                        Udcconv_p(dcbr,iyr:8760*N)=1;
                                                        in=N;
                                                    end
                                                else
                                                    if t2r< 8760*N
                                                        Udcconv_p(dcbr,iyr:t1r)=0;
                                                        Udcconv_p(dcbr,t2r)=1;
                                                        in =floor(t2r/8760);
                                                        iyr=t2r;
                                                    else
                                                        Udcconv_p(dcbr,iyr:8760*N)=0;
                                                        in=N;
                                                    end
                                                end
                                            case 2
                                                if trans_tipo(ipol)==0
                                                    if t2r< 8760*N
                                                        Udcconv_n(dcbr,iyr:t1r)=1;
                                                        Udcconv_n(dcbr,t2r)=0;
                                                        in =floor(t2r/8760);
                                                        iyr=t2r;
                                                    else
                                                        Udcconv_n(dcbr,iyr:8760*N)=1;
                                                        in=N;
                                                    end
                                                else
                                                    if t2r< 8760*N
                                                        Udcconv_n(dcbr,iyr:t1r)=0;
                                                        Udcconv_n(dcbr,t2r)=1;
                                                        in =floor(t2r/8760);
                                                        iyr=t2r;
                                                    else
                                                        Udcconv_n(dcbr,iyr:8760*N)=0;
                                                        in=N;
                                                    end
                                                end
                                            end
                                        end
                                        clear transition_rate StatusC_ini Vs transitio
                                    end
                                end
                            end
                            for i = 1:N
                                Udcconvs_n{i} = Udcconv_n(:,1+(i-1)*8760:8760*i);
                                Udcconvs_p{i} = Udcconv_p(:,1+(i-1)*8760:8760*i);
                            end

                            % PRINT CSV FILES FOR CSV TREE STRUCTURE
                            for in = 1:N
                                folderctgout = [folderout filesep 'YR_' num2str(in)];
                                mkdir(folderctgout);
                                mkdir([folderctgout filesep 'branch'])
                                csvwrite([folderctgout filesep 'branch' filesep 'br_status.csv'],[[1:size(mpc.branch,1)];Uaclines{in}']);
                                mkdir([folderctgout filesep 'gen'])
                                csvwrite([folderctgout filesep 'gen' filesep 'gen_status.csv'],[[1:size(mpc.gen,1)];Uacgens{in}']);
                                mkdir([folderctgout filesep 'branchdc'])
                                Udcbr_r{in}=Udclines_r{in}(idexdc(find(poli==0)),1:8760);
                                csvwrite([folderctgout filesep 'branchdc' filesep 'status_r.csv'],[[idexdc(find(poli==0))];Udcbr_r{in}'])
                                Udcbr_p{in}=Udclines_p{in}(idexdc(find(poli==1)),1:8760);
                                csvwrite([folderctgout filesep 'branchdc' filesep 'status_p.csv'],[[idexdc(find(poli==1))];Udcbr_p{in}'])
                                Udcbr_n{in}=Udclines_n{in}(idexdc(find(poli==2)),1:8760);
                                csvwrite([folderctgout filesep 'branchdc' filesep 'status_n.csv'],[[idexdc(find(poli==2))];Udcbr_n{in}'])
                                mkdir([folderctgout filesep 'convdc'])
                                csvwrite([folderctgout filesep 'convdc' filesep 'status_p.csv'],[[1:size(mpc.convdc,1)];Udcconvs_p{in}']);
                                csvwrite([folderctgout filesep 'convdc' filesep 'status_n.csv'],[[1:size(mpc.convdc,1)];Udcconvs_n{in}']);
                            end
                        end
                        %%% END - GENERATION of N-1 and HVDC unavailabilities
                        warning on
