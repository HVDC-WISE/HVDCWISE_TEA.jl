function [ ok TOLLE errorep LUNGH0 distanzaEstim mpc] = Modello_forze_matp( mpc,flag )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% try

% keyboard
% Tutte method
switch flag
    case 1
        
XL = mpc.branch(:,4).*(mpc.bus(mpc.branch(:,1),10).^2)./100;

W = 200;
L = 200;


area = W*L;
k = sqrt(area/size(mpc.bus,1));

% k = 1;%sqrt(area/length(nodo));

A = zeros(size(mpc.bus,1),size(mpc.bus,1));
X = rand(size(mpc.bus,1),2)*[W 0;0 L];
iterations0 = 5000;

coola = 37;
X0=X;
h = waitbar(0,'... using force oriented graph theory','Name','calculating geo localisation of grid nodes');

stima_X = 0.5; % 0.5 ohm per km per il 60% lunghezza

if max(mpc.bus(:,8))>220
    xdefault=0.5;
else
   xdefault=0.5; 
end

if mode(mpc.bus(:,8))<100
    LUNGH0 = (XL./xdefault);
else
    LUNGH0 = ceil(XL./xdefault);
end

i_trafi = find(mpc.branch(:,9)>0);
i_linee = find(mpc.branch(:,9)==0);

LUNGH0(i_trafi)=0.01;

Linee_lungh = LUNGH0(i_linee);

flag = 1;it=0;

for it =1:iterations0 
    TOL=0;%it=it+1;
    waitbar(it/iterations0,h);
    for in = 1:size(mpc.bus,1)
        L1=find(mpc.branch(i_linee,1)==in)';L2=find(mpc.branch(i_linee,2)==in)';
        linee_connesse = union(find(mpc.branch(i_linee,1)==in),find(mpc.branch(i_linee,2)==in));
        
        T1=find(mpc.branch(i_trafi,1)==in)';T2=find(mpc.branch(i_trafi,2)==in)';
        tra_connessi = union(find(mpc.branch(i_trafi,1)==in),find(mpc.branch(i_trafi,2)==in));
        
        estremi_lontaniT = [mpc.branch(i_trafi(T2),1)',mpc.branch(i_trafi(T1),2)'];
        estremi_lontani = [mpc.branch(i_linee(L2),1)',mpc.branch(i_linee(L1),2)'];
              
        vertice(in).disp = [0 0];
        
        gruppo = setdiff([1:size(mpc.bus,1)],in);
        
        kL = k*ones(1,length(gruppo));
        
        indici_li = []; indici_tra = [];
     
        for il = 1:length(estremi_lontani)
        indici_li(il) = find(ismember(gruppo,estremi_lontani(il)));
        end
        for itt = 1:length(estremi_lontaniT)
        indici_tra(itt) = find(ismember(gruppo,estremi_lontaniT(itt)));
        end
        try
        kL(indici_li) = [[Linee_lungh(L2)]',[Linee_lungh(L1)]'];
        kL(indici_tra) = 0.01*ones(1,length(indici_tra));
        catch err
            keyboard
        end
        delta(1,:) = X(in,1)-X(gruppo,1);
        delta(2,:) = X(in,2)-X(gruppo,2);
        
        absdelta = sqrt((delta(1,:)).^2+(delta(2,:)).^2);
        
        vertice(in).disp = (delta*(kL.^2./(absdelta.^2))')';
        clear delta absdelta
    end
    for il = 1:length(i_linee)
        
        es1 = mpc.branch(i_linee(il),1);
        es2 = mpc.branch(i_linee(il),2);
        A(es1,es2) = 1;
        delta0(1) = X(es2,1)-X(es1,1);
        delta0(2) = X(es2,2)-X(es1,2);
        k = Linee_lungh(il);
       vertice(es2).disp = vertice(es2).disp -((delta0/norm(delta0))*(norm(delta0))^2/k);
       vertice(es1).disp = vertice(es1).disp +((delta0/norm(delta0))*(norm(delta0))^2/k);
    TOL = TOL +   abs(norm(delta0)-k);
%        vertice(es2).disp = vertice(es2).disp +((delta0/norm(delta0))*k^2/(norm(delta0)));
%        vertice(es1).disp = vertice(es1).disp +((delta0/norm(delta0))*k^2/(norm(delta0)));
       
    end
    for itr = 1:length(i_trafi)
        es1 = mpc.branch(i_trafi(itr),1);
        es2 = mpc.branch(i_trafi(itr),2);
        A(es1,es2) = 1;
        delta0(1) = X(es2,1)-X(es1,1);
        delta0(2) = X(es2,2)-X(es1,2);
       vertice(es2).disp = vertice(es2).disp -((delta0/norm(delta0)).*(norm(delta0))^2/0.01);
       vertice(es1).disp = vertice(es1).disp +((delta0/norm(delta0)).*(norm(delta0))^2/0.01);
%         vertice(es2).disp = vertice(es2).disp +((delta0/norm(delta0))*0.1^2/(norm(delta0)));
%        vertice(es1).disp = vertice(es1).disp +((delta0/norm(delta0))*0.1^2/(norm(delta0)));
       TOL = TOL +   abs(norm(delta0)-0.01);
    end
    
    for in = 1:size(mpc.bus,1)
        X(in,:) = X(in,:) + vertice(in).disp.*(coola)/(norm(vertice(in).disp)); 
         
    end
    
    FORZA1(it) = norm(vertice(min(26,size(mpc.bus,1))).disp);
    FORZA2(it) = norm(vertice(min(26,size(mpc.bus,1))).disp);
    
    prova(it).MATR = X;
    coola = max(coola*(1-it/5e4),0.1);
    TOLLE(it)=TOL;
    
    
%     TOL
%     TOL
    
end
delete(h);
quale = find(TOLLE == min(TOLLE));
X = prova(quale).MATR;
% figure
% close all
% scatter(X(:,1),X(:,2),'bo'),hold on,scatter(X0(:,1),X0(:,2),'rs')
% 
% figure
% plot(FORZA1)
% hold on
% plot(FORZA2,'r')
for ilt = 1:length(i_trafi)
     idx1 = mpc.branch(i_trafi(ilt),1);
        idx2 = mpc.branch(i_trafi(ilt),2);
     Media =( X(idx1,:)+X(idx2,:))*0.5;
     
     X(idx1,:)= Media + [sign(X(idx1,1)-Media(1)) sign(X(idx1,2)-Media(2))].* [0.001 0.001];
     X(idx2,:)= Media + [sign(X(idx2,2)-Media(2)) sign(X(idx2,2)-Media(2))].* [0.001 0.001];
     
end


for il = 1:size(mpc.branch,1)
   beta_nl(il,mpc.branch(il,1))=1;
   beta_nl(il,mpc.branch(il,2))= 1;
   
end
% beta_nl_ = create_parameter('beta_nl',beta_nl,2);        %Upward generation cost of DGUs
% distanza_ = create_parameter('distanza',LUNGH0,1);
% 
% wgdx('input_grafo_rete.gdx',bu,l,distanza_,beta_nl_);
% cd ..
% file_gams = [pwd filesep 'GAMSforISAP' filesep 'costruisci_grafo_rete'];
% cd bin
% gams([file_gams ' lo=2']);    %rse.gms is the GAMS file with the optimization problem
% 
% outstri = 'grafo_coord.gdx';
% 
% Xs = read1p1('X',zeros(1,nr_buses),outstri);
% Ys = read1p1('Y',zeros(1,nr_buses),outstri);
% 
% figure,scatter(Xs,Ys)
nr_bus = size(beta_nl,2);
[zf fvalt] = fminunc(@(z)costruisci_grafo(z,LUNGH0,beta_nl),reshape(X,1,2*size(X,1)));

Xs=zf(1:nr_bus);
Ys = zf(1+nr_bus:end);

for ili = 1:size(mpc.branch,1)
   idx1 = mpc.branch(ili,1);
   idx2 = mpc.branch(ili,2);
    A(idx1,idx2)=1;
  %  A(idx2,idx1)=1;
   distanzaEstim(ili) = sqrt((Xs(idx1)-Xs(idx2))^2+(Ys(idx1)-Ys(idx2))^2);

   errore(ili) = distanzaEstim(ili)-LUNGH0(ili);
    errorep(ili) = (distanzaEstim(ili)-LUNGH0(ili))*100/LUNGH0(ili);
end

X = [Xs' Ys'];

figure
gplot(A,X)

fid = fopen('testCoo.txt','wt+');
for inn = 1:size(mpc.bus,1)
   %nodo(inn).coord = X(inn,:);
   fprintf(fid,[num2str(inn) ' ' num2str(X(inn,1)) ' ' num2str(X(inn,2)) '\n']);
   hold on
   text(X(inn,1),X(inn,2),num2str(inn))
end
fclose(fid);
save coordProva.mat X
mpc.coord = [X];
    case 2
   XL = mpc.branchdc(:,3).*(mpc.busdc(mpc.branchdc(:,1),5).^2)./100;

W = 200;
L = 200;


area = W*L;
k = sqrt(area/size(mpc.busdc,1));

% k = 1;%sqrt(area/length(nodo));

A = zeros(size(mpc.busdc,1),size(mpc.busdc,1));
X = rand(size(mpc.busdc,1),2)*[W 0;0 L];
iterations0 = 5000;

coola = 37;
X0=X;
h = waitbar(0,'... using force oriented graph theory','Name','calculating geo localisation of grid nodes');

stima_X = 0.5; % 0.5 ohm per km per il 60% lunghezza

if max(mpc.busdc(:,5))>220
    xdefault=0.5;
else
   xdefault=0.5; 
end

if mode(mpc.busdc(:,5))<100
    LUNGH0 = (XL./xdefault);
else
    LUNGH0 = ceil(XL./xdefault);
end

% i_trafi = find(mpc.branch(:,9)>0);
i_linee = [1:size(mpc.branchdc,1)];
i_trafi=[];
LUNGH0(i_linee)=0.01;

Linee_lungh = LUNGH0(i_linee);

flag = 1;it=0;

for it =1:iterations0 
    TOL=0;%it=it+1;
    waitbar(it/iterations0,h);
    for in = 1:size(mpc.busdc,1)
        L1=find(mpc.branchdc(i_linee,1)==in)';L2=find(mpc.branchdc(i_linee,2)==in)';
        linee_connesse = union(find(mpc.branchdc(i_linee,1)==in),find(mpc.branchdc(i_linee,2)==in));
        
        T1=find(mpc.branchdc(i_trafi,1)==in)';T2=find(mpc.branchdc(i_trafi,2)==in)';
        tra_connessi = union(find(mpc.branchdc(i_trafi,1)==in),find(mpc.branchdc(i_trafi,2)==in));
        
        estremi_lontaniT = [mpc.branchdc(i_trafi(T2),1)',mpc.branchdc(i_trafi(T1),2)'];
        estremi_lontani = [mpc.branchdc(i_linee(L2),1)',mpc.branchdc(i_linee(L1),2)'];
              
        vertice(in).disp = [0 0];
        
        gruppo = setdiff([1:size(mpc.busdc,1)],in);
        
        kL = k*ones(1,length(gruppo));
        
        indici_li = []; indici_tra = [];
     
        for il = 1:length(estremi_lontani)
        indici_li(il) = find(ismember(gruppo,estremi_lontani(il)));
        end
        for itt = 1:length(estremi_lontaniT)
        indici_tra(itt) = find(ismember(gruppo,estremi_lontaniT(itt)));
        end
        try
        kL(indici_li) = [[Linee_lungh(L2)]',[Linee_lungh(L1)]'];
        kL(indici_tra) = 0.01*ones(1,length(indici_tra));
        catch err
            keyboard
        end
        delta(1,:) = X(in,1)-X(gruppo,1);
        delta(2,:) = X(in,2)-X(gruppo,2);
        
        absdelta = sqrt((delta(1,:)).^2+(delta(2,:)).^2);
        
        vertice(in).disp = (delta*(kL.^2./(absdelta.^2))')';
        clear delta absdelta
    end
    for il = 1:length(i_linee)
        
        es1 = mpc.branchdc(i_linee(il),1);
        es2 = mpc.branchdc(i_linee(il),2);
        A(es1,es2) = 1;
        delta0(1) = X(es2,1)-X(es1,1);
        delta0(2) = X(es2,2)-X(es1,2);
        k = Linee_lungh(il);
       vertice(es2).disp = vertice(es2).disp -((delta0/norm(delta0))*(norm(delta0))^2/k);
       vertice(es1).disp = vertice(es1).disp +((delta0/norm(delta0))*(norm(delta0))^2/k);
    TOL = TOL +   abs(norm(delta0)-k);
%        vertice(es2).disp = vertice(es2).disp +((delta0/norm(delta0))*k^2/(norm(delta0)));
%        vertice(es1).disp = vertice(es1).disp +((delta0/norm(delta0))*k^2/(norm(delta0)));
       
    end
    for itr = 1:length(i_trafi)
        es1 = mpc.branchdc(i_trafi(itr),1);
        es2 = mpc.branchdc(i_trafi(itr),2);
        A(es1,es2) = 1;
        delta0(1) = X(es2,1)-X(es1,1);
        delta0(2) = X(es2,2)-X(es1,2);
       vertice(es2).disp = vertice(es2).disp -((delta0/norm(delta0)).*(norm(delta0))^2/0.01);
       vertice(es1).disp = vertice(es1).disp +((delta0/norm(delta0)).*(norm(delta0))^2/0.01);
%         vertice(es2).disp = vertice(es2).disp +((delta0/norm(delta0))*0.1^2/(norm(delta0)));
%        vertice(es1).disp = vertice(es1).disp +((delta0/norm(delta0))*0.1^2/(norm(delta0)));
       TOL = TOL +   abs(norm(delta0)-0.01);
    end
    
    for in = 1:size(mpc.busdc,1)
        X(in,:) = X(in,:) + vertice(in).disp.*(coola)/(norm(vertice(in).disp)); 
         
    end
    
    FORZA1(it) = norm(vertice(min(26,size(mpc.busdc,1))).disp);
    FORZA2(it) = norm(vertice(min(26,size(mpc.busdc,1))).disp);
    
    prova(it).MATR = X;
    coola = max(coola*(1-it/5e4),0.1);
    TOLLE(it)=TOL;
    
    
%     TOL
%     TOL
    
end
delete(h);
quale = find(TOLLE == min(TOLLE));
X = prova(quale).MATR;
% figure
% close all
% scatter(X(:,1),X(:,2),'bo'),hold on,scatter(X0(:,1),X0(:,2),'rs')
% 
% figure
% plot(FORZA1)
% hold on
% plot(FORZA2,'r')
for ilt = 1:length(i_trafi)
     idx1 = mpc.branchdc(i_trafi(ilt),1);
        idx2 = mpc.branchdc(i_trafi(ilt),2);
     Media =( X(idx1,:)+X(idx2,:))*0.5;
     
     X(idx1,:)= Media + [sign(X(idx1,1)-Media(1)) sign(X(idx1,2)-Media(2))].* [0.001 0.001];
     X(idx2,:)= Media + [sign(X(idx2,2)-Media(2)) sign(X(idx2,2)-Media(2))].* [0.001 0.001];
     
end

for il = 1:size(mpc.branchdc,1)
   beta_nl(il,mpc.branchdc(il,1))=1;
   beta_nl(il,mpc.branchdc(il,2))= 1;
   
end
% beta_nl_ = create_parameter('beta_nl',beta_nl,2);        %Upward generation cost of DGUs
% distanza_ = create_parameter('distanza',LUNGH0,1);
% 
% wgdx('input_grafo_rete.gdx',bu,l,distanza_,beta_nl_);
% cd ..
% file_gams = [pwd filesep 'GAMSforISAP' filesep 'costruisci_grafo_rete'];
% cd bin
% gams([file_gams ' lo=2']);    %rse.gms is the GAMS file with the optimization problem
% 
% outstri = 'grafo_coord.gdx';
% 
% Xs = read1p1('X',zeros(1,nr_buses),outstri);
% Ys = read1p1('Y',zeros(1,nr_buses),outstri);
% 
% figure,scatter(Xs,Ys)
nr_bus = size(beta_nl,2);
[zf fvalt] = fminunc(@(z)costruisci_grafo(z,LUNGH0,beta_nl),reshape(X,1,2*size(X,1)));

Xs=zf(1:nr_bus);
Ys = zf(1+nr_bus:end);

for ili = 1:size(mpc.branchdc,1)
   idx1 = mpc.branchdc(ili,1);
   idx2 = mpc.branchdc(ili,2);
    A(idx1,idx2)=1;
  %  A(idx2,idx1)=1;
   distanzaEstim(ili) = sqrt((Xs(idx1)-Xs(idx2))^2+(Ys(idx1)-Ys(idx2))^2);

   errore(ili) = distanzaEstim(ili)-LUNGH0(ili);
    errorep(ili) = (distanzaEstim(ili)-LUNGH0(ili))*100/LUNGH0(ili);
end

X = [Xs' Ys'];

figure
gplot(A,X)

fid = fopen('testCoo.txt','wt+');
for inn = 1:size(mpc.busdc,1)
   %nodo(inn).coord = X(inn,:);
   fprintf(fid,[num2str(inn) ' ' num2str(X(inn,1)) ' ' num2str(X(inn,2)) '\n']);
   hold on
   text(X(inn,1),X(inn,2),num2str(inn))
end
fclose(fid);
save coordProva.mat X
mpc.busdc_coord = [X]
     
        
end
% %Eades method
% return
% A = zeros(length(nodo),length(nodo));
% for jtempt = 1:Ntempt
%     TOLT=0;
% for in = 1:length(nodo)
%     L1=find(mpc.branch(:,1)==in);L2=find(mpc.branch(:,2)==in);
%    linee_connesse = union(find(mpc.branch(:,1)==in),find(mpc.branch(:,2)==in)); 
%    
%     T1=find([trafo2.estremo1_ID]==in);T2=find([trafo2.estremo2_ID]==in);
%    tra_connessi = union(find([trafo2.estremo1_ID]==in),find([trafo2.estremo2_ID]==in)); 
%    
%       estremi_lontaniT = union([trafo2.estremo1_ID(T2)],trafo2.estremo2_ID(T1));
%    estremi_lontani = union([linea.estremo1_ID(L2)],[linea.estremo2_ID(L1)]);
%    A(in,estremi_lontani) = 1;
%    A(in,estremi_lontaniT) = 1;
%    for is = 1:length(estremi_lontani)
%        TOLT = TOLT + abs(norm(abs(X(in,:)-X(estremi_lontani(is),:)))-linea.lungh(linee_connesse(is)));
%        
%        Fsx2(is) = -ks * (norm(abs(X(in,:)-X(estremi_lontani(is),:)))-linea.lungh(linee_connesse(is))) * (X(in,1)-X(estremi_lontani(is),1))/(norm(abs(X(in,:)-X(estremi_lontani(is),:))));
%        Fsy2(is) = -ks * (norm(abs(X(in,:)-X(estremi_lontani(is),:)))-linea.lungh(linee_connesse(is))) * (X(in,2)-X(estremi_lontani(is),2))/(norm(abs(X(in,:)-X(estremi_lontani(is),:))));
%    end
%    
%    for is = 1:length(estremi_lontaniT)
%        TOLT = TOLT + abs(norm(abs(X(in,:)-X(estremi_lontaniT(is),:)))-0.2);
%        Fsxt2(is) = -kst * (norm(abs(X(in,:)-X(estremi_lontaniT(is),:)))-0.2) * (X(in,1)-X(estremi_lontaniT(is),1))/(norm(abs(X(in,:)-X(estremi_lontaniT(is),:))));
%        Fsyt2(is) = -kst * (norm(abs(X(in,:)-X(estremi_lontaniT(is),:)))-0.2) * (X(in,2)-X(estremi_lontaniT(is),2))/(norm(abs(X(in,:)-X(estremi_lontaniT(is),:))));
%    end
%         
%            gruppo = setdiff([1:length(nodo)],[estremi_lontani estremi_lontaniT]);
%            d = sqrt((X(in,1)-X(setdiff(gruppo,in),1)).^2+(X(in,2)-X(setdiff(gruppo,in),2)).^2)';
%           Frx2 = kr./(d.^2).*(X(in,1)-X(setdiff(gruppo,in),1))'./d;
%           Fry2 = kr./(d.^2).*(X(in,2)-X(setdiff(gruppo,in),2))'./d;
%           
%       Fx2(in) = sum(Fsx2)+sum(Frx2)+sum(Fsxt2);
%       Fy2(in) = sum(Fsy2)+sum(Fry2)+sum(Fsyt2);
%       
%     nodo(in).coord = [X(in,1) X(in,2)];  
% %     keyboard  
%     Fsx=[];Fsy=[];
% Fsxt=[];Fsyt=[];    
% end
%     verifica(jtempt).MATR = X;
%     
% %     X(:,1) = X(:,1) + c4*Fx';
% %     X(:,2) = X(:,2) + c4*Fy';
%      X(:,1) = X(:,1) + dt*Fx2';
%     X(:,2) = X(:,2) + dt*Fy2';
%     TOL(jtempt) = TOLT;%sum(sqrt(Fx.^2+Fy.^2));
%     TOL(end)
%     
%     if TOL(jtempt) < thres
%         break
%     end
%     
% end
% 
% quale = find(TOL==min(TOL));
% 
% X = verifica(quale).MATR;
% close all



ok=1;

function FO = costruisci_grafo(z,distanza,beta_nl)
        FO=0;
        nr_b = size(beta_nl,2);
        x = z(1:nr_b);
        y = z(1+nr_b:end);
        for li = 1:length(distanza)
            quali = find(beta_nl(li,:)>0);
            FO = FO+abs(distanza(li) - norm([diff(x(quali)) diff(y(quali))]));
        end
    end

end

