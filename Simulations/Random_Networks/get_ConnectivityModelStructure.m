function Model = get_ConnectivityModelStructure(N, density, flagConnectedNodes)

%%%% input:
% N:                  Nodes Number to obtain a NxN upper triangular matrix
% density:            connections density selected for the model. The maximum
%                     possible density is 0,50.
% flagConnectedNodes: --> 0: some nodes are disconnected
%                     --> 1: all the nodes are connected  
%%%% output:
% Model-->            Connectivity Model Structure. It is a NxN matrix in which 
%                     entries are equal to 1 in correspondence to a link

% Author: Manuela Petti
% Date: 05/06/2014
% version: 1.0.0

RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock))); %cosa rappresenta questo 100*clock?

if density==0.5
    Model = triu(ones(N,N),+1); %   TRIU(X) is the upper triangular part of X. xchè ci serve la parte triangolare superiore. 
    %questo dipende dal fatto che avendo una sorgente imposta cioè un canale dell'EEG reale 
    %posso conoscere solo i legami di causalità da x1 verso gli altri nodi e non viceversa xchè non conosco il passato di x1
                                 
                                  
   return
end

tic
while ~exist('Model')
    
    t = toc;
    if t > 300
        error('The search for the right combination of connections for selected parameters (N, density) requires more than 5 minutes')
    end
    
    Model = triu(ones(N,N),+1)*2;%senza il +1 metterebbe degli elementi anche sulla diagonale principare che a noi nn interessa e dunque viene  fuori una matrice triangolare superiore di tutti 2 tranne che sulla diagonale principale e sotto la diagonale suddetta
    ind = find(Model); %I = FIND(X) returns the linear indices corresponding to 
                       %   the nonzero entries of the array X. le posizioni
                       %   degli elementi diversi da zero nella matrice del
                       %   modello
    N_conn = round(density*N*(N-1));%ci dà il numero delle connessioni nella rete in funzione del valore di densità impostata 
    
    if (N_conn < N-1) && flagConnectedNodes
        error('The selected density does not allow to connect all the nodes')
    end

    mixed = randperm(length(ind)); %P = RANDPERM(N) returns a vector containing a random permutation of the
                                    %integers 1:N, dunque fa una
                                    %permutazione random di tutti gli
                                    %elementi diversi da 0 del modello
    selected = ind(mixed(1:N_conn));
    Model(selected) = 1; Model(Model==2) = 0;
    
    D = distance_bin(Model);                  %calcola la matrice delle distanze (collegamento più breve tra due nodi) a partire dalla matrice di adiacenza per ogni coppia di nodi
    indInf = find(D(1,:)==Inf);               %find restituisce le posizioni di una matrice o vettore in cui i valori sono diversi da zero 
    if flagConnectedNodes && ~isempty(indInf) %isemmpty restituisce 0 od 1 a seconda se l'array è vuoto o meno
        clear Model
        continue
    end
    
    if ~isempty(indInf)
        for ii=1:length(indInf)
            wrongConn = find(Model(indInf(ii),:));
            if ~isempty(wrongConn)
                for wC=1:length(wrongConn)
                    Model(indInf(ii),wrongConn(wC)) = 0;
                end
            end
        end
    end
    
    updated_N_conn = length(find(Model));
    
    if updated_N_conn == N_conn
        break
   
    elseif updated_N_conn < N_conn
        existingConn = length(find(Model));
        missingConn = N_conn-existingConn;
        N_ImpLnks = 0;
        for ii=1:length(indInf)
            N_ImpLnks = N_ImpLnks+N-indInf(ii);
        end
        if length(ind)-N_ImpLnks-existingConn < missingConn
            clear Model
            continue
        end
        Pos_ImpLnks = [];
        for ii=1:length(indInf)
            if indInf(ii) == N
                break
            end
            start = indInf(ii)+1;
            for nn=start:N
                Pos_ImpLnks = [Pos_ImpLnks (N*indInf(ii))+indInf(ii)+N*(nn-start)];
            end
        end
        Pos_extLinks = find(Model);
        NewPossible = setdiff(ind,union(Pos_ImpLnks, Pos_extLinks));
        
        New_mixed = randperm(length(NewPossible));
        newselected = NewPossible(New_mixed(1:missingConn));
        Model(newselected) = 1;
        
        if length(find(Model)) == N_conn
            break
        end
        
    end
    
end


