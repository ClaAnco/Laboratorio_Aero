%% traccia Hess Smith (2024)

clc
close all
clear 

addpath mat_functions

%% Input

U_inf = 1;  % Velocità all'infinito [m/s]
alpha = 1;   % Angolo di attacco [°]
U_inf_x = U_inf * cos(deg2rad(alpha));
U_inf_y = U_inf * sin(deg2rad(alpha));

U_inf = [U_inf_x; U_inf_y];
U_inf_normal = [-U_inf(2); U_inf(1)];
U_inf_normal = U_inf_normal ./ norm(U_inf_normal);

TestCase = 0;

CodiceProfilo = '0012';
Chord = 1;
NPannelli = 101;

LE_X_Position = 0;
LE_Y_Position = 0;

%% Creazione profilo

% numero profilo:
% [x,y]=createProfile(CodiceProfilo,NPannelli,Chord);

Corpo = importXfoilProfile(strcat('NACA_', CodiceProfilo, '.dat')); %strcat concatena le stringe in modo da ottenere il nome del file richiamato
% Prima flippa i vettori
x = flipud(Corpo.x);
y = flipud(Corpo.y);
Corpo.x = x.*Chord;
Corpo.y = y.*Chord;

figure;
plot(x, y, 'o-')
axis equal

%% Creazione di una struttura di pannelli

[Centro, Normale, Tangente, Estremo_1, Estremo_2, alpha, lunghezza, L2G_TransfMatrix, G2L_TransfMatrix] = CreaStrutturaPannelli(Corpo);
        
%% Inizializzazione matrici e vettori

% Ora che ho i pannelli, posso inizializzare la matrice ed i vettori

NCols = sum(NPannelli) + 1;
NRows = NCols;
matriceA = zeros(NRows, NCols);
TermineNoto = zeros(NRows, 1);

%% Creazione della matrice quadrata As

%Inizializzo i vettori delle velocità da usare nel calcolo del Cp
ViSorgente_Cp = [];
ViVortice_Cp = []   ;

for i = 1:NPannelli
    index_i = i; % riga

    Centro_qui = Centro(i, :)';
    Normale_qui = Normale(i, :)';

    indexStart_colonna = 0;

        for j = 1:NPannelli
            index_j = indexStart_colonna + j;  % Colonna

            Estremo_1_qui = Estremo_1(j, :)';
            Estremo_2_qui = Estremo_2(j, :)';

            L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(j, :, :));
            G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(j, :, :));
            
            %aggiungo un modo di popolare due matrici con tutte le velocità
            %indotte di sorgenti e vortici, mi serviranno più avanti per il
            %calcolo del Cp. Tali matrici contengono le due componenti di
            %ognuno dei 101 contributi di velocità sul pannello i
            ViSorgente_Cp = [ViSorgente_Cp; (ViSorgente(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui)')];
                
            ViVortice_Cp = [ViVortice_Cp; (ViVortice(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui)')];
            
            
            matriceA(index_i, index_j) = dot(ViSorgente(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Normale_qui);

            matriceA(index_i, sum(NPannelli)+1) = matriceA(index_i, sum(NPannelli)+1) + dot(ViVortice(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Normale_qui);


        end

end


%% Creazione delle componenti dei vettori a_v, c_s e c_v


Centro_Start = Centro(1, :)';
Tangente_Start = Tangente(1, :)';

Centro_End = Centro(end, :)';
Tangente_End = Tangente(end, :)';


b = 0;
for j = 1:NPannelli(1)

    index_j = j;

    Estremo_1_qui = Estremo_1(j, :)';
    Estremo_2_qui = Estremo_2(j, :)';
    L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(j, :, :));
    G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(j, :, :));

    a = dot(ViSorgente(Centro_Start, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_Start);
    b = b + dot(ViVortice(Centro_Start, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_Start);

    a = a + dot(ViSorgente(Centro_End, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_End);
    b = b + dot(ViVortice(Centro_End, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_End);


    matriceA(sum(NPannelli) + 1, index_j) = a;

end

matriceA(sum(NPannelli) + 1, sum(NPannelli) + 1) = b;



%% Creazione del termine noto

for j = 1:NPannelli

    Normale_qui = Normale(j, :)';

    index = j;

    TermineNoto(index) = - dot(U_inf, Normale_qui);
end

Tangente_1 = Tangente(1, :)';
Tangente_end = Tangente(end, :)';
TermineNoto(sum(NPannelli) + 1) = - dot(U_inf, (Tangente_1 + Tangente_end));

%% Risoluzione sistema lineare
Soluzione = linsolve(matriceA,TermineNoto);


%% Calcolo del cp e della velocità sui pannelli
%le computazioni seguenti riproducono i calcoli suggeriti nella slide 20
SumSorgenti = zeros(1, 2); %vettore riga
SumVortici = zeros(1, 2); %vettore riga
gamma = Soluzione(102);
Cp = [];
for i = 1:NPannelli
    for j = 1 : NPannelli
        Componente_Sommatoria = ViSorgente_Cp(j, :)*Soluzione(j); %vettore riga
        Sum_sorgenti = SumSorgenti + Componente_Sommatoria; %vettore riga %sommatoria sorgenti
        
        SumVortici = SumVortici + ViVortice_Cp(j, :); %vettore riga %sommatoria vortici
    end
    
    Contr_Pann_i = U_inf + SumSorgenti' + (SumVortici*gamma)'; %vettore colonna
    Scal_tau_i = dot(Contr_Pann_i, (Tangente(i, :)));
   
    Cp = [Cp; 1 - ((Scal_tau_i)/(norm(U_inf)))^2];
end

%plot dei Cp su dorso e ventre
Cp_ventre = flipud(Cp(1:51));
Cp_dorso = Cp(51:101);
%x1 = 0:Chord/50:1;
%y1 = Cp;
y1 = Cp_ventre;
y2 = Cp_dorso;

hold on
plot(flipud(Centro(1:51, 1)),y1, Centro(51:101, 1),y2)
%set(gca,'Ydir','reverse')


%figure;
%plot(Centro(:,1),Cp,'r')



%%


