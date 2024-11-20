%% traccia Hess Smith (2024)

clc
close all
clear 
format long

addpath mat_functions


%% Input

U_inf = 1;  % Velocità all'infinito [m/s]
alpha = 1;   % Angolo di incidenza [°]
U_inf_x = U_inf * cos(deg2rad(alpha));  %scomposiozione della velocità nelle due componenti
U_inf_y = U_inf * sin(deg2rad(alpha));

U_inf = [U_inf_x; U_inf_y];
U_inf_normal = [-U_inf(2); U_inf(1)];  %definizione della velocità indisturbata normale 
U_inf_normal = U_inf_normal ./ norm(U_inf_normal);   %versore 


TestCase = 0;

CodiceProfilo = '0012';
Chord = 1;
NPannelli = 101;

LE_X_Position = 0;  %SDR centrato nel LE del profilo
LE_Y_Position = 0;

%% Creazione profilo

% numero profilo:
% [x,y]=createProfile(CodiceProfilo,NPannelli,Chord);

Corpo = importXfoilProfile(strcat('NACA_', CodiceProfilo, '.dat')); %crea profilo NACA scelto tramite una funzione che usa Xfoil
% Prima flippa i vettori
x = flipud(Corpo.x);
y = flipud(Corpo.y);
Corpo.x = x.*Chord;  %crea il profilo alare (due vettori che contengono le coordinate)
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
for i = 1:NPannelli
    index_i = i;  %riga

    Centro_qui = Centro(i, :)';
    Normale_qui = Normale(i, :)';

    indexStart_colonna = 0;

        for j = 1:NPannelli
            index_j = indexStart_colonna + j;  % Colonna

            Estremo_1_qui = Estremo_1(j, :)';
            Estremo_2_qui = Estremo_2(j, :)';

            L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(j, :, :));
            G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(j, :, :));

            matriceA(index_i, index_j) = dot(ViSorgente(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Normale_qui);

            matriceA(index_i, sum(NPannelli)+1) = matriceA(index_i, sum(NPannelli)+1) + dot(ViVortice(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Normale_qui);
           

        end

end

%% CALCOLO U_SORGENTI UNITARIO NEL RIFERIMENTO GLOBALE
matrice_u_s = zeros(NPannelli,2);

for i=1:NPannelli

    Centro_qui = Centro(i, :)';
    Estremo_1_qui = Estremo_1(i, :)';
    Estremo_2_qui = Estremo_2(i, :)';
    L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(i, :, :));
    G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(i, :, :));
    
    matrice_u_s(i,:) = ViSorgente(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui);
end


%% Creazione delle componenti dei vettori a_v, c_s e c_v


Centro_Start = Centro(1, :)';
Tangente_Start = Tangente(1, :)';

Centro_End = Centro(end, :)';
Tangente_End = Tangente(end, :)';

b = 0;
for j = 1:NPannelli(1)
    index_j = j;

    
    Centro_qui = Centro(j, :)';
    Estremo_1_qui = Estremo_1(j, :)';
    Estremo_2_qui = Estremo_2(j, :)';
    L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(j, :, :));
    G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(j, :, :));

    a = dot(ViSorgente(Centro_Start, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_Start);
    b = b + dot(ViVortice(Centro_Start, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_Start);

    a = a + dot(ViSorgente(Centro_End, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_End);
    b = b + dot(ViVortice(Centro_End, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_End);


    matriceA(sum(NPannelli) + 1, index_j) = a;


    %matrice_u_v(index_j,:)=ViVortice(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui);
end

matriceA(sum(NPannelli) + 1, sum(NPannelli) + 1) = b;

%% CALCOLO U_VORTICE UNITARIO NEL RIF GLOBALE
matrice_u_v=zeros(NPannelli,2);
for i=1:NPannelli

    Centro_qui = Centro(i, :)';
    Estremo_1_qui = Estremo_1(i, :)';
    Estremo_2_qui = Estremo_2(i, :)';
    L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(i, :, :));
    G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(i, :, :));

    matrice_u_v(i,:)=ViVortice(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui);
end




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
%Soluzione2 = matriceA\TermineNoto;






%% CALCOLO U_INF TANGENTE IN RIF GLOBALE

U_inf_tang = zeros(NPannelli,1);

for i=1:NPannelli
    %U_inf(i,:)=[U_inf(1),U_inf(2)];
    U_inf_tang(i) = U_inf(1)*Tangente(i,1) + U_inf(2)*Tangente(i,2);
end
%U_inf_tang = U_inf(:,1).*Tangente(:,1)+U_inf(:,2).*Tangente(:,2);

%% CALCOLO U_SORGENTE TANGENTE IN RIF GLOBALE

u_s=zeros(NPannelli,2);  %definisco velocità indotta da sorgente
for i=1:NPannelli
    u_s(i,:) =Soluzione(i)*matrice_u_s(i,:);
end



%% CALCOLO U_VORTICE TANGENTE IN RIF GLOBALE

u_v = matrice_u_v .* Soluzione(NPannelli+1); %velocità indotta da vortice

u_v_tang = zeros(NPannelli,1);


%% CALCOLO IL CAMPO DI VELOCITA' TOTALE INTORNO AL PROFILO PER OGNI PANNELLO
U_tot=zeros(NPannelli,2);
U_infinito=zeros(101,2);
for i=1:NPannelli
    U_infinito(i,:)=[U_inf(1),U_inf(2)];
end
U_tot= U_infinito + u_s + u_v;

%% CALCOLO VELOCITA' TANGENTE
V_tang=zeros(NPannelli,1);
for i=1:NPannelli
V_tang(i)= U_tot(i,1)*Tangente(i,1)+U_tot(i,2)*Tangente(i,2);
end
%% DEFINISCO SIGMA E GAMMA
sigma=Soluzione(1:NPannelli);
gamma=Soluzione(NPannelli+1);


%% CALCOLO VELOCITA' PANNELLI 2 POSSIBILITA'
V_x=zeros(NPannelli,1);
V_y=zeros(NPannelli,1);


for i = 1:NPannelli
   
    index_i = i;  %riga
    V_x(i)=U_inf(1);
    V_y(i)=U_inf(2);

    Centro_qui = Centro(i, :)';
    Normale_qui = Normale(i, :)';

    indexStart_colonna = 0;

        for j = 1:NPannelli
            index_j = indexStart_colonna + j;  % Colonna

            Estremo_1_qui = Estremo_1(j, :)';
            Estremo_2_qui = Estremo_2(j, :)';

            L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(j, :, :));
            G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(j, :, :));

            Sorgente_unitario=ViSorgente(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui);
            Vortice_unitario=ViVortice(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui);

            V_x(i)= V_x(i)+sigma(j)*Sorgente_unitario(1,1)+gamma*Vortice_unitario(1,1);
            V_y(i)= V_y(i)+sigma(j)*Sorgente_unitario(2,1)+gamma*Vortice_unitario(2,1);
        end
end


V_t = V_x.*Tangente(:,1) + V_y.*Tangente(:,2);


%% CALCOLO Cp SECONDA POSSIBILITA'
Cp_2=1-V_t.^2/(U_inf(1)^2+U_inf(2)^2);


%% CALCOLO Cp
Cp = zeros(NPannelli,1);
for i=1:NPannelli
    %Cp(i) = 1 - ((U_inf_tang(i) + u_s_tang(i) + u_v_tang(i))^2)/ (norm(U_inf))^2;
    Cp(i)= 1-(V_tang(i)^2)/(norm(U_inf))^2;

end 
hold on
%plot(Centro(1:50,:),Cp(1:50),'r');
%hold on
%plot(Centro(51:end,:),Cp(51:end),'g')
%set(gca,'Ydir','reverse')

plot(Centro(:,1),Cp_2,'r')
set(gca,'Ydir','reverse')





