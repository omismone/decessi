clc
clear
close all

%% carico i decessi
decessi=readtable('iss_bydate_italia_deceduti.csv');
giorniD=table2array(decessi(225:375,1));    %dati ott2020-feb2021 nelle righe 225-375
casiD=table2array(decessi(225:375,2));
media7ggD=table2array(decessi(225:375,3));
%media7ggD_sfasati_D=table2array(decessi((225+D):(375+D),3));

%% plotto i decessi
%{
figure(1);
subplot(2,1,1);
plot(giorniD,casiD);
title('casi puntuali decessi');
xlabel('data');
ylabel('decessi');
grid on
hold on
subplot(2,1,2);
plot(giorniD,media7ggD);
title('media 7 giorni dei decessi');
xlabel('data');
ylabel('media decessi');
grid on
%}
%% carico i positivi
positivi=readtable('iss_bydate_italia_positivi.csv');
giorniP=table2array(positivi(247:397,1));   %dati ott2020-feb2021 nelle righe 247-397
casiP=table2array(positivi(247:397,2));
media7ggP=table2array(positivi(247:397,3));
%% stimo D a caso
D=15;                                       %ritardo dei morti rispetto ai rispettivi positivi
media7ggP_sfasati_D=table2array(positivi((247-D):(397-D),3));

%% plotto prime figure
%{
figure(2);
subplot(2,1,1);
plot(giorniP,casiP);
title('casi puntuali positivi');
xlabel('data');
ylabel('positivi');
grid on
hold on
subplot(2,1,2);
plot(giorniP,media7ggP);
title('media 7 giorni dei positivi');
xlabel('data');
ylabel('media positivi');
grid on
%}
%% fattore scala con lscov
phi=media7ggP_sfasati_D;
fattore_scala=lscov(phi,media7ggD);     %x far combaciare ampiezze

%% plotto assieme positivi e decessi (D=15)
figure(3);
plot(giorniP,media7ggP_sfasati_D.*fattore_scala);
grid on
hold on
plot(giorniD,media7ggD,'x','Color','r');
title('modello con guadagno e ritardo');
legend('positivi scalati e moltiplicati','decessi');
ssr2=(media7ggD - media7ggP_sfasati_D.*fattore_scala)' *(media7ggD - media7ggP_sfasati_D.*fattore_scala);

%% modello esponenziale con somma
%{
sommaExp=0;
lambda=0.091;
D=7;                                    %andrebbe modificato il valore di D (ssr migliore con D=7)
%fattore_scala=0.03;                    %va modificato anche lui?
for k=0:150
    ingressi=table2array(positivi((247-D-k):(397-D-k),3));
    g=(fattore_scala*exp(-lambda*k))*lambda;
    sommaExp=sommaExp+g.*ingressi;
end
%}
%% plotto modello exp con decessi
%{
figure(4);
plot(giorniP,sommaExp);
grid on
hold on
plot(giorniD,media7ggD,'x','Color','r');
title('modello esponenziale con sommatoria');
legend('positivi x formula esponenziale','decessi');
ssr3=(media7ggD - sommaExp)' * (media7ggD - sommaExp);
%}

%% trovo D e lambda migliori modello exp
D_vett=linspace(0,19,20);
l_vett=linspace(0,1,100);

%faccio una sorta di prodotto cartesiano x trovare ssr minore
for i=1:20
    for j=1:100
        ssr_matrix(i,j)=calculateSSR(positivi,D_vett(i),fattore_scala,l_vett(j),media7ggD);
    end
end

ssr_min = min(min(ssr_matrix));                                     %minimo trovato in matrice
[indice_D_min, indice_lambda_min] = find(ssr_matrix == ssr_min);    %indici x cui è min
D=D_vett(indice_D_min);                                             %valori x quegli indici
lambda=l_vett(indice_lambda_min);

%% convoluzione
ingressi=table2array(positivi((247):(397),3));
conv=convoluz(ingressi,D,fattore_scala,lambda);
figure(5)
plot(giorniP,conv);
hold on
grid on
plot(giorniD,media7ggD,'x','Color','r');
title('modello esponenziale con convoluzione');
legend('positivi x formula esponenziale','decessi');
%ssr è ssr_min
