clearvars -except tExp ApExp
%Area calculation included
%clustering included

%%  Revisar siempre que el ajuste de parámetros esté bien

%% Parameter definition -----------------------------------------------------------------------
%Simulation parameters
dt=0.00002;
tmax=2000;
%n=0;
cont = 0;
t1=300;
%t2=600;
areat=[];
Lt=[];
clt=[];

vc=200; %Niveles de clustering o número de células
CM=jet(vc);

for i=1:vc
    
%System parameters
kdm = 0.00015;%para una adaptacion de 200 s mientras mas pequeno, mas se demora el sistema en adaptar.
km = 0.05;
kbn = 0.008;
kunb = 150;%para que haya unas 20000 ligaciones y deligaciones por s.

kb=0.001;
pb=1;
k2=8*0.9;%constante de fosfotilacion de cheA por los receptores activos (#receptores-1*s-1)
k1=0.05*0.9;%constante pra los otros tipos de receptores
k3=0.001*0.9;%constante para los receptores menos activos y=20;%contantes de forsorilacion del CHEY (#molCheaA-1s-1)
Cy=700;%CONSTANTE DE FOSFORILACION DE Bp por CheA (#molCheB-1s-1)
gamma=1/2400;%dilucion constant (s-1)

ktdm=20; %constante de desligacion de compuesto intermedio de CheR-Receptor (s-1)
ktm=ktdm/60; %constante de ligacion de compuesto intermedio de CheR-Receptor (s-1) tal que la constante de saturacion Michaelis menten=60 receptores


alphaR = 130*gamma; %tasa de creacion CheR (#proteinas/s)
alphaB = 240*gamma; %tasa de creacion CheR (#proteinas/s)
alphaA = 6700*gamma; %tasa de creacion CheR (#proteinas/s)
alphaRe = 10000*gamma; %tasa de creacion CheR (#proteinas/s)

% c=3;
% K=36042.5;
% b=6619.134;



R = alphaR/gamma; %130
BT = alphaB/gamma; %230
AT = alphaA/gamma; %6700
TT = alphaRe/gamma; %1000


%% Condiciones iniciales
Tm = 0.067*TT;%importa para el nivel inicial
Bp = round(0.715*BT);
Ap = 0.36*AT;%0.056
RTdm = 0.79*R;



%inicialización variables
resultsAp=[];
resultsBp=[];
resultsT1=[];
resultsT2=[];
resultsT3=[];
resultsT4=[];
resultst=[];
resultsfood=[];
resultsRTdm=[];

n=0;
newAp=[];
t=0;
clust=0.0;

tarrayint=[];
aparrayint=[];
%% clustering
a=39000;
%y=10000;
y=abs(normrnd(39000,21000));
% m=4.5; %coeff de Hill
% K=72.15; %en uM
% B=1.02; %actvidad máxima
% a1=0.4654; %ajuste gaussiano
% b1=3.236e+04; %ajuste gaussiano
% c1=2.136e+04; %ajuste gaussiano
% 
% 
% if(a>=y)
%     
%     DifCl=0;
% else
%     nCl=4.5961e-06*(y-a)+0.0972; %nivel de clustering
%     CE=(1-nCl); %efecto del clustering
%     DifCl=(a1*exp(-((y-b1)/c1)^2))*CE;
% 
% end 
% A=B/1+(y/(K*600))^m; %actividad sin clustering,concentración en número de moléculas
% 
% CL=(A-DifCl)/A;%factor que disminuye la actividad según el nivel de clustering
CL=1.0;
%%
while (t<tmax)
    %definición de la comida
    
if(t<t1)
    L=a;
    
else
    %L=y(i);
    L=y;
end
%L= foodv(t);

%valor de reducción en actividad debido al clustering

if(t<t1)
    clust=1.0;
else
    clust=CL;
end
% Creacion y dilucion de proteinas


    dAT = (alphaA-gamma*AT)*dt;
    AT = AT + dAT;
    dBT = (alphaB-gamma*BT)*dt;
    BT = BT + dBT;
    dTT = (alphaRe-gamma*TT)*dt;
    TT = TT + dTT;
   %calculo de numero de receptores (promediado en el tiempo pues es muy rapido)
    T1=(TT-Tm)*(1/(1+kbn*L/kunb));
    T2=Tm*(1/(1+kbn*L/kunb));
    
    T3=(TT-Tm)*(L/(L+kunb/kbn));
    T4=Tm*(L/(L+kunb/kbn));
    % Fosforilacion promedio de CheA por los receptores metilados
    ACT =(k1*(T1+T4)+k2*(T2)+ k3*T3)*clust;%activated
    AVGT= T2;
    
    %formacion de compuesto transiente de CheR y Receptor
    
    dRTdm = (ktm*(TT-Tm)*(R-RTdm)-ktdm*RTdm)*dt;
    RTdm = RTdm + dRTdm;
    %Metilación del receptor por el compuesto trasiente
    dTm = (km*RTdm-kdm*Bp*AVGT)*dt;
    Tm = Tm+dTm;
    %Fosoforilacion de Bp
    dBp = (kb*Ap*(BT-Bp)- pb*Bp)*dt;
    Bp = Bp+dBp;
    
    %Fosforilacion de CheA
    dAp = (ACT*(AT-Ap)-(Cy+kb*(BT-Bp))*Ap)*dt;
    Ap = Ap+dAp;
    
    t = t+dt;
    n = n+1;
    %cont = cont+1;
    if(n==tmax/(1000*dt))
        resultsAp=[resultsAp ; Ap];
        resultst=[resultst ; t];
        resultsT1=[resultsT1 ; T1];
        resultsT2=[resultsT2 ;T2];
        resultsT3=[resultsT3; T3];
        resultsT4=[resultsT4; T4];
        resultsBp=[resultsBp ; Bp];
        resultsfood=[resultsfood;L];
        resultsRTdm=[resultsRTdm;RTdm];
        n=0;
    end
    if (t> t1 && t<= t1+4)
        tarrayint=[tarrayint,t];
        aparrayint=[aparrayint,Ap];
    end
end
%Normalized results
rnAp=resultsAp/resultsAp(10);% se va a normalizar con respecto a la simulación
% de las condiciones de referencia Frank & Vaknin (2013), fig 2

%área bajo la curva
%Otra normalización para poder calcular el área deseada (0)
r2Ap= rnAp-1;


%vector de 4s
aparrayint=aparrayint/mean(resultsAp(1:100))-1;


area1=abs(trapz(tarrayint,aparrayint));
areat=[areat;area1];

%%
% 
% hold on
% subplot(2,1,2)
% plot(resultst,resultsfood)
% %plot(resultst,resultsfood,'Color',CM(i,:))
% xlabel('Tiempo (s)','fontsize',11);
% ylabel('Número de moléculas de atractante','fontsize',7);
% xlim([0 tmax]);
% hold off
% 
% 
% hold on
% subplot(2,1,1)
% plot(resultst,rnAp,'LineWidth',2);
% %plot(resultst,rnAp,'Color',CM(i,:));
% xlabel('Time (s)','fontsize',11);
% ylabel('Number of molecules of Ap','fontsize',11);
% xlim([0 tmax]);
% title('Clustering included');
% hold off
% 
% %para ajustar parámetros
% 
% hold on
% %subplot(2,1,1)
% %plot(tExp,ApExp);
%  hold off

%%

%Registro de ligandos 
Lt=[Lt;L];
clt=[clt,clust];
end
%Legend=cell(2,1);
% for j=1:vc
%   
%  Legend{j}=num2str(clt(j)) ;
% end

%para el fitting

% Legend{1}="determinist" ;
% Legend{2}="experimental data" ;
% 
% 
% legend(Legend,'FontSize',7);
