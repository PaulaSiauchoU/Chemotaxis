clear all
%Código para simulación estocástica de Quimiotaxis en E.coli, poblacional.
%Cálculo de Información Mutua añadido.
%sin clustering
%Se simula menos tiempo con respecto al determinista para tener mayor
%resolución de los datos pero no más tiempo computacional.
%La resolución de los datos de aumenta despues del salto del ligando.
%algunos cáculos para comprobar el funcionamiento de los números aleatorios
%añadidos.
%concentraciones iniciales variables, mismo tamaño de salto

global kdm kunb km kbn kb pb k2 k1 k3 Cy gamma ktdm ktm alphaR alphaB alphaA alphaRe
tic
cells=20;
moda=[];
% frec=[];
% modaA=[];
% frecA=[];
% tot=[];
areat=[];
% M1=0;
% F1=0;
% M2=0;
% F2=0;
Lt=[];
clt=[];
% nums=[];
count3=0;
%for x=1:3;
a=39000;
for i=1:cells;
    
    count3=count3+1;
    if count3==1000 %Visualización del progreso de la simulación.
        i
        count3=0;
    end
    
    
    cont=0;
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
    
    %número promedio de proteínas en célula
    R = alphaR/gamma; %130
    BT = alphaB/gamma; %230
    AT = alphaA/gamma; %6700
    TT = alphaRe/gamma; %1000
    
    
    %% Condiciones iniciales
    Tm = 0.067*TT;%importa para el nivel inicial
    Bp = round(0.715*BT);
    Ap = 0.36*AT;%0.056
    RTdm = 0.79*R;
    
    y=0;
    t1=300;
   
    
    L=[];
    t=0;
    resultsT1=[];
    resultsT2=[];
    resultsT3=[];
    resultsT4=[];
    resultsAt=[];
    resultsBt=[];
    resultsTt=[];
    resultsAp=[];
    resultsA=[];
    resultsBp=[];
    resultsRTdm=[];
    
    resultst=[];
    resultsfood=[];
    
    
    c=zeros(4,1);
    w=zeros(4,1);
    
    
    tend=500;
    y=abs(normrnd(a,5000));
    y2=abs(y+10000);
    tarrayint=[];
    aparrayint=[];
    
     %clustering
    clust=1.0;
    
    while(t<tend)
        
        %Definición de la entrada.
        
        %L=foodStep(t);
        if(t<t1)
            L=y;
            
        else
            %L=y(i);
            L=y2;
        end
        %calculo de numero de receptores (promediado en el tiempo pues es muy rapido)
        T1=(TT-Tm)*(1/(1+kbn*L/kunb));
        T2=Tm*(1/(1+kbn*L/kunb));
        
        T3=(TT-Tm)*(L/(L+kunb/kbn));
        T4=Tm*(L/(L+kunb/kbn));
        %Fosforilacion promedio de CheA por los receptores metilados
        ACT =(k1*(T1+T4)+k2*(T2)+ k3*T3)*clust;
        AVGT=(T2);
        %Formación compuesto trasciente
        
        RTdm=ktm*(TT-Tm)*R/(ktm*(TT-Tm)+ktdm);
        %Fosforilacion de CheA
        
        Ap=AT*ACT/(ACT+kb*(BT-Bp)+Cy);
        
        
        w(1) = km*RTdm;
        w(2) = kdm*Bp*AVGT;
        w(3) = kb*Ap*(BT-Bp);
        w(4) = pb*Bp;
        
        
        
        
        c(1)=w(1);
        for j=2:length(c);
            c(j)=c(j-1)+w(j);
        end
        
        ct=c(length(c));
        z1=rand();
        z2=rand();
        tau=(-log(z1))/ct;
        uct=z2*ct;
        
        if (uct<=c(1))  % w1
            Tm = Tm + 1;
            %R= R+1;
            %RTdm=RTdm - 1;
            
        elseif ((uct>c(1)) && (uct<=c(2)))  % w2
            Tm=Tm-1;
            
        elseif ((uct>c(2)) && (uct<=c(3)))  % w3
            Bp=Bp+1;
            %Ap=Ap-1;
        elseif ((uct>c(3)) && (uct<=c(4)))  % w4
            Bp=Bp-1;
            
        end
        t=t+tau;
        cont = cont+1;
        if(cont==1000)%guardar cada n números
            resultst=[resultst ; t ];
            resultsAp=[resultsAp;Ap];
            resultsfood=[resultsfood ; L];
            cont=0;
        end
        
        if (t> t1 && t<= t1+4)
            tarrayint=[tarrayint,t];
            aparrayint=[aparrayint,Ap];
        end
    end
    %% Procesamiento de los datos
    %Resultados normalizados a los niveles de adaptación perfecta
    rnAp=resultsAp/mean(resultsAp(1:100));
    %área bajo la curva
    %Otra normalización para poder calcular el área deseada (0)
    r2Ap= rnAp-1;
    
    
    %vector de 4s
    aparrayint=aparrayint/mean(resultsAp(1:100))-1;
    
    
    area1=abs(trapz(tarrayint,aparrayint));
    areat=[areat;area1];
    
    %Ligand record
    Lt=[Lt;L-10000];
    clt=[clt,clust];
    %random number record
    %         nums=[nums;y];
    hold on
    subplot(2,1,1)
    plot(resultst,rnAp);
    xlabel('Time','fontsize',15);
    ylabel('Number of Ap','fontsize',15);
    xlim([0, tend]);
    legend('Ap');
    hold off
    hold on
    subplot(2,1,2)
    plot(resultst,resultsfood);
    xlabel('Time','fontsize',15);
    ylabel('Number of attractant molecules','fontsize',10);
    xlim([0 tend]);
    legend('L')
    hold off
end
% [M1,F1]=mode(nums);
% [M2,F2]=mode(areat);
% moda=[moda,M1];
% frec=[frec,F1];
% modaA=[modaA,M2];
% frecA=[frecA,F2];
% snums=sort(nums);
% tot=[tot,snums];
toc
% dataout = [areat, Lt];
% filename = ['Owl11_100k', '.dat'];
% save(filename, 'dataout', '-ascii');
%end
