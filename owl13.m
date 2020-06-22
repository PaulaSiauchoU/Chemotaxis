clear all
%code for stochastic simulation in E.coli. Single cell.
%clustering added
%Mutual Information calculation added
%Less time simulated to obtain more data resolution and less computation
%time. Data resolution was enhanced after the ligand jump.
%concentraciones iniciales variables, tamaño de salto proporcional al valor
%inicial
global kdm kunb km kbn kb pb k2 k1 k3 Cy gamma ktdm ktm alphaR alphaB alphaA alphaRe
tic
cells=100000;

areat=[];
Lt=[];
clt=[];
count3=0;
a=39000;
for i=1:cells;
    count3=count3+1;
    if count3==1000
        i
        count3=0;
    end
    %Tar receptor-
    
    cont=0;
    kdm = 0.00015;%para una adaptacion de 200 s mientras mas pequeno, mas se demora el sistema en adaptar.
    km = 0.05;
    kbn = 0.02;
    kunb = 300;%para que haya unas 20000 ligaciones y deligaciones por s.
    
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
    
    %parámetros clustering
    % n=3.0;
    % K=36042.5;
    % b=5261;
    
    R = alphaR/gamma; %130
    BT = alphaB/gamma; %230
    AT = alphaA/gamma; %4450
    TT = alphaRe/gamma; %4450
    
    
    
    %% Condiciones iniciales
    Tm = 0.0777*TT;%importa para el nivel inicial
    Bp = round(0.715*BT);
    Ap = 0.36*AT;%0.056
    RTdm = 0.79*R;
    
    y=0;
    t1=300;
    
    %Inicialización de variables
    L=[];
    t=0;
    resultsT1=[];
    resultsT2=[];
    resultsT3=[];
    resultsT4=[];
    resultsTm=[];
    resultsAt=[];
    resultsBt=[];
    resultsTt=[];
    resultsAp=[];
    resultsA=[];
    resultsBp=[];
    resultsB=[];
    resultsR=[];
    
    resultst=[];
    resultsfood=[];
    
    
    c=zeros(4,1);
    w=zeros(4,1);
    
    
    tend=400;
    tarrayint=[];
    aparrayint=[];
    
    %clustering
    y=abs(normrnd(a,5000));
    y2=abs(y+(y*0.1));
    m=4.5; %coeff de Hill
    K=72.15; %en uM
    B=1.02; %actvidad máxima
    a1=0.4654; %ajuste gaussiano
    b1=3.236e+04; %ajuste gaussiano
    c1=2.136e+04; %ajuste gaussiano
    
    if(y>=y2)
        
        DifCl=0;
    else
        nCl=4.5961e-06*(y-a)+0.0972; %nivel de clustering
        CE=(1-nCl); %efecto del clustering
        DifCl=(a1*exp(-((y-b1)/c1)^2))*CE;
        
    end
    A=B/1+(y/(K*600))^m; %actividad sin clustering,concentración en número de moléculas
    
    CL=(A-DifCl)/A;%factor que disminuye la actividad según el nivel de clustering
    
    %%
    while(t<tend)
        
        %food definition
        %valor de reducción en actividad debido al clustering
        
        %L=foodStep(t);
        if(t<t1)
            L=y;
            clust=1.0;
        else
            %L=y(i);
            L=y2;
            clust=CL;
        end
        
        
        
        
        %calculo de numero de receptores (promediado en el tiempo pues es muy rapido)
        T1=(TT-Tm)*(1/(1+kbn*L/kunb));
        T2=Tm*(1/(1+kbn*L/kunb));
        
        T3=(TT-Tm)*(L/(L+kunb/kbn));
        T4=Tm*(L/(L+kunb/kbn));
        %Fosforilacion promedio de CheA por los receptores metilados
        ACT =(k1*(T1+T4)+k2*(T2)+ k3*T3)*clust;
        AVGT=(T2);
        %Formaci�n compuesto trasciente
        %RTdm=ktm*(TT-Tm)*R/(ktm*(TT-Tm)+ktdm);
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
        if(cont==1000)
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
    %%
    %Normalized results
    rnAp=resultsAp/mean(resultsAp(1:100));
    %area under the curve
    
    %normalization (0)
    r2Ap= rnAp-1;
    
    
    %4s vector
    aparrayint=aparrayint/mean(resultsAp(1:100))-1;
    
    
    area1=abs(trapz(tarrayint,aparrayint));
    areat=[areat;area1];
    
    %Ligand record
    Lt=[Lt;y];
    clt=[clt,clust];
    
%     hold on
%     subplot(2,1,1)
%     plot(resultst,rnAp);
%     xlabel('Time','fontsize',15);
%     ylabel('Number of Ap','fontsize',15);
%     xlim([0, tend]);
%     legend('Ap');
%     hold off
%     hold on
%     subplot(2,1,2)
%     plot(resultst,resultsfood);
%     xlabel('Time','fontsize',15);
%     ylabel('Number of attractant molecules','fontsize',10);
%     xlim([0 tend]);
%     legend('L')
%     hold off
end
toc
dataout = [areat, Lt];
filename = ['Owl14_100k', '.dat'];
save(filename, 'dataout', '-ascii');