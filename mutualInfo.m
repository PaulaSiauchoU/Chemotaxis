function I=mutualInfo(areat,Lt)
clear I 
clear x
clear y

%% Parametrer of the distrubution of the data
cr=corr(Lt,areat);
mnA=mean(areat);
mnL=mean(Lt);
nL=(Lt-mnL)/mnL;
nA=(areat-mnA)/mnA;
mnnL=mean(nL);
mnnA=mean(nA);
%%
hn=[];
hd=[];

x=nL;
y=nA;
%for k=1:5
    tic
%     %% If the data is a sample of a desired distribution
%      n=10^k;
%      r=cr; %correlation coefficient
%      prom= [mnnL,mnnA]; %mean [x,y]
%      cov= [1,r;r,1]; %cov matrix
%      x=[];
%      y=[];
%      for i=1:n %1000 points with that distribution
%          z=mvnrnd(prom,cov); %to generate a bivariate gaussian distribution
%          %f.write(str(z[0])+','+str(z[1])+'\n')#escribirla en un documento
%          x=[x;z(1)];
%          y=[y;z(2)];
%      end
%%
    nx=zeros(length(x),1);
    ny=zeros(length(y),1);
    
    mx=max(x)-min(x);
    my=max(y)-min(y);
    
    mm=max(mx,my);
    
    
    for i=1:length(x)
        dx=zeros(length(x),1);
        dy=zeros(length(y),1);
        d= zeros(length(y),1);
        for j=1:length(y)
            dx(j)=abs(x(i)-x(j));
            dy(j)=abs(y(i)-y(j));
            d(j)=max(dx(j),dy(j));
        end
        dd=sort(d);
        mini=dd(2);
        for j=1:length(x)
            if dx(j)<mini
                nx(i)=nx(i)+1;
            end
            if dy(j)<mini
                ny(i)=ny(i)+1;
            end
        end
        
        
    end
    phix=0;
    phiy=0;
    for i =1:length(nx)
        phix=phix+psi(nx(i)+1)/length(nx);
        phiy=phiy+psi(ny(i)+1)/length(ny);
    end
     I=psi(1)+psi(length(x))- phix - phiy;
%     teo=-(0.5*log(1-r*r));
%     diff=abs(I-abs(teo));
%     hn=[hn,n];
%     hd=[hd,diff]; 
    toc 
%end

end
