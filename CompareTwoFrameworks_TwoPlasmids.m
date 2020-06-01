% This code calculates the discrepancy between SCF and PCF for
% two-plasmid communities, which corresponds to Fig. S2B in the
% supplementary information.
clear;
clc;
global mu eta kappa D lambda plasmid species lambdack;
species=1;
plasmid=2;
kappa=0.001*ones(species,plasmid);
D=0.005;
mu=0.3*ones(1,species);
timespan=0:1000;
eta=[0.01;0.01];

interv=0.05;
lambdas=-0.4:interv:1;
lambdacks=-0.4:interv:1;

for i=1:length(lambdas)
    i
    lambda=[lambdas(i) lambdas(i)];
    for j=1:length(lambdacks)
        j
        lambdack=lambdacks(j);
        initial1=[0.04 0.02 0.02];
        [t,y1]=ode45(@two_plasmid_new,timespan,initial1);
        frac1=y1(end,2)/y1(end,1);
        initial2=[0.01 0.01 0.01 0.01];
        [t,y2]=ode45(@two_plasmid_old,timespan,initial2);
        frac2=(y2(end,2)+y2(end,4))/(y2(end,1)+y2(end,2)+y2(end,3)+y2(end,4));
        z(j,i)=frac1/frac2;
    end
end

lambdass=[lambdas lambdas(end)+interv];
lambdackss=[lambdacks lambdacks(end)+interv];
z(:,length(lambdass))=z(:,length(lambdas));
z(length(lambdackss),:)=z(length(lambdacks),:);

surf(lambdass,lambdackss,z);hold on;
view([0 90]);shading flat;
set(gca,'fontsize',24);
xlabel('\lambda_1 and \lambda_2','fontsize',24);
ylabel('\lambda_{1,2}','fontsize',24);
%title('discrepancy','fontsize',30,'fontweight','normal');
axis([min(lambdass) max(lambdass) min(lambdackss) max(lambdackss)]);
colormap jet;
colorbar;
set(gcf,'position',[100 100 500 450])
contour3(lambdass,lambdackss,z,[1 1],'k--','linewidth',5);

function dydt=two_plasmid_new(t,y)
    global mu eta kappa D lambda  plasmid species
    dydt(species*(1+plasmid),1)=0;
    pt=0;
    for i=1:species
        pt=pt+y(i);
    end
    for i=1:species
        sum=0;
        for j=1:plasmid
        sum=sum+lambda(i,j)*y(species+plasmid*(i-1)+j);
        end
        
        alphai=y(i)/(sum+y(i));
        dydt(i,1)=alphai*mu(i)*y(i)*(1-pt)-D*y(i);
        
        for j=1:plasmid
        betaij=(1+lambda(i,j))/(1+lambda(i,j)+(sum-lambda(i,j)*y(species+plasmid*(i-1)+j))/y(i));
        muij=mu(i)/(1+lambda(i,j));
            summ=0;
            for k=1:species
                summ=summ+eta(j,k,i)*y(species+(k-1)*plasmid+j);
            end
        dydt(species+(i-1)*plasmid+j,1)=betaij*muij*y(species+(i-1)*plasmid+j)*(1-pt)+(y(i)-y(species+(i-1)*plasmid+j))*summ-(kappa(i,j)+D)*y(species+(i-1)*plasmid+j);
        end
    end
    
    
end

function dydt=two_plasmid_old(t,y)
global lambda lambdack mu eta  kappa D;
p=y(1)+y(2)+y(3)+y(4);
dydt=[mu*y(1)*(1-p)-eta(1)*y(1)*y(2)+eta(2)*y(1)*y(3)-eta(1)*y(1)*y(4)-eta(2)*y(1)*y(4)+kappa(1)*y(2)+kappa(2)*y(3)-D*y(1);
    mu/(1+lambda(1))*y(2)*(1-p)+eta(1)*y(1)*y(4)+eta(1)*y(2)*y(1)-eta(2)*y(2)*y(4)-eta(2)*y(3)*y(2)+kappa(2)*y(4)-kappa(1)*y(2)-D*y(2);
    mu/(1+lambda(2))*y(3)*(1-p)+eta(2)*y(1)*y(4)+eta(2)*y(3)*y(1)-eta(1)*y(3)*y(4)-eta(1)*y(2)*y(3)+kappa(1)*y(4)-kappa(2)*y(3)-D*y(3);
    mu/(1+lambdack)*y(4)*(1-p)+eta(1)*y(2)*y(3)+eta(2)*y(2)*y(3)+eta(1)*y(2)*y(4)+eta(2)*y(3)*y(4)-kappa(1)*y(4)-kappa(2)*y(4)-D*y(4);];
end