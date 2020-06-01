% This code calculates the discrepancy between SCF and PCF for
% single-plasmid communities, which corresponds to Fig. S2A in the
% supplementary information.
clear;
clc;
global mu eta kappa D lambda plasmid species capa;
plasmid=1;

lambdas=[-0.5 -0.2 -0.1 0 0.5 1];
[ha,pos]=tight_subplot(2,length(lambdas)/2,0.01,0.1,0.1);
for jj=1:length(lambdas)
    for rep=1:1000
        rep
        species=fix(1+50*rand);
        lambda=lambdas(jj)*(0.8+0.4*rand(species,plasmid));
        kappa=0.002*(0.8+0.4*rand(species,plasmid));
        D=0.02*(0.8+0.4*rand);
        
        capa=0;
        pre_groupnumber=max(1,fix(species/2));
        species_index=0;
        if pre_groupnumber==1
            groupnumber=1;
            species_index=ones(1,species);
            capa=ones(1,species);
        else
            species_tape=rand(1,species);
            species_insert=[-1 sort(rand(1,pre_groupnumber-1)) 10];
            pin=1;
            for i=1:pre_groupnumber+1
                logi=0;
                for j=1:species
                    if (species_tape(1,j)>species_insert(i))&&(species_tape(1,j)<=species_insert(i+1))
                    species_index(j)=pin;
                    logi=1;
                    end
                end
                if logi==1
                    pin=pin+1;
                end
            end
            groupnumber=max(species_index);
            pre_capa0=rand(1,groupnumber);
            pre_capa=pre_capa0/sum(pre_capa0);
            for i=1:species
                capa(i)=pre_capa(species_index(i));
            end
        end
    
        mu=0.3+0.5*rand(1,species);
        clear initial;
        initial(species*(1+plasmid))=0;
            for i=1:species
                initial(i)=0.01*capa(i);
                for j=1:plasmid
                    initial(species+plasmid*(i-1)+j)=0.01*capa(i)*0.1/plasmid;
                end
            end
        timespan=0:1000;
        eta=0.02*rand*rand(plasmid,species,species); 
        [t,y_new]=ode45(@multi_plasmid_new,timespan,initial);
        [t,y_old]=ode45(@multi_plasmid_old,timespan,initial);
        xx(j,rep)=sum(y_new(end,1+species:species*2))/sum(y_new(end,1:species));
        yy(j,rep)=sum(y_old(end,1+species:species*2))/sum(y_old(end,1:species));
    end
    axes(ha(jj));
    plot(xx,yy,'k.','markersize',15);hold on;
    if jj==1
    text(0.1,0.8,'$\bar{\lambda}=-0.5$','Interpreter','Latex','fontsize',30);
    end
    if jj==2
    text(0.1,0.8,'$\bar{\lambda}=-0.2$','Interpreter','Latex','fontsize',30);
    end
    if jj==3
    text(0.1,0.8,'$\bar{\lambda}=-0.1$','Interpreter','Latex','fontsize',30);
    end
    if jj==4
    text(0.1,0.8,'$\bar{\lambda}=0$','Interpreter','Latex','fontsize',30);
    end
    if jj==5
    text(0.1,0.8,'$\bar{\lambda}=0.5$','Interpreter','Latex','fontsize',30);
    end
    if jj==6
    text(0.1,0.8,'$\bar{\lambda}=1$','Interpreter','Latex','fontsize',30);
    end
    %set(gca,'fontsize',24);
    xticks([]);
    yticks([]);
    plot([0 1],[0 1],'--k');hold on;
    axis([0 1 0 1]);
end
axes(ha(length(lambdas)/2+1));
ylabel('SCF','fontsize',30,'fontweight','normal');
xlabel('PCF','fontsize',30,'fontweight','normal');
%title('plasmid abundance','fontsize',30,'fontweight','normal');

set(gcf,'position',[50 50 800 550]);


function dydt=multi_plasmid_new(t,y)
    global mu eta kappa D lambda  plasmid species capa
    dydt(species*(1+plasmid),1)=0;
    for i=1:species
        
        alphai=y(i)/(lambda(i)*y(i+species)+y(i));
        dydt(i,1)=alphai*mu(i)*y(i)*(capa(i)-y(i))-D*y(i);
        
        for j=1:plasmid
        betaij=1;
        muij=mu(i)/(1+lambda(i,j));
            summ=0;
            for k=1:species
                summ=summ+eta(j,k,i)*y(species+(k-1)*plasmid+j);
            end
        dydt(species+(i-1)*plasmid+j,1)=betaij*muij*y(species+(i-1)*plasmid+j)*(capa(i)-y(i))+(y(i)-y(species+(i-1)*plasmid+j))*summ-(kappa(i,j)+D)*y(species+(i-1)*plasmid+j);
        end
    end
    
    
end

function dydt=multi_plasmid_old(t,y)
    global mu eta kappa D lambda  plasmid species capa
    dydt(species*(1+plasmid),1)=0;
    for i=1:species
        
        
        alphai=(y(i)-y(i+species)+y(i+species)/(1+lambda(i)))/y(i);
        dydt(i,1)=alphai*mu(i)*y(i)*(capa(i)-y(i))-D*y(i);
        
        for j=1:plasmid
        betaij=1;
        muij=mu(i)/(1+lambda(i,j));
            summ=0;
            for k=1:species
                summ=summ+eta(j,k,i)*y(species+(k-1)*plasmid+j);
            end
        dydt(species+(i-1)*plasmid+j,1)=betaij*muij*y(species+(i-1)*plasmid+j)*(capa(i)-y(i))+(y(i)-y(species+(i-1)*plasmid+j))*summ-(kappa(i,j)+D)*y(species+(i-1)*plasmid+j);
        end
    end
    
    
end