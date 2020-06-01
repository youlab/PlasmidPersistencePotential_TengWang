% This code simulates the steady-state distribution of plasmids in a single community.
% It can be adapted to generate Fig.1C and Fig.S1B in the manuscript of "the persistence potential of plamids".
clear all;
clc;
global mu eta kappa D lambda plasmid species capa species_index;
% Here, 'mu', 'eta', 'kappa', 'D' and 'lambda' are the kinetic parameters defined in
% the main text. 'plasmid' is the total number of the plasmids in the
% community. 'species' is the total number of 'species'. 'capa' is the
% vector of the maximum carrying capacities. 'species_index' is the vector
% of which niche the species belongs to.

timespan=[1:1:35000];
% 'timespan' defines the time range that the simulation runs.

species=20;
plasmid=20;

% In the following paragraph, we divide the community into multiple niches
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
    groupnumber=max(species_index);% Here, 'groupnumber' represents the total number of niches;
    
    % In the following paragraph, we randomize the maximum carrying
    % capacity of each niche.
    pre_capa0=rand(1,groupnumber);
    pre_capa=pre_capa0/sum(pre_capa0);
    for i=1:species
        capa(i)=pre_capa(species_index(i));
    end
end
    
% In the followinf praragraph, we define the initial state of the
% community.
    initial(species*(1+plasmid))=0;
    for i=1:species
        pin2=0;
        for j=1:species
            if species_index(j)==species_index(i)
                pin2=pin2+1;
            end
        end       
        initial(i)=0.1*capa(i)/pin2;
        for j=1:plasmid
            initial(species+plasmid*(i-1)+j)=0.1*rand*capa(i)/pin2*rand;
        end
    end
    qt=0;
    
    
    % Here, we randomize the kinetic parameters.
    mu=1*rand(1,species);
    D=0.01*rand;
    kappa=0.05*rand(species,plasmid);
    lambdaa=0.2;
    for i=1:plasmid
        lambda(:,i)=lambdaa*(1-2*rand)*rand(species,1);
    end
    eta=0*ones(plasmid,species,species);
    for jk=1:plasmid
        etamax=0.1*rand;
        connec=fix(rand*species^2);
        R1=rand(species,species);
        R2=R1(:);
        R3=sort(R2);
        if connec<1
            crit=R3(end)+1;
        else
            crit=R3(length(R3)-(connec-1));
        end
        for i=1:species
            for j=1:species
                if(R1(i,j)>=crit)
                    eta(jk,i,j)=etamax*rand;
                end
            end
        end
    end
    
[t,y]=ode45(@multi_plasmid,timespan,initial);
    
% Here, we calculate the relative abundance of each species 
for i=1:species
    abun_rela(i,:)=y(end,species+(species-1)*plasmid+1:species+(species-1)*plasmid+plasmid)./y(end,i);
end

% Here, we calculate the relative abundance of each plasmid
for i=1:plasmid
y_p(:,i)=sum(y(:,species+i:plasmid:species+plasmid*(species-1)+i),2);
end

figure (1);% The distribution of plasmids across different species at the end of simulation
[X1,Y1]=meshgrid([1:species+1],[1:plasmid+1]);
for i=1:species
    for j=1:plasmid
        Z1(j,i)=y(end,species+plasmid*(i-1)+j)./y(end,i);
    end
end
Z1(:,species+1)=Z1(:,species);
Z1(plasmid+1,:)=Z1(plasmid,:);
surf(X1,Y1,Z1);shading flat;
set(gca,'fontsize',18);
xlabel('species','fontsize',24);
ylabel('plasmid','fontsize',24);view([0 90]);
axis([1 species+1 1 plasmid+1]);
c = jet;
colormap(c);
colorbar northoutside;
set(gcf,'position',[100 100 400 400]);

figure (2);% The relative abundance of each species at the end of the simulation
bar(1:species,y(end,1:species)/sum(y(end,1:species),2),1,'stack','edgecolor','none','facecolor','k');
axis([0.5 species+0.5 0 1.2*max(y(end,1:species))]);
set(gcf,'position',[100 100 500 150]);

figure (3);% The relative abundance of each plasmid at the end of the simulation
bar(1:plasmid,y_p(end,plasmid:-1:1)/sum(y(end,1:species),2),1,'stack','edgecolor','none','facecolor','k');
axis([0.5 plasmid+0.5 0 1.2*max(y_p(end,1:plasmid))]);
set(gcf,'position',[100 100 500 150]);

figure (4);% The distribution of plasmids across different species at the initial state
[X2,Y2]=meshgrid([1:species+1],[1:plasmid+1]);
for i=1:species
    for j=1:plasmid
        Z2(j,i)=y(1,species+plasmid*(i-1)+j)./y(1,i);
    end
end
Z2(:,species+1)=Z2(:,species);
Z2(plasmid+1,:)=Z2(plasmid,:);
surf(X2,Y2,Z2);shading flat;
set(gca,'fontsize',18);
xlabel('species','fontsize',24);
ylabel('plasmid','fontsize',24);view([0 90]);
axis([1 species+1 1 plasmid+1]);
c = jet;
colormap(c);
colorbar northoutside;
set(gcf,'position',[100 100 400 400]);

figure (5);% The relative abundance of each species at the initial state
bar(1:species,y(1,1:species)/sum(y(1,1:species),2),1,'stack','edgecolor','none','facecolor','k');
axis([0.5 species+0.5 0 1.2*max(y(end,1:species))]);
set(gcf,'position',[100 100 500 150]);

figure (6);% The relative abundance of each plasmid at the initial state
bar(1:plasmid,y_p(1,plasmid:-1:1)/sum(y(1,1:species),2),1,'stack','edgecolor','none','facecolor','k');
axis([0.5 plasmid+0.5 0 1.2*max(y_p(end,1:plasmid))]);
set(gcf,'position',[100 100 500 150]);

function dydt=multi_plasmid(t,y)
        global mu eta kappa D lambda  plasmid species capa species_index
        dydt(species*(1+plasmid),1)=0;
        for i=1:species
            sum=0;
            for j=1:plasmid
            sum=sum+lambda(i,j)*y(species+plasmid*(i-1)+j);
            end

            alphai=y(i)/(sum+y(i));
            pt=0;
            for j=1:species
                if species_index(j)==species_index(i)
                    pt=pt+y(j);
                end
            end
            dydt(i,1)=alphai*mu(i)*y(i)*(capa(i)-pt)-D*y(i);

            for j=1:plasmid
            betaij=(1+lambda(i,j))/(1+lambda(i,j)+(sum-lambda(i,j)*y(species+plasmid*(i-1)+j))/y(i));
            muij=mu(i)/(1+lambda(i,j));
                summ=0;
                for k=1:species
                    summ=summ+eta(j,k,i)*y(species+(k-1)*plasmid+j);
                end
            dydt(species+(i-1)*plasmid+j,1)=betaij*muij*y(species+(i-1)*plasmid+j)*(capa(i)-pt)+(y(i)-y(species+(i-1)*plasmid+j))*summ-(kappa(i,j)+D)*y(species+(i-1)*plasmid+j);
            end
        end


end