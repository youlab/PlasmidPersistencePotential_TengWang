% This code generates the plasmid-abundance as a function of persistence
% potential, which corresponds to Fig. 3C in the manuscript of "the persistence potential of plamids".
clear;
clc;
global mu eta kappa D lambda plasmid species capa species_index;
% Here, 'mu', 'eta', 'kappa', 'D' and 'lambda' are the kinetic parameters defined in
% the main text. 'plasmid' is the total number of the plasmids in the
% community. 'species' is the total number of 'species'. 'capa' is the
% vector of the maximum carrying capacities. 'species_index' is the vector
% of which niche the species belongs to.

timespan=0:30000;
% 'timespan' defines the time range that the simulation runs.

repp=1;

for rep=1:20 % Here, we define how many communities to simulate
    rep
    speciess(rep)=5+fix(95*rand);
    species=speciess(rep);% species number randomized from 5 to 100
    plasmids(rep)=1+fix(49*rand);
    plasmid=plasmids(rep);% plasmid number randomized from 1 to 50
    
    % In the following paragraph, we divide the community into multiple niches
    pre_groupnumber=max(1,fix(species*rand));
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
    clear initial t y;
    initial(species*(1+plasmid))=0;
    for i=1:species
        pin2=0;
        for j=1:species
            if species_index(j)==species_index(i)
                pin2=pin2+1;
            end
        end       
        initial(i)=0.01*capa(i)/pin2;
        for j=1:plasmid
            initial(species+plasmid*(i-1)+j)=0.01*capa(i)/pin2*0.1/plasmid;
        end
    end
    qt=0;
    
    % Here, we randomize the kinetic parameters.
    mu=0.4+0.4*rand(1,species);
    DD(rep)=0.001+0.004*rand;
    D=DD(rep);
    kappaa(rep)=0.002*rand;
    kappa=kappaa(rep)*(0.9+0.2*rand(species,plasmid));
    lambdaa(rep)=0.2*rand;
    lambda=lambdaa(rep)*(0.9+0.2*rand(species,plasmid));
    eta=0*ones(plasmid,species,species);
    for ety=1:plasmid
        connec=fix(rand*species^2);
        connecc(rep)=connec;
        R1=rand(species,species);
        R2=R1(:);
        R3=sort(R2);
        if connec<1
            crit=R3(end)+1;
        else
            crit=R3(length(R3)-(connec-1));
        end
        
        etamax=0.02*rand;
        ADJ=0*ones(species,species);
        for i=1:species
            for j=1:species
                if(R1(i,j)>=crit)
                    eta(ety,i,j)=etamax*rand;
                    ADJ(i,j)=1;
                end
            end
        end
    end
    
    [t,y]=ode45(@multi_plasmid,timespan,initial);

    for iii=1:plasmid
        qt(repp)=sum(y(end,species+iii:plasmid:species+plasmid*(species-1)+iii));% calculating the total abundance of each plasmid
        pt(repp)=sum(y(end,1:1:species));% calculating the total population size
        fract(repp)=qt(repp)/pt(repp);% calculating the relative abundance of each plasmid

        mu_ave(repp)=0;%calculating average value of mu
        for i=1:species
            mu_ave(repp)=mu_ave(repp)+mu(i)*y(end,i)/pt(repp);
        end

        kappa_ave(repp)=0;%calculating average value of kappa
        for i=1:species
            kappa_ave(repp)=kappa_ave(repp)+kappa(i,iii)*y(end,i)/pt(repp);
        end

        lambda_ave(repp)=0;%calculating average value of lambda
        for i=1:species
            lambda_ave(repp)=lambda_ave(repp)+lambda(i,iii)*y(end,i)/pt(repp);
        end

        etasum(repp)=0;%calculating average value of eta
        for i=1:species
            for j=1:species
            etasum(repp)=etasum(repp)+eta(iii,i,j)*y(end,i)*y(end,j)/(pt(repp))^2;
            end
        end
        
        % Next, we calculate the diversity of the community, defined by the
        % Shannon effective number
        diversity=exp(sum(-y(end,1:1:species)/sum(y(end,1:1:species)).*log(y(end,1:1:species)/sum(y(end,1:1:species)))));
        
        % Now we can calculate the omega value
        omega(repp)=etasum(repp)/(mu_ave(repp)/(mu_ave(repp)-diversity*D)*(D+kappa_ave(repp)-D/(1+lambda_ave(repp))));
        repp=repp+1;
    end
end

for i=1:length(omega)
    plot(omega(i),fract(i),'r.','markersize',10);hold on;
end
axis([0 5 0 1]);
xlabel('\omega','fontsize',30);
ylabel('q_T/p_T','fontsize',30);
set(gca,'fontsize',24);
set(gcf,'position',[200 200 400 400]);



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