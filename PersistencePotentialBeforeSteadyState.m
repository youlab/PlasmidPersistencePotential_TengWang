clear;
clc;
global mu eta kappa D lambda plasmid species capa species_index;
timespan=0:2000;

repp=1;

for rep=1:50000
    rep
    speciess(rep)=1+fix(10*rand);
    species=speciess(rep);
    plasmids(rep)=1+fix(10*rand);
    plasmid=plasmids(rep);
    
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
        groupnumber=max(species_index);
        pre_capa0=rand(1,groupnumber);
        pre_capa=pre_capa0./sum(pre_capa0);
        for i=1:species
            capa(i)=pre_capa(species_index(i));
        end
    end
    
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
    
    mu=0.1+0.7*rand(1,species);
    
    DD(rep)=0.01+0.04*rand;
    D=DD(rep);
    
    kappaa(rep)=0.1*rand;
    kappa=kappaa(rep)*(0.9+0.2*rand(species,plasmid));
    
    lambdaa(rep)=0.02*rand;
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
        
        etamax=0.8*rand;
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
    
    pttt=sum(y(end,1:1:species));
    xx=y(end,1:1:species)./pttt;
    diver=exp(-sum(xx.*log(xx)));
    
    pttt100=sum(y(100,1:1:species));
    xx100=y(100,1:1:species)./pttt100;
    diver100=exp(-sum(xx100.*log(xx100)));
    
    pttt200=sum(y(200,1:1:species));
    xx200=y(200,1:1:species)./pttt200;
    diver200=exp(-sum(xx200.*log(xx200)));

    pttt500=sum(y(500,1:1:species));
    xx500=y(500,1:1:species)./pttt500;
    diver500=exp(-sum(xx500.*log(xx500)));

    for iii=1:plasmid
        qt(repp)=sum(y(end,species+iii:plasmid:species+plasmid*(species-1)+iii));
        pt(repp)=sum(y(end,1:1:species));
        fract(repp)=qt(repp)/pt(repp);

        mu_ave(repp)=0;
        for i=1:species
            mu_ave(repp)=mu_ave(repp)+mu(i)*y(end,i)/pt(repp);
        end

        kappa_ave(repp)=0;
        for i=1:species
            kappa_ave(repp)=kappa_ave(repp)+kappa(i,iii)*y(end,i)/pt(repp);
        end

        lambda_ave(repp)=0;
        for i=1:species
            lambda_ave(repp)=lambda_ave(repp)+lambda(i,iii)*y(end,i)/pt(repp);
        end

        etasum(repp)=0;
        for i=1:species
            for j=1:species
            etasum(repp)=etasum(repp)+eta(iii,i,j)*y(end,i)*y(end,j)/(pt(repp))^2;
            end
        end
        omega(repp)=etasum(repp)/(mu_ave(repp)/(mu_ave(repp)-diver*D)*(D+kappa_ave(repp)-D/(1+lambda_ave(repp))));
        
        qt100(repp)=sum(y(100,species+iii:plasmid:species+plasmid*(species-1)+iii));
        pt100(repp)=sum(y(100,1:1:species));
        fract100(repp)=qt100(repp)/pt100(repp);

        mu_ave100(repp)=0;
        for i=1:species
            mu_ave100(repp)=mu_ave100(repp)+mu(i)*y(100,i)/pt100(repp);
        end

        kappa_ave100(repp)=0;
        for i=1:species
            kappa_ave100(repp)=kappa_ave100(repp)+kappa(i,iii)*y(100,i)/pt100(repp);
        end

        lambda_ave100(repp)=0;
        for i=1:species
            lambda_ave100(repp)=lambda_ave100(repp)+lambda(i,iii)*y(100,i)/pt100(repp);
        end

        etasum100(repp)=0;
        for i=1:species
            for j=1:species
            etasum100(repp)=etasum100(repp)+eta(iii,i,j)*y(100,i)*y(100,j)/(pt100(repp))^2;
            end
        end
        omega100(repp)=etasum100(repp)/(mu_ave100(repp)/(mu_ave100(repp)-diver100*D)*(D+kappa_ave100(repp)-D/(1+lambda_ave100(repp))));

        qt200(repp)=sum(y(200,species+iii:plasmid:species+plasmid*(species-1)+iii));
        pt200(repp)=sum(y(200,1:1:species));
        fract200(repp)=qt200(repp)/pt200(repp);

        mu_ave200(repp)=0;
        for i=1:species
            mu_ave200(repp)=mu_ave200(repp)+mu(i)*y(200,i)/pt200(repp);
        end

        kappa_ave200(repp)=0;
        for i=1:species
            kappa_ave200(repp)=kappa_ave200(repp)+kappa(i,iii)*y(200,i)/pt200(repp);
        end

        lambda_ave200(repp)=0;
        for i=1:species
            lambda_ave200(repp)=lambda_ave200(repp)+lambda(i,iii)*y(200,i)/pt200(repp);
        end

        etasum200(repp)=0;
        for i=1:species
            for j=1:species
            etasum200(repp)=etasum200(repp)+eta(iii,i,j)*y(200,i)*y(200,j)/(pt200(repp))^2;
            end
        end
        omega200(repp)=etasum200(repp)/(mu_ave200(repp)/(mu_ave200(repp)-diver200*D)*(D+kappa_ave200(repp)-D/(1+lambda_ave200(repp))));

        qt500(repp)=sum(y(500,species+iii:plasmid:species+plasmid*(species-1)+iii));
        pt500(repp)=sum(y(500,1:1:species));
        fract500(repp)=qt500(repp)/pt500(repp);

        mu_ave500(repp)=0;
        for i=1:species
            mu_ave500(repp)=mu_ave500(repp)+mu(i)*y(500,i)/pt500(repp);
        end

        kappa_ave500(repp)=0;
        for i=1:species
            kappa_ave500(repp)=kappa_ave500(repp)+kappa(i,iii)*y(500,i)/pt500(repp);
        end

        lambda_ave500(repp)=0;
        for i=1:species
            lambda_ave500(repp)=lambda_ave500(repp)+lambda(i,iii)*y(500,i)/pt500(repp);
        end

        etasum500(repp)=0;
        for i=1:species
            for j=1:species
            etasum500(repp)=etasum500(repp)+eta(iii,i,j)*y(500,i)*y(500,j)/(pt500(repp))^2;
            end
        end
        omega500(repp)=etasum500(repp)/(mu_ave500(repp)/(mu_ave500(repp)-diver500*D)*(D+kappa_ave500(repp)-D/(1+lambda_ave500(repp))));

        repp=repp+1;
    end
end

figure(1);
subplot(1,4,1);
plot(omega100,fract100,'r.','markersize',5);hold on;
axis([0 5 0 1]);
subplot(1,4,2);
plot(omega200,fract200,'r.','markersize',5);hold on;
axis([0 5 0 1]);
subplot(1,4,3);
plot(omega500,fract500,'r.','markersize',5);hold on;
axis([0 5 0 1]);
subplot(1,4,4);
plot(omega,fract,'r.','markersize',5);hold on;
axis([0 5 0 1]);
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