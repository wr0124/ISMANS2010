function autostar(n,N,L,agentwidth,alfa,timestep,agenttime,Thermo)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%the paramaters of this model%%%%%%%%%%%%%%%%%
% n=3000;     %the number of agent number
% N=300;       %the number of star number
% L=N;    %the range of the star's position
% agentwidth=50;
% alfa=2;
% timestep=8000;           %the total timesteps
% agenttime=2;    %the number of agents' evoltion
% Thermo=L*0;


percenttogive=0.05; %the percent of links to give when a new stars splitting
agentadjacency=zeros(n,n);   %connection of the agent
starage=ones(1,N);
speakelo=randi(10,1,N);   %the eloauence of speak ranged (1,10)
agent2star=zeros(N,n);
pc=0.7;    %probability for inital agent connecting probability
timeintermerg=0;
timeinterspli=0;
rewiretime=0;
breaktime=0;
addtime=0;
linksplit=0;
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%output file%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%the information of star%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd 'C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\'
dos(['mkdir n' num2str(n) 'N' num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' num2str(agenttime)  'Thermo' num2str(Thermo) ])

fid1=fopen(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid1numberofstar''.txt'],'wt');


fid2=fopen(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid2positionofstar' '.txt'],'wt');


fid3=fopen(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid3distanceofstarmerg' '.txt'],'wt');


fid4=fopen(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid4intervalofstarmerg'  '.txt'],'wt');


fid5=fopen(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid5distanceofstarsplit' '.txt'],'wt');


fid6=fopen(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid6intervalofstarsplit' '.txt'],'wt');



fid7=fopen(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid7aveopinion' '.txt'],'wt');


%%%%%%%%%%%%%%%%%the information of agent
fid8=fopen(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid8avedegree' '.txt'],'wt');

fid9=fopen(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid9degreedistribution''.txt'],'wt');


%%%%%%%%%%%%%%%%%the information of their relation
fid10=fopen(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid10numberofagentinstar' '.txt'],'wt');

fid11=fopen(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid11timetaken' '.txt'],'wt');

fid12=fopen(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid12agentaveopinion'  '.txt'],'wt');

fid13=fopen(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid13staraveopinion' '.txt'],'wt');

fid14=fopen(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid14rewiringtime' '.txt'],'wt');

%%%%%%%%%%%%%%%%%%%the evolution information of the agent
fid15=fopen(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid15rewiringdis' '.txt'],'wt');
fid16=fopen(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid16breakingtime' '.txt'],'wt');
fid17=fopen(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid17breakingdis' '.txt'],'wt');
fid18=fopen(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid18addingtime' '.txt'],'wt');
fid19=fopen(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid19addingdis' '.txt'],'wt');

%%%%%%%%%%%%%%%%%%%the evolution information of the agents' inner connect  with star
fid20=fopen(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\avgdegreetostar' '.txt'],'wt');
fid21=fopen(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\edgeconnectstrength' '.txt'],'wt');

fid22=fopen(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\eloposition' '.txt'],'wt');

fid23=fopen(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\agedistribution' '.txt'],'wt');

fid24=fopen(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\link2split' '.txt'],'wt');

fid25=fopen(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\nodeconnectstrength' '.txt'],'wt');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(fid11,'the program begins at:   %.5g / %.5g  / %.5g/   %.5g:   %.5g:   %.5g.   \n',clock );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%inital configuation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%1. agent inital connected with each other randomly with Poisson distribution
for i=1:n-1
    agentadjacency(i,1+i)=1;
    agentadjacency(1+i,i)=1;
end
leftlink=round(n*(n-1)*pc/2-(n-1));
co=0;
while(1)
    ina=randi(n);
    inb=randi(n);
    if (agentadjacency(ina, inb)~=1) && (ina~=inb)
        agentadjacency(ina,inb)=1;
        agentadjacency(inb,ina)=1;
        co=co+1;
    end
    if co==leftlink
        break;
    end
end

%%%%2. stars randomly distributed in the range [1,L]with intergrate number
%get the strength connection adjacency
for i=1:L                                  %%this part can optimaze!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    for j=1:L
        strength1=power( abs((i-j)/(L-1)) , -alfa);
        maxstrength=power(1, -alfa);
        strength=(strength1-maxstrength);
        strengthadjacency(i,j)=strength;
        strengthadjacency(j,i)=strength;
    end
end


%for the agent2star information stor
initialagent = round( L/2+ agentwidth.*randn(n,1));
%%%%%%%%%%%%userandperm to general no same position star ranged [1,L]
norepeatrstar=randperm(L);
[numagentinStar,starposition]=hist(initialagent,sort( norepeatrstar(1:N) ));  %the number N must small than L, otherwise error.


% % the uniform distribution
% initialagent=randi(N,n,1);
% norepeatrstar=randperm(L);
% [numagentinStar,starposition]=hist(initialagent,sort( norepeatrstar(1:N) ));


%%%inital configuration
% fprintf(fid2,'%.5g    \n',starposition);
averposition=unique(starposition);
fprintf(fid7,'%.5g  \n', sum(averposition())/length(averposition));
fprintf(fid8,'%.5g  \n', sum(sum(agentadjacency))/n );
agentavgopin=0;
for i=1:length(numagentinStar)
    agentavgopin=agentavgopin+numagentinStar(i)*starposition(i)/n;
end
fprintf(fid12,'%.5g  \n', agentavgopin );
fprintf(fid13,'%.5g  \n', sum(starposition)/length(starposition));


for i=1:length(starposition)
    fprintf(fid10,'%.5g   %.5g \n ', starposition(i), numagentinStar(i));
end


getagent(1)=0;
getagent(2:length(numagentinStar)+1)=cumsum(numagentinStar);
for i=1:length(numagentinStar)
    agent2star(i, getagent(i)+1:getagent(i+1))=1;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%the dynamics of the system%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:timestep
    starage=starage+1;
    fprintf('t is  %.5g \n ',t );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%1.dynamics of the agent%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %%%%%%%%a. randomly chose a agent randa and add a link to other agent with probability p_alink
    randa=randi(n);
    neighbranda=find(agentadjacency(randa,:)==1); 
    inda1=find(agentadjacency(randa,:)~=1);
    if length(inda1)>1
        choseneighborranda=[];
        choseneighborranda(1:length(inda1))=inda1; %the node can be connected
        choseneighborranda(:,choseneighborranda()==randa)=[] ;%deleted its self
        
        
        starofranda=find( agent2star(:,randa)==1);
        randaposition=starposition(starofranda);
        
        multiply=[];
        chosennodeset=[];
         
        for i=1:length(choseneighborranda )
           
            neighbchoseneighborranda=[];
            
            degreenode=length( find(agentadjacency(choseneighborranda(i),:)==1)  );
            neighbchoseneighborranda=find(agentadjacency(choseneighborranda(i),:)==1);
            commonneighb=length( intersect( neighbchoseneighborranda,neighbranda )); %%the common neighb of the the agent pair. 
            
            starofrandaneig=find(agent2star(:,choseneighborranda(i))==1);
            starinthisposition=[];
            temp1=[];
            agentinthisposition=[];
            randaneigposition=starposition(starofrandaneig);
            starinthisposition=find(starposition==randaneigposition);
            agentinthisposition=find(agent2star(starinthisposition,:)'==1);  
            temp1=mod(agentinthisposition,n);
            temp1(find(temp1==0))=n;
            commonagentforthisposition=length(intersect( neighbranda, temp1));
           
            neighbfact=max(1,commonneighb+commonagentforthisposition );
            multiply(i)=degreenode*strengthadjacency( starposition(starofranda),starposition(starofrandaneig) )*neighbfact;
            if ( starposition(starofranda)==starposition(starofrandaneig) )   %%%find the agent in the same star to connect
                chosennodeset(i)=choseneighborranda(i);
            end
        end
        chosennodeset(find(chosennodeset==0))=[];%%get ride of the case filled by 0
        normalizepalink=sum(multiply);
        
        %%%%which node to chose to add link
        %%case 1. the normalizepalink is infinity which means their existed the nodes that belong to the same star as randa
        if normalizepalink==Inf
            chosennode=chosennodeset(randi(length(chosennodeset)));
            agentadjacency(randa,chosennode)=1;  %adding a link
            agentadjacency(chosennode,randa)=1;
            fprintf(fid19, '%.5g  \n ',  abs(starposition(find(agent2star(:,randa)==1))- starposition(find(agent2star(:,chosennode)==1))) );
        end
        
        %%case 2. the normalizepalink is not infinity which chose a node according to the probability
        if normalizepalink~=Inf
            for i=1:length(choseneighborranda )
                if  rand < multiply(i)/( normalizepalink )
                    agentadjacency( randa  , choseneighborranda(i)  )=1; %adding a link
                    agentadjacency( choseneighborranda(i),   randa  )=1;
                    break;
                end
            end
            fprintf(fid19, '%.5g  \n ', abs(starposition(find(agent2star(:,randa)==1))- starposition(find(agent2star(:, choseneighborranda(i))==1))) );
        end
        fprintf(fid18, '%.5g  \n ', t-addtime);
        addtime=t;
    end
    
    %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%  %%%%
    %%%%%% b. randomly chose a agent and breaking one of its link with probability p_blink
    inda3=[];
    neighborofstrength=[];
    neighbofrandb=[];
    randb=randi(n);
    inda3=find(agentadjacency(randb,:)==1);
    if length(inda3) >2    %to advoice the isolated node with no link to break
        neighbofrandb(1:length(inda3))=inda3; %the neighbor node of randb
        starofrandb=find(agent2star(:,randb)==1);  % the star that randb belong to
        starpositionrandb=starposition(starofrandb )  ;  %get the position of starrandb
        
        
        chosenbreak=[];
        for i=1:length(neighbofrandb)
            starpositionrandbneig=starposition( find(agent2star(:,neighbofrandb(i))==1 ));
            neighofrandbi=find(agentadjacency(:,neighbofrandb(i))==1);
            commonneighb=length( intersect(neighofrandbi,neighbofrandb));
            
             starinthispostition=[];
             agentinthisposition=[];
            starinthispostition=find(starposition==starpositionrandbneig);
            agentinthisposition=find(agent2star(starinthispostition,:)'==1);  
            temp1=mod(agentinthisposition,n);
            temp1(find(temp1==0))=n;
            commonagentforthisposition=length(intersect(neighbofrandb, temp1));
            
            neighbfact=max(1,commonneighb+commonagentforthisposition );
            
            if strengthadjacency( starpositionrandb, starpositionrandbneig )==0
                chosenbreak(i)=neighbofrandb(i);
            end
            neighborofstrength(i) = 1/( (strengthadjacency( starpositionrandb, starpositionrandbneig))* neighbfact); %the normalizepblink can also be Inf
        end
        chosenbreak(find(chosenbreak==0))=[];
        
        normalizepblink=sum( neighborofstrength );
        
        % %         case 1. the normalizepblink is infinity which chose the node far away
        if normalizepblink==Inf
            choseb=chosenbreak(randi(length(chosenbreak)));
            agentadjacency(randb  , choseb )=0; %breaking a link
            agentadjacency(choseb,   randb )=0; %breaking a link
            fprintf(fid17, '%.5g  \n ', abs( starposition(find(agent2star(:,randb)==1)) - starposition(find(agent2star(:,choseb)==1)))  );
        end
        
        % %         case 2. the normalizepblink is not infinity which chose the node
        if normalizepblink~=Inf
            for i=1:length(neighbofrandb)
                if rand < neighborofstrength(i)/(normalizepblink)
                    agentadjacency(randb  , neighbofrandb(i) )=0; %breaking a link
                    agentadjacency(neighbofrandb(i),   randb )=0; %breaking a link
                    fprintf(fid17, '%.5g  \n ', abs( starposition(find(agent2star(:,randb)==1)) - starposition(find(agent2star(:,neighbofrandb(i))==1)))  );
                    break;
                end
            end
        end
        fprintf(fid16, '%.5g  \n ', t-breaktime);
        breaktime=t;
    end %%the end of the b
    
    
    %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
    %%%%%% c. randomly chose a agent and rewiring its link with probability p_rlink
    starofrandc=[];
    positioncanbechosen=[];
    inda4=[];
    multiply=[];
    averP=[];
    chosenposition=[];
    inda4=[];
    chosenposition=[];
    neighbrandc=[];
    
    
    averP=unique(starposition);
    averageO=sum(averP)/ length(averP);
    
    randc=randi(n);
    starofrandc=find( agent2star(:,randc)==1 );
    positionofstarrandc=starposition(starofrandc );
    
    positioncanbechosen=starposition;
    positioncanbechosen(find(positioncanbechosen==positionofstarrandc))=[];
    
    inda4=unique(positioncanbechosen);  %the chosen star is the stars that not in the same position
    chosenposition=inda4;
    
    inda5=[];
    positioncontainstar=[];
    
    neighbrandc=find(agentadjacency(:,randc)==1);
    
    
    if length(chosenposition)>0
        
        for i=1:length(chosenposition)%%get the probability for each position
            %find the agent number in this position stars
            inda5=find(starposition==chosenposition(i));
            positioncontainstar(1:length(inda5) )=inda5;
            %%%%%           
    
            clear starofagent temp1
            starofagent=find(agent2star(positioncontainstar,:)'==1);
            temp1=mod(starofagent,n);
            temp1(find(temp1==0))=n;
            
            agentinthisposition=length(intersect( neighbrandc, temp1));
            neighbfact=max(1, agentinthisposition);
            %%%%%%%%%%
            stardegree=sum(sum(agent2star(positioncontainstar,:),2));
            elo=max(speakelo(positioncontainstar) );
            thermoeff=exp((Thermo-averageO)*chosenposition(i)/((averageO*Thermo)));
            
            %  multiply(i)=stardegree*strengthadjacency(positionofstarrandc,chosenposition(i) )*elo*thermoeff;
            multiply(i)=stardegree*strengthadjacency(positionofstarrandc,chosenposition(i) )*elo*neighbfact;
        end
        
        normalizepclink=sum( multiply );
        starcancon=[];
        for i=1:length(chosenposition)
            if rand < ( multiply(i)/normalizepclink )
                starofrandc=find(agent2star(:,randc)==1);
                agent2star( :, randc)=0;
                starcancon=find(starposition==chosenposition(i));
                agent2star(starcancon( randi( length( starcancon) )),randc)=1;
                fprintf(fid15, '%.5g  \n ', -starposition(starofrandc)+starposition(starcancon( randi( length( starcancon) ))) );
                break;
            end
        end
        
        
        fprintf(fid14, '%.5g  \n ', t-rewiretime);
        rewiretime=t;
        
    end
    
    if length(find(sum(agent2star)~=1))>0
        fprintf('the agent2star is error in line 305');
        return
    end
    
    %%%%%%%%%%%%%%%%%%%%%%the end of the dynamics agenttime%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%1.dynamics of the star%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%get the information of star merging
    
    if(mod(t,agenttime)==0)
        positionstarSP1=[];
        positionstarSP1=unique(starposition);
        averageO=sum(positionstarSP1)/ length(positionstarSP1);
        
        %%%%%%%%%a. randomly chose a star position and merging to other star position with probability p_mstar
        if length(positionstarSP1)>1
            degreeofpositionSP1=[];
            mergingstar=[];
            for i=1:length(positionstarSP1) %get the information of all star position
                mergingstar= find( starposition==positionstarSP1(i) );%the number of stars in this position
                positiondegreeSP1=sum(sum(agent2star( mergingstar,:),2));
                degreeofpositionSP1(i)=positiondegreeSP1;
            end
            
            randstarposition=positionstarSP1( randi(length(positionstarSP1))); %chose a star position
            degreeofrandstarposition= degreeofpositionSP1( find(positionstarSP1==randstarposition)  );
            
            probmerg=[];
            for i=1:length(positionstarSP1)
                if strengthadjacency( randstarposition ,positionstarSP1(i) )~=Inf
                    probmerg(i)=strengthadjacency( randstarposition , positionstarSP1(i) )*( degreeofpositionSP1(i)+degreeofrandstarposition);
                else
                    probmerg(i)=0;
                end
            end
            normalizepmstar=sum( probmerg );
            
            for i=1:length(positionstarSP1)
                if rand < probmerg(i)/normalizepmstar
                    stpage=starage(find(starposition==positionstarSP1(i)));
                    strage=starage(find(starposition==randstarposition ) );
                    fprintf(fid23,'%.5g  \n', stpage(1));
                    fprintf(fid23,'%.5g  \n', strage(1));
                    
                    fprintf(fid3,'%.5g  %.5g  \n',t, positionstarSP1(i)- randstarposition );
                    fprintf(fid4,'%.5g  \n',t-timeintermerg);
                    timeintermerg=t;
                    %change the position of staradjacency
                    % %                     if ( mod( (positionstarSP1(i)+ randstarposition)/2, 1)==0 )
                    % %                         newposition=(positionstarSP1(i)+ randstarposition)/2 ;%change the position of staradjacency
                    % %                     else
                    % %                         setleft=rand;
                    % %                         if round(setleft)==1
                    % %                             newposition=(positionstarSP1(i)+ randstarposition-1)/2;
                    % %                         else
                    % %                             newposition=(positionstarSP1(i)+ randstarposition+1)/2;
                    % %                         end
                    % %                     end
                    % %
                    
                    % % for the new version of the position
                    totaldegree=degreeofpositionSP1(i)+degreeofrandstarposition;
                    
                    if positionstarSP1(i)<randstarposition
                        newposition=positionstarSP1(i)+round( (randstarposition-positionstarSP1(i))* degreeofpositionSP1(i)/totaldegree );
                    end
                    
                    if positionstarSP1(i) > randstarposition
                        newposition=randstarposition+round( (positionstarSP1(i)-randstarposition)* degreeofrandstarposition/totaldegree );
                    end
                    
                    starposition(find(starposition==positionstarSP1(i)))=newposition;
                    starposition(find(starposition==randstarposition ))=newposition;
                    
                    
                    starage( find(starposition==positionstarSP1(i)) )=1;%change the age of star
                    starage( find(starposition==randstarposition )  )=1;
                    
                    break;
                end
            end
        end
        
        if length(find(sum(agent2star)~=1))>0
            fprintf('the agent2star is error in line 305');
            return
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%b.one star splitting into two stars j and k with probability p_sstar%%%%%%%%
        
        averageO=sum(positionstarSP1)/ length(positionstarSP1);
        splitposition=[];
        positionstarSP2=[];
        indb1=[];
        positionstarSP2=unique(starposition);
        
        k=1;
        for i=1:length(positionstarSP2)    %get  the position contain more than two stars
            indb1=find(starposition==positionstarSP2(i));
            if(length(indb1)>1)
                splitposition(k)=positionstarSP2(i);
                k=k+1;
            end
        end
        
        degreeofpositionSP2=[];
        probsplit=[];
        strengthnewposition=[];
        splitpositionstar=[];
        if length(splitposition)>0
            for i=1:length(splitposition) %get the probability of each position
                splitpositionstar= find( starposition==splitposition(i) );
                splitpositionage=starage(splitpositionstar(1)); %the age of stars in this position !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! wrong? wrong until checked by 23/06/2010
                degreespositionSP2=sum(sum(agent2star(splitpositionstar,:),2));
                degreeofpositionSP2(i)=degreespositionSP2;
                probsplit(i)= degreeofpositionSP2(i)*splitpositionage;
            end
            
            normalizeps=sum(probsplit);
            
            for i=1:length(splitposition)
                if rand < probsplit(i)/normalizeps;  %spliting beginning at the position of splitposition(i)
                    motherstarposition=splitposition(i);  %the motherstarposition
                    %get the new position porbibility from motherposition
                    for ii=1:L
                        if ii==motherstarposition
                            strengthnewposition(ii)=0;
                        else
                            % strengthnewposition(ii)=strengthadjacency(motherstarposition,ii)*exp(ii*(Thermo-averageO)/((Thermo*averageO)) );
                            strengthnewposition(ii)=strengthadjacency(motherstarposition,ii);
                        end
                    end
                    
                    normalizepositionprob=sum( strengthnewposition( ) );
                    a=0;
                    while(a<1)   %get the new position
                        for ii=1:L
                            if rand < strengthnewposition(ii)/normalizepositionprob %the new position is ii
                                fprintf(fid5,'%.5g  %.5g  \n ',t,ii-motherstarposition );
                                fprintf(fid6,'%.5g  \n',t-timeinterspli);
                                timeinterspli=t;
                                newbornposition=ii;
                                a=a+1;
                                break;
                            end
                        end
                    end
                    
                    
                    %%chose one of the stars in motherposition(splitposition(i)) as the newbornstar
                    splitstarset=[];
                    splitstarset=find( starposition==splitposition(i) );
                    newbornstar=splitstarset(randi(length(splitstarset))); %the star that need to change position
                    
                    starposition(newbornstar)=newbornposition;  %the star changes position
                    fprintf(fid23,'%.5g  \n',starage(newbornstar));
                    fprintf(fid23,'%.5g  \n',starage(splitstarset(1)));
                    starage(newbornstar )=1; %the star changes age
                    starage(splitstarset)=1;  %the simpler way
                    
                    
                    
                    %%%%%get off the agent in the newborn star
                    splitstarset1=[];
                    agentinnewbornstar=[];
                    
                    
                    splitstarset1=splitstarset;   %%splitstarset1 is the star in the position of motherposition
                    splitstarset1(:,splitstarset==newbornstar)=[];
                    starreconnagent=splitstarset1(randi(length(splitstarset1)));
                    agentinnewbornstar=find(agent2star(newbornstar,:)==1);
                    for ii=1:length(agentinnewbornstar)
                        agent2star(starreconnagent,agentinnewbornstar(ii))=1;
                        agent2star(newbornstar,agentinnewbornstar(ii))=0;
                    end
                    y=[];
                    x=[];
                    [x,y]=find(agent2star(splitstarset1,:)==1);
                    agentnum=unique(y);
                    if length(agentnum) ~= length(y)
                        fprintf('The   length(agentnum)< length(y) \n\n');
                    end
                    numerlink2rece=round(sum(sum(agent2star(splitstarset1,:),2))*percenttogive);
                    if numerlink2rece >1
                        for t123=1: numerlink2rece
                            res=find(agent2star(:,y(t123))==1);
                            agent2star(res,y(t123))=0;
                            agent2star(newbornstar, y(t123))=1;
                        end
                    end % numerlink2rece >1
                    linksplit=linksplit+numerlink2rece;
                    
                    
                    %%%%%end of get off the agent in the newborn star
                    positionstarSP21=[];
                    
                    positionstarSP21=positionstarSP2;
                    positionstarSP21(:,find(positionstarSP2==motherstarposition))=[];
                    
                    %the star get links form other stars positioned from
                    %starnowposition(must have more than one links)
                    elo=speakelo(newbornstar);
                    for ii=1:length(positionstarSP21)
                        indexstar=[];;
                        indexstar=find(starposition==positionstarSP21(ii));
                        if rand < (strengthadjacency(positionstarSP21(ii),newbornposition)*elo*sum(sum(agent2star(indexstar,:),2)))/(n*10*sum(strengthadjacency(1,2:L)))
                            for iii=1:length(indexstar)
                                numberlik2g=round(length(find(agent2star(indexstar(iii),:)==1))*percenttogive);
                                if numberlik2g > 1
                                    link=[];
                                    agentremove=[];
                                    link=find(agent2star(indexstar(iii),:)==1);
                                    agentremove=link(1: numberlik2g);
                                    for iiii=1:length(agentremove)
                                        agent2star(indexstar(iii), agentremove(iiii))=0;
                                        agent2star(newbornstar, agentremove(iiii))=1;
                                    end
                                end
                                linksplit=linksplit+numberlik2g;
                            end% for each star in indexstar (each star in the same position)to give its link if possible
                        end %%for if condition
                    end %%for each position
                    fprintf(fid24,'%.5g  \n', linksplit );
                    linksplit=0;
                    break; % no continue of the line 455 for i=1:length(splitposition)
                end %%for the end of the spliting(i)
            end %for the star spliting
        end%for the star spliting if there is the position have more than 2 stars
    end %the end of the star dynamics
    
    
    if length(find(sum(agent2star)~=1))>0
        fprintf('the agent2star is error in line 505');
        return
    end
    
    
    
    fprintf(fid1,'  %.5g\n',length(unique(starposition)));
    fprintf(fid2,'%.5g    \n',starposition);
    positionstarSP1=[];
    positionstarSP1=unique(starposition);
    averageO=sum(positionstarSP1)/ length(positionstarSP1);
    fprintf(fid7,'%.5g  \n', averageO);
    fprintf(fid8,'%.5g  \n', sum(sum(agentadjacency))/n );
    if (mod(t,agenttime)==0 && t~=0)
        degreeall=sum(agentadjacency);
        al=length(find(degreeall~=0));
        for i=1:n
            fprintf(fid9,'%.5g   %.5g  \n',i,length(find(degreeall==i))/al  );
        end
    end
    
    agentavgopin=0;
    degree=sum(agentadjacency);
    for st=1:length(positionstarSP1)
        starnumSP1=find(starposition==positionstarSP1(st));
        positiondegreeSP1=sum(sum(agent2star(starnumSP1,:),2));
        fprintf(fid10,' %.5g  %.5g \n ', positionstarSP1(st), positiondegreeSP1);
        agentavgopin=agentavgopin+positionstarSP1(st)*positiondegreeSP1/(n);
        spagent=[];
        for sp=1:length(starnumSP1)
            spagent(length(spagent)+1:length(spagent)+length(find(agent2star(starnumSP1(sp),:)==1)))=find(agent2star(starnumSP1(sp),:)==1);
        end
        %         fprintf(fid20,'%.5g  %.5g \n', positionstarSP1(st), (sum(degree(spagent))/length(spagent))/(sum(degree)/n )  );
        
        countinterconnect=0;
        countouterconnect=0;
        x=[];
        neigsp=[];
        neigspU=[];
        [x,neigsp]=find(agentadjacency(spagent,:)==1);
        neigspU=unique(neigsp);
        %         fprintf(fid25,'%.5g  %.5g \n', positionstarSP1(st), length(spagent)/length(neigspU )  );
        
        %         for connet1=1:length(spagent)
        %             for connet2=connet1+1:length(spagent)
        %            indxnumin=(spagent(connet1)-1)*n+spagent(connet2);
        %             end
        %         end
        %         countinterconnect=length(find(agentadjacency(indxnumin)==1));
        %         for connet3=1:length(neigspU)
        %             for connet4=connet3+1:length(neigspU)
        %                 indxnumall=(neigspU(connet3)-1)*n+neigspU(connet4);
        %             end
        %         end
        %          countallconnect=length(find(agentadjacency(indxnumall)==1));
        %         fprintf(fid21,'%.5g  %.5g \n', positionstarSP1(st), countinterconnect/countallconnect );
        
        fprintf(fid22,'%.5g  %.5g \n', positionstarSP1(st), max(speakelo(starnumSP1)));
    end
    
    
    for ttt=length(positionstarSP1)+1:N
        fprintf(fid10,'%.5g  %.5g\n ',0,0);
        %         fprintf(fid20,'%.5g  %.5g \n',0,0);
        %         fprintf(fid21,'%.5g  %.5g \n',0,0);
        fprintf(fid22,'%.5g  %.5g \n', 0,0);
        %         fprintf(fid25,'%.5g  %.5g \n', 0,0);
    end
    
    fprintf(fid13,'%.5g  \n', sum(starposition)/length(starposition));
    fprintf(fid12,' %.5g\n ',agentavgopin);
    
    
    
    
    
    
end %the end of the whole dynamics

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% degree,richness distribution %%%%%%%%%%%%%%%%%%%%%%
fclose(fid1);
fclose(fid2);
fclose(fid3);
fclose(fid4);
fclose(fid5);
fclose(fid6);
fclose(fid7);
fclose(fid8);
fclose(fid9);
fclose(fid10);
fclose(fid12);
fclose(fid13);
fclose(fid14);
fclose(fid15);
fclose(fid16);
fclose(fid17);
fclose(fid18);
fclose(fid19);
fclose(fid20);
fclose(fid21);
fclose(fid22);
fclose(fid23);
fclose(fid24);
fclose(fid25);
toc
fprintf(fid11,'The paramaters are set as fowllowing (with newlink from mother star and  pc=0.7; the rescaled S function \n\n');

fprintf(fid11,'n=  %.5g  (agent number) \n', n);
fprintf(fid11,'N=  %.5g  (star number, we set range of star''s position is equal to L (N=L) ) \n', N);
fprintf(fid11,'Alfa=  %.5g  (the strength exponent ) \n', alfa);
fprintf(fid11,'timestep=  %.5g  (the total timesteps and we set agent evolution time is t=50 ) \n',timestep);
%fprintf(fid11,'Thermo=  %.5g  (the temperature of the Thermo ) \n',Thermo);
fprintf(fid11,'No Thermo  effect (the temperature of the Thermo ) \n');
fprintf(fid11, 'agentwidth=  %.5g  \n', agentwidth);
fprintf(fid11, 'agenttime=  %.5g (the number of agents''evoltion) \n', agenttime);
fprintf(fid11, 'percenttogive=  %.5g (the percent of links to give when a new stars splitting) \n',percenttogive);
fprintf(fid11,'the program running time is :  %.5g  hours \n\n', toc/3600);
fprintf(fid11,'the program finishes at:   %.5g / %.5g  / %.5g/   %.5g:   %.5g:   %.5g.   \n',clock );
fclose(fid11);


cd 'C:\WANGRU\10.3.1ISMANS\work\2popinion model1\11.27\'
plotfig(n,N,L,alfa,timestep,agentwidth,agenttime,Thermo)
