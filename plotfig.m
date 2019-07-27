function plotfig(n,N,L,alfa,timestep,agentwidth,agenttime,Thermo)
figtime=10;


cd 'C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring'
dos(['mkdir n' num2str(n) 'N' num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo)])



%%%%%%%%%fid1
numstar=load(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid1numberofstar''.txt']);  %the evolution of number of stars with timestep of star

y1=reshape(numstar,agenttime,[]);
figure
hold on
plot(0,N,'k:o','MarkerSize',2);
plot(1:size(y1,2), y1(agenttime,:) ,'k:o','MarkerSize',2);
xlabel(['time of star *' num2str(agenttime)  ] )

ylabel('number of star')
str=['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\1fid1numberofstar_n' num2str(n) 'N' num2str(N) 'ts' num2str(timestep) '.jpg'];
saveas(gcf,str);
close all


%fid2
positstar=load(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid2positionofstar.txt']);                   %the evolution of number of stars with timestep of star
xx=reshape(positstar,N,[]);
for i=1:100*agenttime:size(xx,2)
    figure;
    for j=1:N
        xu=unique(xx(:,i));
        p=hist(xx(:,i),xu);
        plot(xu,p,'ko-','MarkerSize',2);
    end
    
    xlabel('star position')
    ylabel('frequence')
    title(['t=',num2str(i-1)])
    str=['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
        num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
        num2str(agenttime)  'Thermo' num2str(Thermo) '\fid2positionofstar_n' num2str(n) 'N' num2str(N) 'ts' num2str(timestep) 't=' num2str(i) '.jpg'];
    saveas(gcf,str);
    hold on
end
close all
% 
% 
% 
% %%%%%%%%%%fid3
distanceofmergstar=load(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid3distanceofstarmerg' '.txt']);                   %the distrance of stars can merging
dis=distanceofmergstar(:,2);
disrange=unique(dis);
disranges=sort(disrange);
for i=1:length(disranges)
    y3(i)=length(find(dis==disranges(i)));
    x3(i)=disranges(i);
end
figure;
hold on
plot(x3/L,y3/length(dis),'k:o','MarkerSize',2)
xlabel('rescaled star merge distance(selected - rand)')
ylabel('probability ')
str=['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid3distanceofstarmerg_n' num2str(n) 'N' num2str(N) 'ts' num2str(timestep) '.jpg'];
saveas(gcf,str);

close all



%%%%%%%%%%%fid4
intervalofmergstar=load(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid4intervalofstarmerg'  '.txt']);     %the interval time of star can merging
intervalrange=unique(intervalofmergstar);
 y4=[];
for i=1:length(intervalrange)
    y4(i)=length(find(intervalofmergstar==intervalrange(i)));
    x4(i)=intervalrange(i);
end

figure
hold on
plot(x4/agenttime,y4/length(intervalofmergstar),'k:o','MarkerSize',2)
xlabel('star merging time interval')
ylabel('probability ')
str=['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid4intervalofstarmerg_n' num2str(n) 'N' num2str(N) 'ts' num2str(timestep) '.jpg'];
saveas(gcf,str);
close all


%%%%%%%%%%%fid5
distanceofsplitstar=load(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid5distanceofstarsplit' '.txt']);                   %the distrance of stars can merging
dis5=distanceofsplitstar(:,2);
disrange5=unique(dis5);
disranges5=sort(disrange5);
for i=1:length(disranges5)
    y5(i)=length(find(dis5==disranges5(i)));
    x5(i)=disranges5(i);
end
figure
hold on
plot(x5/L,y5/length(dis),'k:o','MarkerSize',2)
xlabel('rescaled star spliting distance(newborn-mother)')
ylabel('probability')
str=['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid5distanceofstarsplit_n' num2str(n) 'N' num2str(N) 'ts' num2str(timestep) '.jpg'];
saveas(gcf,str);
close all


%%%%%%%%%%%fid6
intervalofsplitstar=load(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid6intervalofstarsplit' '.txt']);     %the interval time of star can merging
intervalrange=unique(intervalofsplitstar);
for i=1:length(intervalrange)
    y6(i)=length(find(intervalofsplitstar==intervalrange(i)));
    x6(i)=intervalrange(i);
end
figure
hold on
plot(x6/agenttime,y6/length(intervalofsplitstar),'k:o','MarkerSize',2)
xlabel('star spliting time interval')
ylabel('probability')
str=['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid6intervalofstarsplit_n' num2str(n) 'N' num2str(N) 'ts' num2str(timestep) '.jpg'];
saveas(gcf,str);
close all




%%%%%%%%%%fid7
averopin=load(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid7aveopinion' '.txt']);                   %the evolution of average opinion with timestep of star
figure
hold on
plot(0,averopin(1)/L,'k:o','MarkerSize',2);
averopin(1)=[];
y7=reshape(averopin,agenttime,[]);
plot(1:size(y7,2), y7(agenttime,:)/L ,'k:o','MarkerSize',2);
xlabel(['time of star *' num2str(agenttime)  ] )
ylabel('position average Opinion')
str=['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\2fid7aveopinion_n' num2str(n) 'N' num2str(N) 'ts' num2str(timestep) '.jpg'];
saveas(gcf,str);
close all


%%%%%%%%%%fid13 staraverage
averopin=load(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid13staraveopinion' '.txt']);                   %the evolution of average opinion with timestep of star
figure
hold on
plot(0,averopin(1)/L,'k:o','MarkerSize',2);
averopin(1)=[];
y7=reshape(averopin,agenttime,[]);
plot(1:size(y7,2), y7(agenttime,:)/L ,'k:o','MarkerSize',2);
xlabel(['time of star *' num2str(agenttime)  ] )
ylabel('star average Opinion')
str=['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\2fid13staraveopinion_n' num2str(n) 'N' num2str(N) 'ts' num2str(timestep) '.jpg'];
saveas(gcf,str);
close all


%%%%%%%%%%fid12 agentaverage
averopin=load(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid12agentaveopinion'  '.txt']);                   %the evolution of average opinion with timestep of star
figure
hold on
plot(0,averopin(1)/L,'k:o','MarkerSize',2);
averopin(1)=[];
y7=reshape(averopin,agenttime,[]);
plot(1:size(y7,2), y7(agenttime,:)/L ,'k:o','MarkerSize',2);
xlabel('time')
ylabel('agent average Opinion')
str=['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\2fid12agentaveopinion_n' num2str(n) 'N' num2str(N) 'ts' num2str(timestep) '.jpg'];
saveas(gcf,str);
close all


%%%%%%%%%%%fid8
averagedegree=load(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid8avedegree' '.txt']);  % the evolution of averagedegree with timestep of agent
y8=averagedegree(:,1);
x8=(1:length(averagedegree))';
figure
hold on
plot(x8,y8,'k:o','MarkerSize',2);
xlabel('time')
ylabel('average degree')
str=['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\3fid8avedegree_n' num2str(n) 'N' num2str(N) 'ts' num2str(timestep) '.jpg'];
saveas(gcf,str);
close all



%%%%%%%fid9
degree=load(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid9degreedistribution''.txt']);       %degree  probability distribution
x=degree(:,1);
y=degree(:,2);
xx=reshape(x,n,[]);
yy=reshape(y,n,[]);
for i=1:50*figtime:size(xx,2)
    figure;
    plot(xx(:,i),yy(:,i),'k:o','MarkerSize',2);
    xlabel('degree')
    ylabel('probability')
    title(['t= (*2)',num2str(i)])
    hold on
    str=['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
        num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
        num2str(agenttime)  'Thermo' num2str(Thermo) '\fid9degreedistribution_n' num2str(n) 'N' num2str(N) 'ts' num2str(timestep) 't=' num2str(i) '.jpg'];
    saveas(gcf,str);
    close all
end




%%%%%%%%%%fid10
agentinstar=load(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid10numberofagentinstar'  '.txt']);
eachtimestar=reshape(agentinstar,N,[]);

for i=1:50*figtime:size(eachtimestar,2)/2
    figure;
    ind=find(eachtimestar(:,i)~=0);
    
    plot( eachtimestar(ind,i)/N, eachtimestar(ind,i+size(eachtimestar,2)/2)/n,'k:o')
    xlabel('rescaled star position');
    ylabel('percent of agent in this star');
    title(['t=',num2str(i-1)]);
    hold on
    str=['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
        num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
        num2str(agenttime)  'Thermo' num2str(Thermo) '\fid10numberofagentinstar_n' num2str(n) 'N' num2str(N) 'ts' num2str(timestep) 't=' num2str(i-1) '.jpg'];
    saveas(gcf,str);
    close all
end



% agentinstar=load(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
%     num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
%     num2str(agenttime)  'Thermo' num2str(Thermo) '\fid10numberofagentinstar_n'  num2str(n) ...
%     'N' num2str(N) 'ts' num2str(timestep) '.txt']);
% eachtimestar=reshape(agentinstar,N,[]);
% for i=size(eachtimestar,2)/2+1:10*figtime:size(eachtimestar,2)
%     figure;
%     ind=find(eachtimestar(:,i)~=0 );
%     nonzeroeachtimestar=eachtimestar(ind,i);
%     for j=1:N
%         xnonzeroeachtimestar=unique(nonzeroeachtimestar);
%         pro=hist(nonzeroeachtimestar,xnonzeroeachtimestar);
%         plot(xnonzeroeachtimestar,pro,'ko-','MarkerSize',2);
%     end
%     xlabel('number of agent in one star(votes)');
%     ylabel('number of star(number of candidates)');
%     title(['t=',num2str(i-1-size(eachtimestar,2)/2)]);
%     hold on
%     str=['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
%         num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
%         num2str(agenttime)  'Thermo' num2str(Thermo) '\fid10-1probability of number in star_n' num2str(n) 'N' num2str(N) 'ts' num2str(timestep) 't=' num2str(i-(size(eachtimestar,2)/2)-1) '.jpg'];
%     saveas(gcf,str);
% end
% close all


%%%%%%%%%%fid14
intervalofrewire=load(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid14rewiringtime' '.txt']);                   %the distrance of stars can merging
interrange=unique(intervalofrewire);
for i=1:length(interrange)
    y14(i)=length(find(intervalofrewire==interrange(i)));
end
x14=(1:length(interrange));
figure;
hold on
plot(x14,y14/length(intervalofrewire),'k:o','MarkerSize',2)
xlabel('agent rewiring time inter')
ylabel('probability')
str=['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid14rewiringtime''.jpg'];
saveas(gcf,str);
close all


%%%%%%%%%%fid15
disrewire=load(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid15rewiringdis' '.txt']);   
dis=unique(disrewire);
diss15=sort(dis);
for i=1:length(diss15)
    yfid15(i)=length(find(disrewire==diss15(i)));
    xfid15(i)=diss15(i);
end
figure
hold on
plot(xfid15/L,yfid15/length(disrewire),'k:o','MarkerSize',2 )
xlabel('rescaled agent rewiring distance(rewiring to - origion)')
ylabel('probability')
str=['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid15rewiringdis_n' num2str(n) 'N' num2str(N) 'ts' num2str(timestep) '.jpg'];
saveas(gcf,str);
close all

%%%%%%%%%%fid16
intervalofbreak=load(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid16breakingtime' '.txt']);                   %the distrance of stars can merging
interb=intervalofbreak;
interrangeb=unique(interb);
for i=1:length(interrangeb)
    y16(i)=length(find(interb==interrangeb(i)));
end
x16=(1:length(interrangeb));
figure;
hold on
plot(x16,y16/length(interb),'k:o','MarkerSize',2)
xlabel('agent breaking link time inter')
ylabel('probability ')
str=['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid16breakingtime''.jpg'];
saveas(gcf,str);
close all


%%%%%%%%%%fid17
disbreak=load(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid17breakingdis' '.txt']);   
disb=unique(disbreak);
disbs17=sort(disb);
for i=1:length(disbs17)
    yfid17(i)=length(find(disbreak==disbs17(i)));
    xfid17(i)=disbs17(i);
end
figure
hold on
plot(xfid17/L,yfid17/length(disbreak),'k:o','MarkerSize',2 )
xlabel('agent breaking link''s star rescaled distance')
ylabel('probability')
str=['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid17breakingdis_n' num2str(n) 'N' num2str(N) 'ts' num2str(timestep) '.jpg'];
saveas(gcf,str);
close all




%%%%%%%%%%fid18
interadd=load(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid18addingtime' '.txt']);   
intera=unique(interadd);
for i=1:length(intera)
    yfid18(i)=length(find(interadd==intera(i)));
end
xfid18=(1:length(intera));
figure
hold on
plot(xfid18,yfid18/sum(yfid18),'k:o','MarkerSize',2 )
xlabel('agent adding link time interval')
ylabel('probability')
str=['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid18addingtime_n' num2str(n) 'N' num2str(N) 'ts' num2str(timestep) '.jpg'];
saveas(gcf,str);
close all


%%%%%%%%%%fid19
disadd=load(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid19addingdis' '.txt']);   
disa=unique(disadd);
c19=0;
for i=1:length(disa)
    yfid19(i)=length(find(disadd==disa(i)));
    c19=c19+yfid19(i);
end
xfid19=(1:length(disa));
figure
hold on
plot(xfid19/L,yfid19/c19,'k:o','MarkerSize',2 )
xlabel('rescaled agent adding link''s star distance')
ylabel('probability of adding')
str=['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\fid19addingdis_n' num2str(n) 'N' num2str(N) 'ts' num2str(timestep) '.jpg'];
saveas(gcf,str);
close all

%%%%%%%%%%fid20
% avgdegstar=load(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
%     num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
%     num2str(agenttime)  'Thermo' num2str(Thermo) '\avgdegreetostar' '.txt']);   
% eachtimeavgstar=reshape(avgdegstar,N,[]);
% for i=1:50*figtime:size(eachtimeavgstar,2)/2
%     figure;
%     ind=find(eachtimeavgstar(:,i)~=0);
%     plot( eachtimeavgstar(ind,i)/L,eachtimeavgstar(ind,i+size(eachtimeavgstar,2)/2),'k:o','MarkerSize',2 )
%     xlabel('rescaled star position');
%     ylabel('average degree in this star');
%     title(['t=',num2str(i-1)]);
%     hold on
%     str=['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
%         num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
%         num2str(agenttime)  'Thermo' num2str(Thermo) '\avgdegreetostar_n' num2str(n) 'N' num2str(N) 'ts' num2str(timestep) 't=' num2str(i-1) '.jpg'];
%     saveas(gcf,str);
%     close all
% end

%%%%%%%%%%fid21
% connet=load(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
%     num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
%     num2str(agenttime)  'Thermo' num2str(Thermo) '\edgeconnectstrength' '.txt']);   
% strengcon=reshape(connet,N,[]);
% for i=1:50*figtime:size(strengcon,2)/2
%     figure;
%     ind=find(strengcon(:,i)~=0);
%     plot( strengcon(ind,i)/L,strengcon(ind,i+size(strengcon,2)/2),'k:o','MarkerSize',2)
%     xlabel('rescaled star position');
%     ylabel('edge connection strength in this position(linkinstar/linkoutstar)');
%     title(['t=',num2str(i-1)]);
%     hold on
%     str=['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
%         num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
%         num2str(agenttime)  'Thermo' num2str(Thermo) '\edgeconnectstrength_n' num2str(n) 'N' num2str(N) 'ts' num2str(timestep) 't=' num2str(i-1) '.jpg'];
%     saveas(gcf,str);
%     close all
% end

%%%%%%%%%%fid22
eloevo=load(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\eloposition' '.txt']);   
elo=reshape(eloevo,N,[]);
for i=1:50*figtime:size(elo,2)/2
    figure;
    ind=find(elo(:,i)~=0);
    
    plot( elo(ind,i)/L,elo(ind,i+size(elo,2)/2),'k:o','MarkerSize',2)
    xlabel('rescaled star position');
    ylabel('elo in this position');
    title(['t=',num2str(i-1)]);
    hold on
    str=['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
        num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
        num2str(agenttime)  'Thermo' num2str(Thermo) '\eloposition_n' num2str(n) 'N' num2str(N) 'ts' num2str(timestep) 't=' num2str(i-1) '.jpg'];
    saveas(gcf,str);
    close all
end

%%%%%%%%%%fid23
agedis=load(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\agedistribution' '.txt']);     %the  
age=unique(agedis);
for i=1:length(age)
    y23(i)=length(find(agedis==age(i)));
end
x23=(1:length(age));
figure
hold on
loglog(x23,y23/length(agedis),'k:o','MarkerSize',2)
xlabel('age')
ylabel('probability of existed age')
str=['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\agedistribution_n' num2str(n) 'N' num2str(N) 'ts' num2str(timestep) '.jpg'];
saveas(gcf,str);
close all

%%%%%%%%%%fid24
link=load(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\link2split' '.txt']);  
linkdis=unique(link);
for i=1:length(linkdis)
    y24(i)=length(find(link==linkdis(1)));
end
x24=(1:length(linkdis));
figure
hold on
plot(x24/n,y24/length(link),'k:o','MarkerSize',2)
xlabel('split link percent(number/n)')
ylabel('probability of split link number ')
str=['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
    num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
    num2str(agenttime)  'Thermo' num2str(Thermo) '\link2split_n' num2str(n) 'N' num2str(N) 'ts' num2str(timestep) '.jpg'];
saveas(gcf,str);
close all

%%%%%%%%%%fid25
% nodeconnet=load(['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
%     num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
%     num2str(agenttime)  'Thermo' num2str(Thermo) '\nodeconnectstrength' '.txt']);   
% nodestrengcon=reshape(nodeconnet,N,[]);
% for i=1:50*figtime:size(nodestrengcon,2)/2
%     figure;
%     ind=find(nodestrengcon(:,i)~=0);
%     plot( nodestrengcon(ind,i)/L,nodestrengcon(ind,i+size(nodestrengcon,2)/2),'k:o','MarkerSize',2)
%     xlabel('rescaled star position');
%     ylabel('node connection strength in this position(nodeinstar/nodeinstarconnectednode)');
%     title(['t=',num2str(i-1)]);
%     hold on
%     str=['C:\Temp\Wang\matlab data\opinionmodel\2010.11.26with neighbor rewiring\n' num2str(n) 'N' ...
%         num2str(N) 'ts' num2str(timestep) 'alfa' num2str(alfa) 'agentwidth'  num2str(agentwidth) 'agenttime' ...
%         num2str(agenttime)  'Thermo' num2str(Thermo) '\nodeconnectstrength_n' num2str(n) 'N' num2str(N) 'ts' num2str(timestep) 't=' num2str(i-1) '.jpg'];
%     saveas(gcf,str);
%     close all
% end



end
