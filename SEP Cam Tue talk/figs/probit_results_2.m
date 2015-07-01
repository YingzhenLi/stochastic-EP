%Australian 0.328±0.0127 0.325±0.0135 0.330±0.0133 -0.634±0.010 -0.631±0.009 -0.631±0.009 
%Breast 0.037±0.0045 0.034±0.0034 0.034±0.0039 -0.100±0.015 -0.094±0.011 -0.093±0.011 
%Crabs 0.062±0.0125 0.040±0.0106 0.048±0.0117 -0.290±0.010 -0.177±0.012 -0.217±0.011
%Ionos 0.126±0.0166 0.130±0.0147 0.131±0.0149 -0.373±0.047 -0.336±0.029 -0.324±0.028 
%Pima 0.242±0.0093 0.244±0.0098 0.241±0.0093 -0.516±0.013 -0.514±0.012 -0.513±0.012 
%Sonar 0.198±0.0208 0.198±0.0217 0.198±0.0243 -0.461±0.053 -0.418±0.021 -0.415±0.021

datanames = {'Australian','Breast','Crabs','Ionos','Pima','Sonar'};

ADF = [-0.634,-0.100,-0.290,  -0.373, -0.516, -0.461];
ADFEB = [0.010, 0.015, 0.010, 0.047, 0.013, 0.053];

SEP = [-0.631, -0.094, -0.177, -0.336, -0.514, -0.418];
SEPEB = [0.009, 0.011, 0.012, 0.029, 0.012, 0.021];

EP = [-0.631 -0.093  -0.217 -0.324  -0.513 -0.415];
EPEB = [0.009 0.011 0.011 0.028 0.012 0.021];


datanames2 = {'Kin8nm','Naval','Power','Protein','Wine','Year'};
ADF2  = [0.896,3.731,-2.837,-2.973,-0.968,-3.603];
SEP2 =  [1.013,4.590,-2.846,-2.961,-0.976,-3.924];
EP2 = [1.005,4.207,-2.852,-2.979,-0.958,-3.929];
ADFEB2  =  [0.006,0.006,0.009,0.003,0.014,0];
SEPEB2 = [0.011,0.014,0.008,0.003,0.013,0];
EPEB2 = [0.007,0.011,0.008,0.003,0.011,0];


tits= {'N=8192,D=8, 58MB', 'N= 11934, D=16, 147MB','N=9568, D=4, 37MB','N=45730 ,D=9, 694MB','N = 1599, D=11, 14MB','N=515340, D=90, 65107MB'};

tits= {'58MB', '147MB','37MB','694MB','14MB','65107MB'};


All = [ADF2;SEP2;EP2];
min_All = min(All);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FontName = 'Arial';
FSsm = 16;
FSmed = 16;
FSlg = 22;

lw1 = 2;
lw2 = 2;
lw3 = 2;
ms = 14;


left = 0.1;
right = 0.05;
top = 0.07;
bottom = 0.12;
hspace = 0.05;

width = (1-left-5*hspace-right)/6;
height = (1-top-bottom);

across = [hspace+width,0,0,0];
pos1 = [left,bottom,width,height];
pos2 = pos1+across;
pos3 = pos2+across;
pos4 = pos3+across;
pos5 = pos4+across;
pos6 = pos5+across;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 figure1 = figure;
    
    PP = [0,0,20.00,12.40]*0.7; %*** paper position in centimeters
    PS = PP(end-1:end); % paper size in centimeters

set(figure1,'paperpositionmode','manual','paperposition', ...
        PP,'papersize',PS, 'paperunits','centimeters');

% So the figure is the same size on the screen as when it is printed:
pu = get(gcf,'PaperUnits');
pp = get(gcf,'PaperPosition');
set(gcf,'Units',pu,'Position',pp)

PR = PS(1)/PS(2);

order = [1, 5,2,4,6,3];
order2 = [5,6,3,4,1,2];

ax1 = axes('position',pos1);
hold on

    ylabel('test log-likelihood /nats','FontName',FontName,'FontSize',FSlg)
  set(gca,'FontName',FontName,'FontSize',FSsm)
set(gca,'TickDir','out')
    xlabel('datasets','FontName',FontName,'FontSize',FSlg)
set(gca,'xtick',[1],'xticklabel',datanames(order(1)))
    
h1=errorbar(ADF(order(1)),ADFEB(order(1))/2,'or','linewidth',lw1,'markersize',ms,'markerfacecolor',[1,0,0]);
h2=errorbar(EP(order(1)),EPEB(order(1))/2,'ok','linewidth',lw2+2,'markersize',ms+2,'markerfacecolor',[0,0,0]);
h3 = errorbar(SEP(order(1)),SEPEB(order(1))/2,'om','linewidth',lw3,'markersize',ms,'markerfacecolor',[1,0,1]);

legend([h1,h2,h3],'ADF','EP','SEP')

ax2 = axes('position',pos2);
hold on

 %   ylabel('test log-likelihood /nats','FontName',FontName,'FontSize',FSlg)
  set(gca,'FontName',FontName,'FontSize',FSsm)
set(gca,'TickDir','out')
 %   xlabel('datasets','FontName',FontName,'FontSize',FSlg)
set(gca,'xtick',[1],'xticklabel',datanames(order(2)))
    
h1=errorbar(ADF(order(2)),ADFEB(order(2))/2,'or','linewidth',lw1,'markersize',ms,'markerfacecolor',[1,0,0]);
h2=errorbar(EP(order(2)),EPEB(order(2))/2,'ok','linewidth',lw2,'markersize',ms,'markerfacecolor',[0,0,0]);
h3 = errorbar(SEP(order(2)),SEPEB(order(2))/2,'om','linewidth',lw3,'markersize',ms,'markerfacecolor',[1,0,1]);


ax3 = axes('position',pos3);
hold on

 %   ylabel('test log-likelihood /nats','FontName',FontName,'FontSize',FSlg)
  set(gca,'FontName',FontName,'FontSize',FSsm)
set(gca,'TickDir','out')
 %   xlabel('datasets','FontName',FontName,'FontSize',FSlg)
set(gca,'xtick',[1],'xticklabel',datanames(order(3)))
    
h1=errorbar(ADF(order(3)),ADFEB(order(3))/2,'or','linewidth',lw1,'markersize',ms,'markerfacecolor',[1,0,0]);
h2=errorbar(EP(order(3)),EPEB(order(3))/2,'ok','linewidth',lw2,'markersize',ms,'markerfacecolor',[0,0,0]);
h3 = errorbar(SEP(order(3)),SEPEB(order(3))/2,'om','linewidth',lw3,'markersize',ms,'markerfacecolor',[1,0,1]);



ax4 = axes('position',pos4);
hold on

 %   ylabel('test log-likelihood /nats','FontName',FontName,'FontSize',FSlg)
  set(gca,'FontName',FontName,'FontSize',FSsm)
set(gca,'TickDir','out')
 %   xlabel('datasets','FontName',FontName,'FontSize',FSlg)
set(gca,'xtick',[1],'xticklabel',datanames(order(4)))
    
h1=errorbar(ADF(order(4)),ADFEB(order(4))/2,'or','linewidth',lw1,'markersize',ms,'markerfacecolor',[1,0,0]);
h2=errorbar(EP(order(4)),EPEB(order(4))/2,'ok','linewidth',lw2,'markersize',ms,'markerfacecolor',[0,0,0]);
h3 = errorbar(SEP(order(4)),SEPEB(order(4))/2,'om','linewidth',lw3,'markersize',ms,'markerfacecolor',[1,0,1]);




ax5 = axes('position',pos5);
hold on

 %   ylabel('test log-likelihood /nats','FontName',FontName,'FontSize',FSlg)
  set(gca,'FontName',FontName,'FontSize',FSsm)
set(gca,'TickDir','out')
 %   xlabel('datasets','FontName',FontName,'FontSize',FSlg)
set(gca,'xtick',[1],'xticklabel',datanames(order(5)))
    
h1=errorbar(ADF(order(5)),ADFEB(order(5))/2,'or','linewidth',lw1,'markersize',ms,'markerfacecolor',[1,0,0]);
h2=errorbar(EP(order(5)),EPEB(order(5))/2,'ok','linewidth',lw2,'markersize',ms,'markerfacecolor',[0,0,0]);
h3 = errorbar(SEP(order(5)),SEPEB(order(5))/2,'om','linewidth',lw3,'markersize',ms,'markerfacecolor',[1,0,1]);



ax6 = axes('position',pos6);
hold on

 %   ylabel('test log-likelihood /nats','FontName',FontName,'FontSize',FSlg)
  set(gca,'FontName',FontName,'FontSize',FSsm)
set(gca,'TickDir','out')
 %   xlabel('datasets','FontName',FontName,'FontSize',FSlg)
set(gca,'xtick',[1],'xticklabel',datanames(order(6)))
    
h1=errorbar(ADF(order(6)),ADFEB(order(6))/2,'or','linewidth',lw1,'markersize',ms,'markerfacecolor',[1,0,0]);
h2=errorbar(EP(order(6)),EPEB(order(6))/2,'ok','linewidth',lw2,'markersize',ms,'markerfacecolor',[0,0,0]);
h3 = errorbar(SEP(order(6)),SEPEB(order(6))/2,'om','linewidth',lw3,'markersize',ms,'markerfacecolor',[1,0,1]);

%legend([h1,h2,h3],'ADF','EP','SEP')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 figure2 = figure;
    
    PP = [0,0,20.00,12.40]*0.7; %*** paper position in centimeters
    PS = PP(end-1:end); % paper size in centimeters

set(figure2,'paperpositionmode','manual','paperposition', ...
        PP,'papersize',PS, 'paperunits','centimeters');

% So the figure is the same size on the screen as when it is printed:
pu = get(gcf,'PaperUnits');
pp = get(gcf,'PaperPosition');
set(gcf,'Units',pu,'Position',pp)

PR = PS(1)/PS(2);

ax1 = axes('position',pos1);
hold on

    ylabel('test log-likelihood /nats','FontName',FontName,'FontSize',FSlg)
  set(gca,'FontName',FontName,'FontSize',FSsm)
set(gca,'TickDir','out')
    xlabel('datasets','FontName',FontName,'FontSize',FSlg)
set(gca,'xtick',[1],'xticklabel',datanames(order2(1)))
    
title(tits{order2(1)},'FontName',FontName,'FontSize',FSlg-4)

h1=errorbar(ADF2(order2(1)),ADFEB2(order2(1))/2,'or','linewidth',2,'markersize',ms,'markerfacecolor',[1,0,0]);
h2=errorbar(EP2(order2(1)),EPEB2(order2(1))/2,'ok','linewidth',2,'markersize',ms,'markerfacecolor',[0,0,0]);
h3 = errorbar(SEP2(order2(1)),SEPEB2(order2(1))/2,'om','linewidth',2,'markersize',ms,'markerfacecolor',[1,0,1]);

legend([h1,h2,h3],'ADF','EP','SEP')

ax2 = axes('position',pos2);
hold on
title(tits{order2(2)},'FontName',FontName,'FontSize',FSlg-4)

 %   ylabel('test log-likelihood /nats','FontName',FontName,'FontSize',FSlg)
  set(gca,'FontName',FontName,'FontSize',FSsm)
set(gca,'TickDir','out')
 %   xlabel('datasets','FontName',FontName,'FontSize',FSlg)
set(gca,'xtick',[1],'xticklabel',datanames2(order2(2)))
    
h1=errorbar(ADF2(order2(2)),ADFEB2(order2(2))/2,'or','linewidth',2,'markersize',ms,'markerfacecolor',[1,0,0]);
h2=errorbar(EP2(order2(2)),EPEB2(order2(2))/2,'ok','linewidth',2,'markersize',ms,'markerfacecolor',[0,0,0]);
h3 = errorbar(SEP2(order2(2)),SEPEB2(order2(2))/2,'om','linewidth',2,'markersize',ms,'markerfacecolor',[1,0,1]);


ax3 = axes('position',pos3);
hold on
title(tits{order2(3)},'FontName',FontName,'FontSize',FSlg-4)

 %   ylabel('test log-likelihood /nats','FontName',FontName,'FontSize',FSlg)
  set(gca,'FontName',FontName,'FontSize',FSsm)
set(gca,'TickDir','out')
 %   xlabel('datasets','FontName',FontName,'FontSize',FSlg)
set(gca,'xtick',[1],'xticklabel',datanames2(order2(3)))
    
h1=errorbar(ADF2(order2(3)),ADFEB2(order2(3))/2,'or','linewidth',2,'markersize',ms,'markerfacecolor',[1,0,0]);
h2=errorbar(EP2(order2(3)),EPEB2(order2(3))/2,'ok','linewidth',2,'markersize',ms,'markerfacecolor',[0,0,0]);
h3 = errorbar(SEP2(order2(3)),SEPEB2(order2(3))/2,'om','linewidth',2,'markersize',ms,'markerfacecolor',[1,0,1]);



ax4 = axes('position',pos4);
hold on
title(tits{order2(4)},'FontName',FontName,'FontSize',FSlg-4)

 %   ylabel('test log-likelihood /nats','FontName',FontName,'FontSize',FSlg)
  set(gca,'FontName',FontName,'FontSize',FSsm)
set(gca,'TickDir','out')
 %   xlabel('datasets','FontName',FontName,'FontSize',FSlg)
set(gca,'xtick',[1],'xticklabel',datanames2(order2(4)))
    
h1=errorbar(ADF2(order2(4)),ADFEB2(order2(4))/2,'or','linewidth',2,'markersize',ms,'markerfacecolor',[1,0,0]);
h2=errorbar(EP2(order2(4)),EPEB2(order2(4))/2,'ok','linewidth',2,'markersize',ms,'markerfacecolor',[0,0,0]);
h3 = errorbar(SEP2(order2(4)),SEPEB2(order2(4))/2,'om','linewidth',2,'markersize',ms,'markerfacecolor',[1,0,1]);




ax5 = axes('position',pos5);
hold on
title(tits{order2(5)},'FontName',FontName,'FontSize',FSlg-4)

 %   ylabel('test log-likelihood /nats','FontName',FontName,'FontSize',FSlg)
  set(gca,'FontName',FontName,'FontSize',FSsm)
set(gca,'TickDir','out')
 %   xlabel('datasets','FontName',FontName,'FontSize',FSlg)
set(gca,'xtick',[1],'xticklabel',datanames2(order2(5)))
    
h1=errorbar(ADF2(order2(5)),ADFEB2(order2(5))/2,'or','linewidth',2,'markersize',ms,'markerfacecolor',[1,0,0]);
h2=errorbar(EP2(order2(5)),EPEB2(order2(5))/2,'ok','linewidth',2,'markersize',ms,'markerfacecolor',[0,0,0]);
h3 = errorbar(SEP2(order2(5)),SEPEB2(order2(5))/2,'om','linewidth',2,'markersize',ms,'markerfacecolor',[1,0,1]);



ax6 = axes('position',pos6);
hold on
title(tits{order2(6)},'FontName',FontName,'FontSize',FSlg-4)

 %   ylabel('test log-likelihood /nats','FontName',FontName,'FontSize',FSlg)
  set(gca,'FontName',FontName,'FontSize',FSsm)
set(gca,'TickDir','out')
 %   xlabel('datasets','FontName',FontName,'FontSize',FSlg)
set(gca,'xtick',[1],'xticklabel',datanames2(order2(6)))
    
h1=errorbar(ADF2(order2(6)),ADFEB2(order2(6))/2,'or','linewidth',2,'markersize',ms,'markerfacecolor',[1,0,0]);
h2=errorbar(EP2(order2(6)),EPEB2(order2(6))/2,'ok','linewidth',2,'markersize',ms,'markerfacecolor',[0,0,0]);
h3 = errorbar(SEP2(order2(6)),SEPEB2(order2(6))/2,'om','linewidth',2,'markersize',ms,'markerfacecolor',[1,0,1]);

