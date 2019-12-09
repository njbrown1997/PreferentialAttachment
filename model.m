%% Age-Weight Driver
% This file explores the effects of age-weighting on Barabasi-Albert style network 
% generation. Age functions are considered as follows:

% PART 1: linear, power-law:  f(age) = 1/age^alpha
% PART 2: threshold:          f(age) = (1 if age<threshold, 0 if age>threshold)
% PART 3: normal distribution: f(age) = N(mu,lambda)

%for all networks, the following will be found
% (i) power-law coeff of degree distribution
% (ii) mean shortest path
% (iii) assortativity coefficient
%% PART 1: none/linear/power-law age-weighting
% for a given set of parameters, each of the above will be explored. A distribution 
% will be generated for each parametrization 

%set parameters to be constant across all of part 1:
N=200;
M=N; 
nreps=20;


%% 1A: linear
%generates linear_nets, containing a set of networks generated with
% f(age) linearly decreasing w/ age.

%set parameters
%nreps=20; %SET TO 20+ FOR FULL RUN

tic
linear_nets={};
for net_idx=1:nreps
      
    %define current f(age)
    age_fcn=@(x) 1-(x-1)/(N-1);

    %run network generation with relevant input function
    A=build_network(N,age_fcn);

    linear_nets{net_idx}=A;
% % %     G = graph(A);
% % %     figure
% % %     plot(G,'layout','force');
% % %     title('linear')
end
toc


%% 1B: constant + power-law
%generates powerlaw_nets, containing a set of networks generated with
% f(age) proportional to age^-a, for various a values (0=constant age).


%set parameters
%nreps=20; %SET TO 20+ FOR FULL
alphas=[0 .5 1 2 4]; 


powerlaw_nets={};
tic
%loop through each alpha
for alpha_idx=1:length(alphas)
    alpha=alphas(alpha_idx);
    
    for net_idx=1:nreps
    
        %define current f(age)
        age_fcn=@(x) x^(-alpha);
        
        %run network generation with relevant input function
        A=build_network(N,age_fcn);

        powerlaw_nets{alpha_idx,net_idx}=A;
% % %         if alpha==0
% % %             G = graph(A);
% % %             figure
% % %             plot(G,'layout','force');
% % %             title(['Power-law: \alpha=-' num2str(alpha)])
% % %         end
    end
end
toc


%% 1C: threshold
%generates threshold_nets, containing a set of networks generated with
% f(age) =1, or =0 after a threshold.

%set parameters
%nreps=20; %SET TO 20+ FOR FULL
thresholds=[5 25 50 100]; %SET TO 0:10:200 FOR FULL RUN


threshold_nets={};

%loop through each alpha
tic
for thresh_idx=1:length(thresholds)
    thresh=thresholds(thresh_idx);
    
    for net_idx=1:nreps

        %define current f(age)
        age_fcn=@(x) x<=thresh;
        
        %run network generation with relevant input function
        A=build_network(N,age_fcn);

        threshold_nets{thresh_idx,net_idx}=A;
% % %         G = graph(A);
% % %         figure
% % %         plot(G,'layout','force');
% % %         title(['Power-law: t_{threshold}=' num2str(thresh)])
    end
end
toc



%% 1D: NORMAL DISTRIBUTION AGING
%generates:
%   (a) norm_midpoint_nets, containing a set of networks generated with f(age)
%        proportional to a normal distribution, sweeping through midpoints
%   (b) norm_width_nets, fixed midpoint but varying widths 

%set parameters
%nreps=20; %SET TO 20+

%for varying midpoint
midpoint_sweep=0:25:100; 
stdev_fixed=20; 

%for varying width
midpoint_fixed=50;
stdev_sweep=[5 10 25 50 100];

tic
% (a) LOOP THROUGH MIDPOINT
norm_midpoint_nets={};
for midpoint_idx=1:length(midpoint_sweep)
    midpoint=midpoint_sweep(midpoint_idx);
    
    for net_idx=1:nreps
        
        %define current f(age)
        age_fcn=@(x) normpdf(x,midpoint,stdev_fixed);
        
        %run network generation with relevant input function
        A=build_network(N,age_fcn);

        norm_midpoint_nets{midpoint_idx,net_idx}=A;
% % %         G = graph(A);
% % %         figure
% % %         plot(G,'layout','force');
% % %         title(['Curve: midpoint=' num2str(midpoint)])
    end
end
toc


tic
% (b) LOOP THROUGH WIDTHS
norm_width_nets={};
for stdev_idx=1:length(stdev_sweep)
    stdev=stdev_sweep(stdev_idx);
    
    for net_idx=1:nreps
        
        %define current f(age)
        age_fcn=@(x) normpdf(x,midpoint_fixed,stdev);
        
        %run network generation with relevant input function
        A=build_network(N,age_fcn);

        norm_width_nets{stdev_idx,net_idx}=A;
% % %         G = graph(A);
% % %         figure
% % %         plot(G,'layout','force');
% % %         title(['Curve: stdev=' num2str(stdev)])
    end
end
toc





%% PART 4: NETWORK STATISTICS

%% 4A: compute CCDF and assortativity for each network set

%...for linear nets
linear_assort=zeros(size(linear_nets,2),1);
linear_CCDF=zeros([size(linear_nets,2) N]);
for i=1:size(linear_nets,2)  
    linear_assort(i)=assortativity(linear_nets{i},0);
    [f, x] = ccdf(sum(linear_nets{i},2));
    linear_CCDF(i,:)=f;
end

%...for powerlaw nets
constant_assort=[];
powerlaw_assort=[];
powerlaw_assort_key=[];
powerlaw_CCDF=zeros([size(powerlaw_nets) N]);
for i=1:size(powerlaw_nets,2)
    for alpha_i=1:size(powerlaw_nets,1)
        if alphas(alpha_i)==0
            constant_assort=[constant_assort assortativity(powerlaw_nets{alpha_i,i},0)];
        else
            powerlaw_assort=[powerlaw_assort assortativity(powerlaw_nets{alpha_i,i},0)];
            powerlaw_assort_key=[powerlaw_assort_key alphas(alpha_i)];
        end
        [f, x] = ccdf(sum(powerlaw_nets{alpha_i,i},2));
        powerlaw_CCDF(alpha_i,i,:)=f;
    end
end

%...for threshold nets
threshold_assort=[];
threshold_assort_key=[];
threshold_CCDF=zeros([size(threshold_nets) N]);
for i=1:size(threshold_nets,2)
    for t_i=1:size(threshold_nets,1)
        threshold_assort=[threshold_assort assortativity(threshold_nets{t_i,i},0)];
        threshold_assort_key=[threshold_assort_key alphas(t_i)];
        [f, x] = ccdf(sum(threshold_nets{t_i,i},2));
        threshold_CCDF(t_i,i,:)=f;
    end
end

%...for midpoint nets
midpoint_assort=[];
midpoint_assort_key=[];
midpoint_CCDF=zeros([size(norm_midpoint_nets) N]);
for i=1:size(norm_midpoint_nets,2)
    for mid_i=1:size(norm_midpoint_nets,1)
        midpoint_assort=[midpoint_assort assortativity(norm_midpoint_nets{mid_i,i},0)];
        midpoint_assort_key=[midpoint_assort_key alphas(mid_i)];
        [f, x] = ccdf(sum(norm_midpoint_nets{mid_i,i},2));
        midpoint_CCDF(mid_i,i,:)=f;
    end
end

%...for width nets
width_assort=[];
width_assort_key=[];
width_CCDF=zeros([size(norm_width_nets) N]);
for i=1:size(norm_width_nets,2)
    for std_i=1:size(norm_width_nets,1)
        width_assort=[width_assort assortativity(norm_width_nets{std_i,i},0)];
        width_assort_key=[width_assort_key alphas(std_i)];
        [f, x] = ccdf(sum(norm_width_nets{std_i,i},2));
        width_CCDF(std_i,i,:)=f;
    end
end

%% 4B: Figure 1: CCDFs
figure

%%% FIG 1a: constant age weight vs linear decrease %%%
subplot(2,3,1)
hold on
%plot constant f(age)
plot(1:N,squeeze(mean(powerlaw_CCDF(find(alphas==0),:,:),2)),...
        'color','k','linewidth',4,'DisplayName','constant','linestyle',':')
        set(gca,'yscale','log')
        set(gca,'xscale','log') 
%plot linear
plot(1:N,squeeze(mean(linear_CCDF,1)),...
        'color','r','linewidth',4,'DisplayName','linear','linestyle','-')
        set(gca,'yscale','log')
        set(gca,'xscale','log') 
hold off
legend()
title('f(age) \propto x','fontsize',18)


%%% FIG 1b: constant age weight vs various powerlaw exponents %%%
subplot(2,3,2)
hold on
%plot constant
plot(1:N,squeeze(mean(powerlaw_CCDF(find(alphas==0),:,:),2)),...
        'color','k','linewidth',4,'DisplayName','constant','linestyle',':')
        set(gca,'yscale','log')
        set(gca,'xscale','log') 
%plot power-law  
colors={[1 0 0],[0 0 1],[0 .8 .3],[0 0 0],[.4 0 0]};
for alpha_i=2:size(powerlaw_nets,1)
    color=colors{alpha_i-1};
    plot(1:N,squeeze(mean(powerlaw_CCDF(alpha_i,:,:),2)),...
        'color',[color .8],'linewidth',3,'DisplayName',['\alpha=' num2str(alphas(alpha_i),'%0.1g')])
    set(gca,'yscale','log')
    set(gca,'xscale','log') 
end
legend()
hold off
title('f(age) \propto x^{-\alpha}','fontsize',18)
%


%%% FIG 1c: constant age weight vs various thresholds %%%
subplot(2,3,3)
hold on
%plot constant
plot(1:N,squeeze(mean(powerlaw_CCDF(find(alphas==0),:,:),2)),...
        'color','k','linewidth',4,'DisplayName','constant','linestyle',':')
        set(gca,'yscale','log')
        set(gca,'xscale','log') 
%plot power-law  
colors={[.4 0 0],[0 0 0],[0 .8 .3],[0 0 1],[1 0 0],};
for t_i=1:size(threshold_nets,1)
    color=colors{t_i+1};
    plot(1:N,squeeze(mean(threshold_CCDF(t_i,:,:),2)),...
        'color',[color .8],'linewidth',3,'DisplayName',['threshold=' num2str(thresholds(t_i),'%2.1d')])
    set(gca,'yscale','log')
    set(gca,'xscale','log') 
end
legend()
hold off
title('f(age) = 1 if x\leq threshold, 0 otherwise','fontsize',18)

%%% FIG 1d: constant age weight vs various midpoints %%%
subplot(2,3,4)
hold on
%plot constant
plot(1:N,squeeze(mean(powerlaw_CCDF(find(alphas==0),:,:),2)),...
        'color','k','linewidth',4,'DisplayName','constant','linestyle',':')
        set(gca,'yscale','log')
        set(gca,'xscale','log') 
%plot power-law  
colors={[0 0 0],[0 .8 .3],[0 0 1],[1 0 0],[.4 0 0],[.4 .4 0]};
for mid_i=1:size(norm_midpoint_nets,1)
    color=colors{mid_i};
    plot(1:N,squeeze(mean(midpoint_CCDF(mid_i,:,:),2)),...
        'color',[color .8],'linewidth',3,'DisplayName',['\mu=' num2str(midpoint_sweep(mid_i),'%2.1d')])
    set(gca,'yscale','log')
    set(gca,'xscale','log') 
end
legend()
hold off
title({['f(age) \propto N(\mu,\lambda=' num2str(stdev_fixed) ')'],'normal dist., variable center'},'fontsize',18)



%%% FIG 1e: constant age weight vs various widths %%%
subplot(2,3,5)
hold on
%plot constant
plot(1:N,squeeze(mean(powerlaw_CCDF(find(alphas==0),:,:),2)),...
        'color','k','linewidth',4,'DisplayName','constant','linestyle',':')
        set(gca,'yscale','log')
        set(gca,'xscale','log') 
%plot power-law  
colors={[0 0 0],[0 .8 .3],[0 0 1],[1 0 0],[.4 0 0],[.4 .4 0]};
for std_i=1:size(norm_width_nets,1)
    color=colors{std_i};
    plot(1:N,squeeze(mean(width_CCDF(std_i,:,:),2)),...
        'color',[color .8],'linewidth',3,'DisplayName',['\lambda=' num2str(stdev_sweep(std_i),'%2.1d')])
    set(gca,'yscale','log')
    set(gca,'xscale','log') 
end
legend()
hold off
title({['f(age) \propto N(\mu=' num2str(midpoint_fixed) ',\lambda)'],'normal dist., variable width'},'fontsize',18)




%% 4B: Figure 2: assortativity by function/parameter
figure

%linear and constant 
subplot(1,5,1)
hold on
boxplot([constant_assort' linear_assort],...
    'Labels',{'constant','linear decrease'})
xlabel('f(age)','fontsize',18)
ylabel('r_{assortativity}','fontsize',18)
ylim([-1 0])
title('Linear','fontsize',18)
pos = get(gca, 'Position');
set(gca, 'Position', [pos(1)-.055 .15 pos(3)+.075 .76])
plot([0 10],[1 1].*mean(constant_assort),'linestyle',':','linewidth',2)
hold off

%power-law
subplot(1,5,2)
hold on
boxplot(powerlaw_assort,powerlaw_assort_key)
title('Power-law','fontsize',18)
xlabel('\alpha','fontsize',18)
ylim([-1 0])
pos = get(gca, 'Position');
set(gca, 'Position', [pos(1) .15 pos(3) .76])
plot([0 10],[1 1].*mean(constant_assort),'linestyle',':','linewidth',2)
hold off

%threshold
subplot(1,5,3)
hold on
boxplot(threshold_assort,threshold_assort_key)
title('Threshold','fontsize',18)
xlabel('t_{threshold}','fontsize',18)
ylim([-1 0])
pos = get(gca, 'Position');
set(gca, 'Position', [pos(1) .15 pos(3) .76])
plot([0 10],[1 1].*mean(constant_assort),'linestyle',':','linewidth',2)
hold off

%normal midpoint
subplot(1,5,4)
hold on
boxplot(midpoint_assort,midpoint_assort_key)
title('N(\mu,\lambda=20)','fontsize',18)
xlabel('\mu','fontsize',18)
ylim([-1 0])
pos = get(gca, 'Position');
set(gca, 'Position', [pos(1) .15 pos(3) .76])
plot([0 10],[1 1].*mean(constant_assort),'linestyle',':','linewidth',2)
hold off

%normal width
subplot(1,5,5)
hold on
boxplot(width_assort,width_assort_key)
title('N(\mu=50,\lambda)','fontsize',18)
xlabel('\lambda','fontsize',18)
ylim([-1 0])
pos = get(gca, 'Position');
set(gca, 'Position', [pos(1) .15 pos(3) .76])
plot([0 10],[1 1].*mean(constant_assort),'linestyle',':','linewidth',2)
hold off


%% up to here as of Monday Dec 9
%TODO NEXT: 
% a. find average shortest paths for each network, 
% b. diameters
% c. and possibly more.

% THEN:
% Plot the above a,b against the parameters of each model (alpha,
% lambda, etc)








%%
%Power-law statistics.
for k = 1:nreps
    for i=1:length(alphas)
        %Degree Distribution.
        G = graph(powerlaw_nets{i,1});
        D = degree(G);
        degDist = histcounts(D,max(D));
        X = [];
        Y = [];
        for j = 1:size(degDist,2)
            if not(degDist(j) == 0)
                X = [X,j];
                Y = [Y,degDist(j)];
            end
        end
        p = polyfit(log(X),log(Y),1);
        plExp = p(1);
        
        %Average Shortest Path.
        sPaths = distances(G);
        avgSPath = mean2(sPaths);
        
        %Assortativity Coefficient
        assort = assortativity(powerlaw_nets{i,k},0);
        
        stats = [plExp,avgSPath,assort];
        powerlaw_stats{i,k} = stats;
    end
end

%Threshold statistics.
for k = 1:nreps
    for i=1:length(thresholds)
        %Degree Distribution.
        G = graph(threshold_nets{i,1});
        D = degree(G);
        degDist = histcounts(D,max(D));
        X = [];
        Y = [];
        for j = 1:size(degDist,2)
            if not(degDist(j) == 0)
                X = [X,j];
                Y = [Y,degDist(j)];
            end
        end
        p = polyfit(log(X),log(Y),1);
        plExp = p(1);
        
        %Average Shortest Path.
        sPaths = distances(G);
        avgSPath = mean2(sPaths);
        
        %Assortativity Coefficient
        assort = assortativity(threshold_nets{i,k},0);
        
        stats = [plExp,avgSPath,assort];
        threshold_stats{i,k} = stats;
    end
end

%%
%scratch code for finding powerlaw slope
%fixed_xmin=6;
%x=sum(powerlaw_nets{2,1},1);
%[alpha, xmin, ~]=plfit(x,'finite') %'xmin',fixed_xmin,
%[p,gof]=plpva(x,xmin,'reps',20,'silent')
%plplot(x,xmin,alpha)

%%
% FOR ALL NETWORKS:
%    -find power-law exponent of degree distribution's CCDF
%    -find average shortest path between all nodes (using distances() func)
%    -find coefficient of assortativity

% plot these variables as a function of the parameters:
%    -for power-law f(age), x axis is alpha
%    -for threshold f(age), x axis is threshold
%    -for normalized poisson f(age), x axis is lambda








% % % % % % % %% PLOT NETS SO FAR
% % % % % % % 
% % % % % % % %plot powerlaw
% % % % % % % figure
% % % % % % % len_a=length(alphas);
% % % % % % % for i=1:len_a
% % % % % % %     subplot(1,len_a,i)
% % % % % % %     plot(graph(powerlaw_nets{i,1}),'layout','force');
% % % % % % %     title(['Power-law: \alpha=-' num2str(alphas(i))])
% % % % % % % end
% % % % % % % 
% % % % % % % %plot threshold
% % % % % % % figure
% % % % % % % len_th=length(thresholds);
% % % % % % % for i=1:len_th
% % % % % % %     subplot(1,len_th,i)
% % % % % % %     plot(graph(threshold_nets{i,1}),'layout','force');
% % % % % % %     title(['Threshold: t=' num2str(thresholds(i))])
% % % % % % % end
% % % % % % % %%
% % % % % % % adj2gephilab('toGephi',A);
% % % % % % % %%
% % % % % % % function EdgeL=adj2gephilab(filename,ADJ,parameters)
% % % % % % %     % Convert ana adjacency matrix of a graph to 2 spreadhseets csv files
% % % % % % %     % one for the edge table and the other for node table.
% % % % % % %     % The files _node.csv and _edge.csv have to be open 
% % % % % % %     %in Gephi via Data Laboratory.
% % % % % % %     % INPUTS:
% % % % % % %     %           filename: string for the prefix name of the two files .csv
% % % % % % %     %           ADJ: the adjacency matrix
% % % % % % %     %           parameters: vector as properties of the node to use as
% % % % % % %     %                              attributes of them.
% % % % % % %     % OUTPUTS:
% % % % % % %     %            two csv spreadsheet files: 
% % % % % % %     %                       filename_node.csv
% % % % % % %     %                       filename_edge.csv
% % % % % % %     %             EdgeL = it returns the edge list corresponing to the
% % % % % % %     %             adjacency matrix. (it can be saved to be open in Gephi too)
% % % % % % %     %
% % % % % % %     % The two files must be open in Gephi via the Data Laboratory
% % % % % % %     nodecsv=[filename,'_node.csv'];
% % % % % % %     edgecsv=[filename,'_edge.csv'];
% % % % % % %     n=size(ADJ,1); % square adjacency matrix
% % % % % % %     if nargin<3 
% % % % % % %         parameters=ones(n,1);% all nodes have the same attributes
% % % % % % %     end
% % % % % % %     ps=parameters;
% % % % % % %     % Node Table:
% % % % % % %     % header for node csv
% % % % % % %     fidN = fopen(nodecsv,'w','native','UTF-8'); % best format for gephi
% % % % % % %     fprintf(fidN,'%s\n','Id;Label;Attribute');
% % % % % % %     % 
% % % % % % %     for i=2:n+1,
% % % % % % %         fprintf(fidN,'%s\n',[ num2str(i-1) ';"Node ' num2str(i-1) '"'...
% % % % % % %             ';' num2str(ps(i-1))]);
% % % % % % %     end
% % % % % % %     fclose(fidN);
% % % % % % %     % Edge Table
% % % % % % %     EdgeL=conv_EdgeList(ADJ);
% % % % % % %     S=EdgeL(:,1); % sources
% % % % % % %     T=EdgeL(:,2); % targets
% % % % % % %     W = EdgeL(:,3); % weights
% % % % % % %     fidE = fopen(edgecsv,'w','native','UTF-8');
% % % % % % %     % header for edge csv
% % % % % % %     fprintf(fidE,'%s\n','Source;Target;Label;Weight');
% % % % % % %     for i=2:length(S)+1,
% % % % % % %           fprintf(fidN,'%s\n',[ num2str(S(i-1)) ';' num2str(T(i-1)) ';'...
% % % % % % %               '"Edge from ' num2str(S(i-1)) ' to ' num2str(T(i-1)) '"' ';'...
% % % % % % %               num2str(W(i-1))]);
% % % % % % %     end
% % % % % % %      fclose(fidE);
% % % % % % % end
% % % % % % % %% Aux function
% % % % % % % 
% % % % % % % function EdgeL=conv_EdgeList(adj)
% % % % % % %     % convert adj matrix to edge list
% % % % % % %     n=size(adj,1); % number of nodes
% % % % % % %     edges=find(adj>0); % indices of all edges
% % % % % % %     n_e=length(edges);
% % % % % % %     EdgeL=zeros(n_e,3);
% % % % % % %     for e=1:n_e
% % % % % % %       [i,j]=ind2sub([n,n],edges(e)); % node indices of edge e  
% % % % % % %       EdgeL(e,:)=[i j adj(i,j)];
% % % % % % %     end
% % % % % % % end
