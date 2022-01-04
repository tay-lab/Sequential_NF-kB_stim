%%load data
clc
clear
load('scmat_sequentialstim.mat')
load('featmat_sequentialstim.mat')

maxidx_all = [1, 8, 15, 22];
aucidx_all = [7, 14, 21, 28];
maxtwoidx_all = [4, 11, 18, 15];
timeidx_all = [2, 9, 16, 23];
widthidx_all = [3, 10, 17, 24];

%%  F1D: representative single cell traces
rng(1)
scmat = scmatcomb_norm;
traces = scmat(scmat(:, 168) == 4 & scmat(:, 169) == 5 & scmat(:, 170) == 6 & scmat(:, 171) == 7, 1:83);
traces = datasample(traces, 10, 1, 'Replace', false);
ylim1 = [-0.2 6.5];
xlim = [0 485];
plot(0:6:492, traces(:, 1:83), 'linewidth', 1.5)
xline(120, '--', 'linewidth', 2.5)
xline(240, '--', 'linewidth', 2.5)
xline(360, '--', 'linewidth', 2.5)
set(gca,'Xlim', xlim, 'Ylim', ylim1);
ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
    
%% Fig 1D: peakheight heatmap by order
rng(2)

scmat = scmatcomb_norm(any(isnan(scmatcomb_norm),2)==0,:);
max_all = featmat(:, maxidx_all);
maxmat = [];
maxall_mat = [];
traceall_mat = [];

treatmat = [];
treatall_mat = [];
condlabs = {'T I L P', 'T I P L', 'T L I P', 'T L P I', 'T P I L', 'T P L I', 'I T L P', 'I T P L', 'I L T P', 'I L P T', ...
    'I P T L', 'I P L T', 'L T I P', 'L T P I', 'L I T P', 'L I P T', 'L P T I', 'L P I T', 'P T L I', 'P T I L', 'P I T L',...
    'P I L T', 'P L T I', 'P L I T'};
ylabs = {'S1', 'S2', 'S3', 'S4'};
conditions = unique(scmat(:, 167));
for aa = [1:24, 26:49, 51:74]
    maxmat = [maxmat; [mean(max_all(scmat(:, 167) == conditions(aa), :)), conditions(aa)] ];
    if aa ~=63:68
        traceall_mat = [traceall_mat; sortrows(datasample( scmat(scmat(:, 167) == conditions(aa), 1:83), 50, 'Replace', false), 5, 'descend') ]; 
    else
        traceall_mat = [traceall_mat; sortrows(datasample( scmat(scmat(:, 167) == conditions(aa), 1:83), 50, 'Replace', false), 10, 'descend') ]; 

    end
    treatmat = [treatmat; scmat(find(scmat(:, 167) == conditions(aa), 1, 'first'), 168:172)];
    treatall_mat = [treatall_mat; scmat(find(scmat(:, 167) == conditions(aa), 50), 168:172)];
    
end
treatall_mat = [repmat(treatall_mat(:,1), [1, 20]), repmat(treatall_mat(:,2), [1, 20]), repmat(treatall_mat(:,3), [1, 20]),...
    repmat(treatall_mat(:,4), [1, 20]),  ones([length(treatall_mat(:,1)), 3])*99,  ];

normmat = [];
normall_mat = [];
traceall_mat = traceall_mat';
maxmat_trim = maxmat(:, 1:4)';
foo = [];

redwhiteblue = [0 0 1];
for aa = 1:99
    temp = redwhiteblue(end, :);
    redwhiteblue = [redwhiteblue; [temp(1) + 0.01, temp(2) + 0.01, 1] ];
end

for aa = 1:98
    temp = redwhiteblue(end, :);
    redwhiteblue = [redwhiteblue; [1, temp(3) - 0.01, temp(3) - 0.01] ];
end
normmat = treatmat(:, 1:4)';
normall_mat = treatall_mat';

for aa = 3 %pick dose 1 = high, 2 = mid, 3 = low
    for bb = 0:3
        peak_l1 = mean( max_all(scmat(:, 172) == aa & scmat(:, 168) == bb + (aa-1)*4));

        normmat(find( [treatmat(:,1:4) == bb + (aa-1)*4]' ))  = maxmat_trim( [treatmat(:,1:4) == bb + (aa-1)*4]')./peak_l1 ;
        maxmat_trim( [treatmat(:,1:4) == bb + (aa-1)*4]')./peak_l1 ;
        
        normall_mat(find( [treatall_mat(:,1:end) == bb + (aa-1)*4]' )) = traceall_mat( [treatall_mat(:,1:end) == bb + (aa-1)*4]')./peak_l1;

    end
    if aa == 1
        normall_mat = normall_mat(1:80, 1:1200)';
    elseif aa == 2
        normall_mat = normall_mat(1:80, 1201:2400)';
    else
        normall_mat = normall_mat(1:80, 2401:3600)';
    end
    figure(1)
    normall_mat(normall_mat < 0) = 0;
    colormap(double( redwhiteblue) )
    imagesc(log10( normall_mat + 0.01), [-2, 0.5]) 
    set(gca,'XColor', 'none','YColor','none')
    colorbar

end

%% Fig 2D: amplitude based on dose and time


maxidx_all = [1, 8, 15, 22];
overall_feat = featmat(:, maxidx_all);

colors = hsv(5);

xnames = { 'S1', 'S2', 'S3', 'S4'};
titles = { 'High', 'Mid', 'Low'};


figure(2)
clf
plotcounter = 1;
allcells = {};
for aa = 1:3
    subplot(1, 3, plotcounter)
    hold on
    counter = 1;
    for bb = 1:4
        d1mean = mean( overall_feat(scmat(:, 172) == aa & ismember(scmat(:, 168), 1:11), 1) );
        tempmat = overall_feat(scmat(:, 172) == aa & ismember(scmat(:, 167+bb), 1:11), bb)./d1mean;
        Violin( (tempmat), bb, 'Bandwidth', 0.2,  'MedianColor', [1,1,1], 'ViolinColor', [0 0 0], 'ViolinAlpha', .3, 'BoxColor',[0 0 0], 'ShowData', false, 'Width', .32, ...
            'BoxWidth', 0.1, 'ShowMean', true); 
        [nr, ~] = size(tempmat);
        allcells{aa, bb} = tempmat;
     end
        set(gca, 'XLim', [0.25, 5], 'YLim', [0, 3])
        ax = gca;
        ax.XAxis.FontSize = 14;
        ax.YAxis.FontSize = 14;
        xticks([1 2 3 4 ])
        xticklabels(xnames)
        title(titles{plotcounter}, 'FontSize', 14)
        xlabel("Stim. order",'FontSize', 14)
        ylabel("Max Response", 'FontSize', 14)
        plotcounter = plotcounter + 1;
        hold off
end
pval_f1e = [];
for aa = 1:3
    for bb = 1:3
        p = ranksum(allcells{aa, bb}, allcells{aa, bb+1});
        pval_f1e(aa, bb) = p*9;
        
    end
    
end

esize_f1e = [];
for aa = 1:3
    for bb = 1:3
        e =  mean(allcells{aa, bb})./mean(allcells{aa, bb+1});
        esize_f1e(aa, bb) = e; 
    end
    
end


%% Fig S1A - PAM
scmatp = scmatcomb_norm(any(isnan(scmatcomb_norm),2)==0,:);
prom = scmatp( scmatp(:, 168) == 12, 1:20);
prom = max(prom');
prom = prctile(prom, 98);
% prom = 1;

load('scmat_PAMtitr.mat')
scmat = scmat_PAMtit_20200630;
interval = 6;       % Time between each frame (in min)
xlim = [0, 150];
ylim1 = [-0.2 6];
howmanytraces = 20;

chambid = unique(scmat(:, 97));



fig = figure(1); clf;
counter = 6;
for aa=2:6 %length(chambid)
    if ~ isempty(find(scmat(:, 97) == chambid(aa), 1))
        counter = counter-1;
        subplot(1,5,counter);
        %subplot(9,9,counter);
        
        temp = scmat(scmat(:, 97) == chambid(aa), 1:48);
        [nr, ~] = size(temp);
        amps = [];
        for bb=1:nr
            amp = sort( findpeaks(temp(bb, 1:20)) , 'descend');
            amps = [amps; amp(1)];
        end
        prcactive = sum(amps > prom)/length(amps)
        if nr > howmanytraces
            Sel = randperm(nr, howmanytraces);   % How many single traces to plot
        else 
            Sel = 1:length(howmanytraces);
        end
        hold on
        
        for bb=1:length(Sel)
            plot(0:6:282, temp(Sel(bb), :)-mean(temp(Sel(bb), 1:2)), '-', 'Color', [0.5 0.5 0.5, .4]);
        end
        %median
        plot(0:6:282, nanmean(temp(:,:)) - mean( mean(temp(:,1:2)) ), '-', 'Color', 'r', 'linewidth', 2);
        figtit = [];
        if aa == 2
            figtit = "0.01 ng";
        elseif aa == 3
            figtit = "0.05 ng";
        elseif aa == 4
            figtit = "0.1 ng";
        elseif aa == 5
            figtit = "0.2 ng";
        else
            figtit = "1 ng";
        end
                
        set(gca, 'Xlim', xlim, 'Ylim', ylim1);
        
    
        ax = gca;
        ax.XAxis.FontSize = 14;
        ax.YAxis.FontSize = 14;
               
        if counter == 1
            ylabel("PAM2CSK4", 'FontSize', 16)
        end
        
        if counter == 3 
            xlabel("Time (min)", 'FontSize', 16)
        end

        
        title(figtit, 'FontSize', 14)
        

        hold off
    end
end

%% Fig S1B IL-1
load('IL-1titr.mat')
xlim = [-30 120];
ylim = [-.2 5];
traces = proctrace_20190828;

counter = 0;

figure(1); clf;
for aa=[1, 3:6]

    counter = counter + 1;
    subplot(1, 5, counter); cla;
    hold on
     if aa == 3
         scmat = scmat_PAMtit_20200630;
         temp = scmat(scmat(:, 97) == chambid(1), 1:48);
         
         amps = [];
         [nr, ~] = size(temp);
         for bb=1:nr
             amp = sort( findpeaks(temp(bb, 1:20)) , 'descend');
             amps = [amps; amp(1)];
         end
         prcactive = sum(amps > prom)/length(amps)
         if nr > howmanytraces
             Sel = randperm(nr, howmanytraces);   % How many single traces to plot
         else 
             Sel = 1:length(howmanytraces);
         end
         
         for bb=1:length(Sel)
             plot(-30:6:252, temp(Sel(bb), :)-mean(temp(Sel(bb), 1:2)), '-', 'Color', [0.5 0.5 0.5, .4]);
         end
         %median
         plot(-30:6:252, nanmean(temp(:,:)) - mean( mean(temp(:,1:2)) ), '-', 'Color', 'r', 'linewidth', 2);
         
     else
     
        temp = traces{aa, 3};
        amp = [] ;
        for bb = 1:length(temp)
            amp = sort( findpeaks(temp{bb}(41:6:120, 2)) , 'descend');
            if isempty(amp)
                amps = [amps; 0];
            else
                amps = [amps; amp(1)];
            end
        end
        for bb = datasample( 1:length(temp), 20, 'Replace', false)
            plot(temp{bb,1}(:,1), temp{bb,1}(:,2)-mean( temp{bb,1}(1:5,2) ), '-', 'Color', [0.7 0.7 0.7, .5]);
        end
        plot(traces{aa, 4}(:,1), traces{aa, 4}(:,2)./traces{aa, 4}(:,3)-mean(traces{aa, 4}(1:5,2)./traces{aa, 4}(1:5,3)),...
            '-', 'Color', 'r', 'linewidth', 2);
        prcactive = [counter, sum(amps > prom)/length(amps)]

     end
     

    if counter == 1
        ylabel("IL-1", 'FontSize', 16)
    end
    
    if counter == 3
        xlabel("Time (min)", 'FontSize', 16)
    end
    
    set(gca, 'XLim', xlim, 'YLim', ylim);
    if aa == 1
        figtit = "3 ng";
    elseif aa == 3
        figtit = "0.2 ng";
    elseif aa == 4
        figtit = "0.1 ng";
    elseif aa == 5
        figtit = "0.03 ng";
    else
        figtit = "0.01 ng";
    end
    ax = gca;
    ax.XAxis.FontSize = 14;
    ax.YAxis.FontSize = 14;
    title(figtit, 'FontSize', 14)

    hold off
end

%% Fig S1C-D - TNF LPS
load('scmat_TNFLPStitr.mat')
scmat = scmat_TNFLPStit_20191015;
interval = 6;       % Time between each frame (in min)
xlim = [0, 180];
ylim1 = [-0.2 5];
howmanytraces = 20;

chambid = unique(scmat(:, 97));

titles = {"90 ng"; "30 ng"; "10 ng"; "3 ng"; "1 ng"; "400 ng"; "100 ng"; "50 ng"; "25 ng"; "12.5 ng"};


fig = figure(3); clf;
counter = 0;
for aa=[1:5, 7,9:12] %length(chambid)
    if ~ isempty(find( chambid(aa, 1) , 1))
        counter = counter+1;
        subplot(2,5,counter);
        amps = [];
        temp = scmat( ismember( scmat(:, 97),  [chambid(aa*2-1), chambid(aa*2)]), 1:48);
        for bb = 1:length(temp)
            amp = sort( findpeaks(temp(bb, 1:20 )) , 'descend');
            if isempty(amp)
                amps = [amps; 0];
            else
                amps = [amps; amp(1)];
            end
        end
        prcactive = sum(amps > prom)/length(amps)
        [nr, ~] = size(temp);
        if nr > howmanytraces
            Sel = randperm(nr, howmanytraces);   % How many single traces to plot
        else 
            Sel = 1:length(howmanytraces);
        end
        hold on
        
        for bb=1:length(Sel)
            plot(0:6:282, temp(Sel(bb), :), '-', 'Color', [0.5 0.5 0.5, .4]);
        end
        %mean
        plot(0:6:282, nanmean(temp(:,:)), '-', 'Color', 'r', 'linewidth', 2);
        figtit = titles{counter};

                
        set(gca, 'Xlim', xlim, 'Ylim', ylim1);
        
    
        ax = gca;
        ax.XAxis.FontSize = 14;
        ax.YAxis.FontSize = 14;
               
        if aa == 1
            ylabel("TNF", 'FontSize', 14)
        elseif aa == 7
            ylabel("LPS", 'FontSize', 14)
        end
        
        if ismember(counter, [3, 8]) 
            xlabel("Time (min)", 'FontSize', 14)
        end
        
        title(figtit, 'FontSize', 14)
        

        hold off
    end
end

%% Fig S1E
ylim = [0.4 1];

prct_active = [0.9042,1, 0.9889,  0.9625;
    0.9042, 0.959, 0.969, 0.8619;
    0.618, 0.55, 0.5875, 0.6889	];

figure(5); clf;
bar(prct_active)
set(gca,'YLim', ylim);
ylabel("Percent activated", 'FontSize', 14)
xticks([1 2 3 ])
xticklabels(["High", "Mid", "Low"])
ax.XAxis.FontSize = 14;
legend('TNF','IL-1', 'LPS', 'PAM')



%% Fig S2 High dose

scmat = scmatcomb_norm;
interval = 6;       % Time between each frame (in min)
xlim = [0, 480];
ylim1 = [-0.2 7];
howmanytraces = 20;

chambid = unique(scmat(:, 167));
input_names = ['T'; 'I'; 'L'; 'P'];


fig = figure(1); clf;
counter = 0;
for aa=1:24 %length(chambid)
    if ~ isempty(find(scmat(:, 167) == chambid(aa), 1))
        counter = counter+1;
        subplot(4,6,counter);
        
        temp = scmat(scmat(:, 167) == chambid(aa), 1:83);
        [nr, ~] = size(temp);
        if nr > howmanytraces
            Sel = randperm(nr, howmanytraces);   % How many single traces to plot
        else 
            Sel = 1:length(howmanytraces);
        end
        hold on
        
        for bb=1:length(Sel)
            plot(0:6:492, temp(Sel(bb), :), '-', 'Color', [0.5 0.5 0.5, .4]);
        end
        %median
        plot(0:6:492, nanmedian(temp(:,:)), '-', 'Color', 'r', 'linewidth', 2);
        xline(120, '--', 'linewidth', 1.5)
        xline(240, '--', 'linewidth', 1.5)
        xline(360, '--', 'linewidth', 1.5)
        figtit = [];
        for cc = 1:4
            inputnum = unique( scmat(scmat(:, 167) == chambid(aa), 167+cc));
            if inputnum < 4
                figtit = strcat( figtit, " ", input_names( inputnum+1, :));
            elseif inputnum < 8 
                figtit = strcat( figtit, " ", input_names( inputnum-3, :));
            elseif inputnum < 12
                figtit = strcat( figtit, " ", input_names( inputnum-7, :));
            else
                figtit = strcat( figtit, " ", "FM");
            end
        end
                
        set(gca, 'Xlim', xlim, 'Ylim', ylim1);
        
    
        ax = gca;
        ax.XAxis.FontSize = 14;
        ax.YAxis.FontSize = 14;
               
        if ismember(aa, [1, 7, 13, 19] )
            ylabel("Nuc/Cyto p65", 'FontSize', 14)
        end

        if ismember(aa, 19:24 )
            xlabel("Time (min)", 'FontSize', 14)
        end
        
        title(figtit, 'FontSize', 14)
        

        hold off
    end
end

%% Fig S3 mid dose

scmat = scmatcomb_norm;
interval = 6;       % Time between each frame (in min)
xlim = [0, 480];
ylim1 = [-0.2 5];
howmanytraces = 20;

chambid = unique(scmat(:, 167));
input_names = ['T'; 'I'; 'L'; 'P'];


fig = figure(1); clf;
counter = 0;
for aa=26:49 %length(chambid)
    if ~ isempty(find(scmat(:, 167) == chambid(aa), 1))
        counter = counter+1;
        subplot(4,6,counter);
        
        temp = scmat(scmat(:, 167) == chambid(aa), 1:83);
        [nr, ~] = size(temp);
        if nr > howmanytraces
            Sel = randperm(nr, howmanytraces);   % How many single traces to plot
        else 
            Sel = 1:length(howmanytraces);
        end
        hold on
        
        for bb=1:length(Sel)
            plot(0:6:492, temp(Sel(bb), :), '-', 'Color', [0.5 0.5 0.5, .4]);
        end
        %median
        plot(0:6:492, nanmedian(temp(:,:)), '-', 'Color', 'r', 'linewidth', 2);
        xline(120, '--', 'linewidth', 1.5)
        xline(240, '--', 'linewidth', 1.5)
        xline(360, '--', 'linewidth', 1.5)
        figtit = [];
        for cc = 1:4
            inputnum = unique( scmat(scmat(:, 167) == chambid(aa), 167+cc));
            if inputnum < 4
                figtit = strcat( figtit, " ", input_names( inputnum+1, :));
            elseif inputnum < 8 
                figtit = strcat( figtit, " ", input_names( inputnum-3, :));
            elseif inputnum < 12
                figtit = strcat( figtit, " ", input_names( inputnum-7, :));
            else
                figtit = strcat( figtit, " ", "FM");
            end
        end
                
        set(gca, 'Xlim', xlim, 'Ylim', ylim1);
        
    
        ax = gca;
        ax.XAxis.FontSize = 14;
        ax.YAxis.FontSize = 14;
               
        if ismember(aa, [26, 32, 38, 44] )
            ylabel("Nuc/Cyto p65", 'FontSize', 14)
        end

        if ismember(aa, 44:50 )
            xlabel("Time (min)", 'FontSize', 14)
        end
        
        title(figtit, 'FontSize', 14)
        

        hold off
    end
end
%% Fig S4: Low dose

scmat = scmatcomb_norm;
interval = 6;       % Time between each frame (in min)
xlim = [0, 480];
ylim1 = [-0.2 2.5];
howmanytraces = 20;

chambid = unique(scmat(:, 167));
input_names = ['T'; 'I'; 'L'; 'P'];


fig = figure(1); clf;
counter = 0;
for aa=51:74 %length(chambid)
    if ~ isempty(find(scmat(:, 167) == chambid(aa), 1))
        counter = counter+1;
        subplot(4,6,counter);
        
        temp = scmat(scmat(:, 167) == chambid(aa), 1:83);
        [nr, ~] = size(temp);
        if nr > howmanytraces
            Sel = randperm(nr, howmanytraces);   % How many single traces to plot
        else 
            Sel = 1:length(howmanytraces);
        end
        hold on
        
        for bb=1:length(Sel)
            plot(0:6:492, temp(Sel(bb), :), '-', 'Color', [0.5 0.5 0.5, .4]);
        end
        %median
        plot(0:6:492, nanmedian(temp(:,:)), '-', 'Color', 'r', 'linewidth', 2);
        xline(120, '--', 'linewidth', 1.5)
        xline(240, '--', 'linewidth', 1.5)
        xline(360, '--', 'linewidth', 1.5)
        figtit = [];
        for cc = 1:4
            inputnum = unique( scmat(scmat(:, 167) == chambid(aa), 167+cc));
            if inputnum < 4
                figtit = strcat( figtit, " ", input_names( inputnum+1, :));
            elseif inputnum < 8 
                figtit = strcat( figtit, " ", input_names( inputnum-3, :));
            elseif inputnum < 12
                figtit = strcat( figtit, " ", input_names( inputnum-7, :));
            else
                figtit = strcat( figtit, " ", "FM");
            end
        end
                
        set(gca, 'Xlim', xlim, 'Ylim', ylim1);
        
    
        ax = gca;
        ax.XAxis.FontSize = 14;
        ax.YAxis.FontSize = 14;
               
        if ismember(aa, [51, 57, 63, 69] )
            ylabel("Nuc/Cyto p65", 'FontSize', 14)
        end

        if ismember(aa, 69:74 )
            xlabel("Time (min)", 'FontSize', 14)
        end
        
        title(figtit, 'FontSize', 14)
        

        hold off
    end
end

%% Fig S5: peakheight heatmap
rng(2)

scmat = scmatcomb_norm(any(isnan(scmatcomb_norm),2)==0,:);
max_all = featmat(:, maxidx_all);
maxmat = [];
maxall_mat = [];
traceall_mat = [];

treatmat = [];
treatall_mat = [];
condlabs = {'T I L P', 'T I P L', 'T L I P', 'T L P I', 'T P I L', 'T P L I', 'I T L P', 'I T P L', 'I L T P', 'I L P T', ...
    'I P T L', 'I P L T', 'L T I P', 'L T P I', 'L I T P', 'L I P T', 'L P T I', 'L P I T', 'P T L I', 'P T I L', 'P I T L',...
    'P I L T', 'P L T I', 'P L I T'};
ylabs = {'TNF', 'IL-1', 'LPS', 'PAM'};
conditions = unique(scmat(:, 167));
for aa = [1:24, 26:49, 51:74]
    maxmat = [maxmat; [mean(max_all(scmat(:, 167) == conditions(aa), :)), conditions(aa)] ];
    if aa ~=63:68
        traceall_mat = [traceall_mat; sortrows(datasample( scmat(scmat(:, 167) == conditions(aa), 1:83), 50, 'Replace', false), 5, 'descend') ]; 
    else
        traceall_mat = [traceall_mat; sortrows(datasample( scmat(scmat(:, 167) == conditions(aa), 1:83), 50, 'Replace', false), 10, 'descend') ]; 

    end
    treatmat = [treatmat; scmat(find(scmat(:, 167) == conditions(aa), 1, 'first'), 168:172)];
    treatall_mat = [treatall_mat; scmat(find(scmat(:, 167) == conditions(aa), 50), 168:172)];
    
end
treatall_mat = [repmat(treatall_mat(:,1), [1, 20]), repmat(treatall_mat(:,2), [1, 20]), repmat(treatall_mat(:,3), [1, 20]),...
    repmat(treatall_mat(:,4), [1, 20]),  ones([length(treatall_mat(:,1)), 3])*99,  ];

normmat = [];
normall_mat = [];
traceall_mat = traceall_mat';
maxmat_trim = maxmat(:, 1:4)';
maxall_mat = maxall_mat';
foo = [];

redwhiteblue = [0 0 1];
for aa = 1:99
    temp = redwhiteblue(end, :);
    redwhiteblue = [redwhiteblue; [temp(1) + 0.01, temp(2) + 0.01, 1] ];
end

for aa = 1:98
    temp = redwhiteblue(end, :);
    redwhiteblue = [redwhiteblue; [1, temp(3) - 0.01, temp(3) - 0.01] ];
end

for aa = 1
    for bb = 0:3
        peak_l1 = mean( max_all(scmat(:, 172) == aa & scmat(:, 168) == bb + (aa-1)*4));
        normmat(:, bb+1)  = maxmat_trim( [treatmat(:,1:4) == bb + (aa-1)*4]')./peak_l1 ;
        normall_mat(:, (bb*20+1):(bb*20+20)) = reshape( traceall_mat( [treatall_mat(:,1:end) == bb + (aa-1)*4]')./peak_l1, [20, 1200])';
        
    end
    figure(1)
    heatmap(normmat + 0.01, 'CellLabelColor','none', 'Colormap',redwhiteblue, 'YDisplayLabels', condlabs, 'XDisplayLabels', ylabs,...
        'ColorScaling','log', 'ColorLimits', [-2 0.5])
    figure(4)
    normall_mat(normall_mat < 0) = 0;
    colormap(double( redwhiteblue) )
    imagesc(log10( normall_mat + 0.01), [-2, 0.5]) 
    set(gca,'XColor', 'none','YColor','none')
    colorbar

end

%% Fig S6: Features based on dose and ligand and time
maxidx_all = [1, 8, 15, 22];
overall_feat = featmat(:, maxidx_all);

ylab =["TNF", "IL-1", "LPS", "PAM"];
titles = { 'High', 'Moderate', 'Low'};

colors = hsv(5);

xnames = { '1', '2', '3', '4', 'NC'};

figure(1)
clf
plotcounter = 1;
zz=1;
xx=1;
position = [1, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12] ;
for aa = 1:3
    
    for bb = 1:4
        subplot(4, 3, position(plotcounter) )
        position(plotcounter)
        hold on
        counter = 1;
        for cc = 1:4
            if cc == 5
                tempmat = overall_feat(scmat(:, 168) == 12, 1);

            else
                tempmat = overall_feat(scmat(:, 172) == aa & ismember(scmat(:, 167+cc), [bb-1 bb+3 bb+7]), ...
                   cc);
                
                %for normalizing to first ligand
                peak_l1 = [ mean( overall_feat(scmat(:, 172) == aa & ismember(scmat(:, 168), [0 4 8])) ), ...
                   mean( overall_feat(scmat(:, 172) == aa & ismember(scmat(:, 168), [1 5 9])) ), ...
                   mean( overall_feat(scmat(:, 172) == aa & ismember(scmat(:, 168), [2 6 10])) ), ...
                   mean( overall_feat(scmat(:, 172) == aa & ismember(scmat(:, 168), [3 7 11])) ) ];
                tempmat = (overall_feat(scmat(:, 172) == aa & ismember(scmat(:, 167+cc), [bb-1 bb+3 bb+7]), ...
                   cc) - mean(overall_feat(scmat(:, 168) == 12, 1)) ) ./ (peak_l1(bb)- mean(overall_feat(scmat(:, 168) == 12, 1)) ) ;
            end
            
            Violin( (tempmat+0.1), cc, 'Bandwidth', 0.15, 'ViolinColor',[0 0 0], 'BoxColor',[0 0 0], 'ShowData', false, 'Width', .32, 'BoxWidth', 0.1, 'ShowMean', true); % This is actually showing median

            ax = gca;
            ax.XAxis.FontSize = 14;
            ax.YAxis.FontSize = 14;
            counter = counter+1;
            
        end
        set(gca, 'YLim', [-0.1, 3], 'XLim', [0.25, 4.5] )
        xticks([1 2 3 4 5])
        xticklabels(xnames)
        
        if ismember( position(plotcounter),  1:3 )
            title(titles(xx), 'FontSize', 16)
            xx = xx+1;
        end
        
        if ismember( position(plotcounter), [1, 4, 7, 10])
            ylabel(ylab(zz), 'FontSize', 16)
            zz=zz+1;
        end
        
        if position(plotcounter) > 9
            xlabel("Treatment position",'FontSize', 16)
        end
        plotcounter = plotcounter + 1;
        hold off
    end
    

end

hold off