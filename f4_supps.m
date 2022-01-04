%% Load data
clear
clc
load('scmat_MydKO.mat')
load('scmat_sequentialstim.mat')
load('featmat_sequentialstim.mat')
featmat_combnorm = featmat;
load("featmat_MydKO.mat")
load('scmat_LPSTS.mat')
load('scmatTS_norm.mat')
load('scmatLOqPCR_norm.mat')
load('featmat_JNK.mat')


maxidx_all = [1, 8, 15, 22];
aucidx_all = [7, 14, 21, 28];
maxtwoidx_all = [4, 11, 18, 15];
timeidx_all = [2, 9, 16, 23];
widthidx_all = [3, 10, 17, 24];


%% Fig 4A TNF permits subsequent MYD88 responses and MyD88 permits TNF

featmat = featmat_combnorm;
scmat = scmatcomb_norm;

overall_feat = featmat(:, maxidx_all);

titles =["D1 TNF"];

colors = parula(3);

xnames = { 'High', 'Mid', 'Low'};

figure(2)
clf
plotcounter = 1;
zz=1;
fun = {};
counter = 1;
allcells = {};


for aa = 1:3 %dose
    hold on
    
    for bb = 1:3 %ligand 1
        tempmat = [];
        peak_l1 = [ mean( overall_feat(scmat(:, 172) == aa & ismember(scmat(:, 168), [0 4 8])) ), ...
            mean( overall_feat(scmat(:, 172) == aa & ismember(scmat(:, 168), [1 5 9])) ), ...
            mean( overall_feat(scmat(:, 172) == aa & ismember(scmat(:, 168), [2 6 10])) ), ...
            mean( overall_feat(scmat(:, 172) == aa & ismember(scmat(:, 168), [3 7 11])) ) ];
        if bb == 1
            for cc = 2:4 %ligand 2
                tempmat2 = overall_feat(find(scmat(:, 172) == aa & ismember(scmat(:, 168), [bb-1 bb+3 bb+7]) &...
                    ismember(scmat(:, 169), [cc-1 cc+3 cc+7])), 2);
                
                tempmat2 = tempmat2 ./ peak_l1(cc) ;
                tempmat = [tempmat; tempmat2];
            end
            [nr, ] = size(tempmat)
            Violin(tempmat, counter, 'Bandwidth', 0.2, 'ViolinColor',[0 0 1 ], 'ViolinAlpha',0.3, 'BoxColor',[0 0 1], 'ShowData', false, 'Width', .35, 'BoxWidth', 0.1, 'ShowMean', true); % This is actually showing median
            fun{counter} = tempmat;
            
        elseif bb == 3
            for cc = 2:4 %ligand 2
                tempmat2 = overall_feat(find(scmat(:, 172) == aa & ismember(scmat(:, 168), [1:3 5:7 9:11]) &...
                    ismember(scmat(:, 169), [cc-1 cc+3 cc+7])), 2);
                tempmat2 = tempmat2 ./ peak_l1(cc) ;
                tempmat = [tempmat; tempmat2];
            end
            [nr, ] = size(tempmat)
            Violin(tempmat, counter, 'Bandwidth', 0.2, 'ViolinColor',[1 0 0 ],'ViolinAlpha',0.3, 'BoxColor',[1 0 0], 'ShowData', false, 'Width', .35, 'BoxWidth', 0.1, 'ShowMean', true); % This is actually showing median
            fun{counter} = tempmat;
        else
            for cc = 2:4
                
                tempmat2 = overall_feat(find(scmat(:, 172) == aa & ismember(scmat(:, 168), [cc-1 cc+3 cc+7]) &...
                    ismember(scmat(:, 169), [0, 4, 8])), 2);
                
                tempmat2 = tempmat2 ./ peak_l1(1) ;
                tempmat = [tempmat; tempmat2];
            end
            [nr, ] = size(tempmat)
            Violin(tempmat, counter, 'Bandwidth', 0.15, 'ViolinColor',[0 0.6 0 ], 'ViolinAlpha',0.3, 'BoxColor',[0 0 1], 'ShowData', false, 'Width', .35, 'BoxWidth', 0.1, 'ShowMean', true); % This is actually showing median
        end
        allcells{aa, bb} = tempmat;
        counter = counter  +1;
        
    end
    counter = counter  +1;
    
    rectangle('Position',[0.4, 1.97, 4.4, 1], 'FaceColor',[1 1 1])

    set(gca, 'YLim', [-0.1, 3], 'XLim', [0, 12.1])
    xticklabels(xnames)
    
    ax = gca;
    ax.XAxis.FontSize = 14;
    ax.YAxis.FontSize = 14;
    xticks([2 6 10 ])
    
    text(1, 2.8, 'S1 TNF \rightarrow S2 MyD88', 'FontSize', 12);
    text(0.5, 2.9, '___', 'Color', 'blue');

    text(1, 2.5, 'S1 MyD88 \rightarrow S2 TNF', 'FontSize', 12);
    text(0.5, 2.6, '___', 'Color', [0 0.6 0]);
    
    text(1, 2.2, 'S1 \rightarrow S2 MyD88', 'FontSize', 12);
    text(0.5, 2.3, '___', 'Color', 'red');
    
    ylabel("S2 max", 'FontSize', 14)
    xlabel("Dose", 'FontSize', 14)

    
    hold off

end
hold off


pval_f1e = [];
for aa = 1:3
    for bb = 1:2
        if bb == 1
            p = ranksum(allcells{aa, bb}, allcells{aa, bb+2});
        else
            p = ranksum(allcells{aa, bb}, allcells{aa, bb+1});
        end
        pval_f1e(aa, bb) = p*6;
        
    end
    
end

esize_f1e = [];
for aa = 1:3
    for bb = 1:2
        if bb == 1
            e =  mean(allcells{aa, bb})./mean(allcells{aa, bb+2});
        else
            e =  mean(allcells{aa, bb})./mean(allcells{aa, bb+1});
        end
        
        esize_f1e(aa, bb) = e; 
    end
    
end

%% Fig 4B IL-1 and PAM exhibit dose dependent cross-inhibition of other Myd88-dep ligands, LPS does not

overall_feat = featmat(:, maxidx_all);


xlabels =["Dose"];

colors = parula(3);

xnames = { 'High', 'Mid', 'Low'};

figure(2)
clf
plotcounter = 1;
counter = 1;

zz=1;
for aa = 1:3 

    peak_l1 = [ mean( overall_feat(scmat(:, 172) == aa & ismember(scmat(:, 168), [0 4 8])) ), ...
        mean( overall_feat(scmat(:, 172) == aa & ismember(scmat(:, 168), [1 5 9])) ), ...
        mean( overall_feat(scmat(:, 172) == aa & ismember(scmat(:, 168), [2 6 10])) ), ...
        mean( overall_feat(scmat(:, 172) == aa & ismember(scmat(:, 168), [3 7 11])) ) ];

    counter2 = 1;
    for bb = [2, 4, 3] 
        tempmat = [];
        hold on

        for cc = 1:3 

            tempmat2 = overall_feat(find(scmat(:, 172) == aa & ismember(scmat(:, 168), [bb-1 bb+3 bb+7]) &...
                ismember(scmat(:, 169), [cc cc+4 cc+8])), 2);

            tempmat2 = tempmat2 ./ peak_l1(cc+1) ;
            tempmat = [tempmat; tempmat2];
        end
        [nre, ~] = size(tempmat)
        if bb == 2
            color = [0 0 1];
        elseif bb == 4
            color = [0 0.6 0];
        else
            color = [1 0 0];
        end

        Violin(tempmat, counter, 'Bandwidth', 0.15, 'ViolinColor',color, 'BoxColor',[0 0 0], 'ShowData', false, 'Width', .4, 'BoxWidth', 0.05, 'ShowMean', true); % This is actually showing median

        counter = counter+1;

        set(gca, 'YLim', [-0.1, 2.5], 'XLim', [0, 12])
        xticks([2 6 10 ])
        xticklabels(xnames)

        allcells{aa, counter2} = tempmat;
        counter2 = counter2+1;


        if plotcounter < 4
            xlabel(xlabels(plotcounter), 'FontSize', 14)
        end


        ylabel("S2 max", 'FontSize', 14)
        %xlabel("Dose",'FontSize', 14)
        ax = gca;
        ax.XAxis.FontSize = 14;
        ax.YAxis.FontSize = 14;

        hold off


    end
    counter = counter + 1;

    rectangle('Position',[0.4, 1.47, 2.3, 1], 'FaceColor',[1 1 1])


    text(1, 2.3, 'S1 IL-1', 'FontSize', 12);
    text(0.5, 2.36, '___', 'Color', 'blue');

    text(1, 2.0, 'S1 PAM', 'FontSize', 12);
    text(0.5, 2.06, '___', 'Color', [0 0.6 0]);

    text(1, 1.7, 'S1 LPS', 'FontSize', 12);
    text(0.5, 1.76, '___', 'Color', 'red');


end

pval_f1e = [];
for aa = 1:3
    for bb = 1:2
        if bb == 1
            p = ranksum(allcells{aa, bb}, allcells{aa, bb+2});
        else
            p = ranksum(allcells{aa, bb}, allcells{aa, bb+1});
        end
        pval_f1e(aa, bb) = p*6;
        
    end
    
end

esize_f1e = [];
for aa = 1:3
    for bb = 1:2
        if bb == 1
            e =  mean(allcells{aa, bb})./mean(allcells{aa, bb+2});
        else
            e =  mean(allcells{aa, bb})./mean(allcells{aa, bb+1});
        end
        
        esize_f1e(aa, bb) = e; 
    end
    
end

%mean(allcells{1, 3})./mean(allcells{1, 1})
hold off

%% Fig 4C means for LPS
condmat = [14, 1, 2, 3, 4, 5, 6;
    15, 7, 8, 9, 10, 11, 13;
    29, 17, 18, 19, 20, 21, 22;
    30, 23, 24, 25, 26, 27, 28];
d2_times = [0, 10, 20, 40, 60, 90, 120];


xlim = [0, 180];
ylim1 = [-0.1 6];
scmat = scmat_LPSTS;


titles = {'12.5'; '25'; '100'; '400' };

counter = 0;
colors = parula(7);

figure(1)
clf

for aa=1:4
    c2 = 0;
    for bb = [2, 4, 5, 6, 7, 1]
        counter = counter+1;
        subplot(4,6,counter);
        temp = scmat(scmat(:, 82) == condmat(aa, bb), 1:81);
        c2  =c2 + 1;
        [nr, ~] = size(temp);
        hold on
        if bb == 1
            plot(0:3:240, nanmean(temp(:,:)), '-', 'Color', [0.6 0.6 0.6], 'linewidth', 2);
        else
            plot(0:3:240, nanmean(temp(:,:)), '-', 'Color', 'r', 'linewidth', 2);
        end
        if ismember( counter, [1:5. 7:11, 13:17, 19:23]) 
            xline(d2_times(bb), '--', 'Color', 'r', 'linewidth', 2 )
        end
        ax = gca;
        ax.XAxis.FontSize = 16;
        ax.YAxis.FontSize = 16;
        set(gca, 'Xlim', xlim, 'Ylim', ylim1);


        title(titles{aa}, 'FontSize', 16)
        if ismember( counter, [1, 7, 13, 19])
            ylabel("Nuc/Cyto RelA", 'FontSize', 16)
        end
    end

    c2 =0;
end

%% Fig 4D superposition of medians for IL-1

scmatIL_1 = scmatTS_norm;
scmatIL_1_2 = scmatLOqPCR_norm;
scmat_og = scmatcomb_norm;
colors = parula(7);


chambersIL_1 = [16, 1, 17, 2, 18, 3, 4, 20];
chambersIL_1_2 = [1, 3, 5];
titles2 = {'3'; '0.2'};
xlim = [0, 180];
ylim1 = [-0.1 6];

figure(1)
clf
d2_times2 = [30, 60, 120, 0];
for aa = [2, 1]
    subplot(1, 2, aa);
    c2 = 1;
    if aa == 1
        for bb = [3, 4, 6, 7, 2]
            temp = scmatIL_1(scmatIL_1(:, 62) == chambersIL_1(bb), 1:61) ;
            [nr, ~] = size(temp)
            hold on
            if bb == 2
                plot(0:3:180, nanmean(temp(:,:)), '-', 'Color', 'r', 'linewidth', 2);
            else
                plot(0:3:180, nanmean(temp(:,:)), '-', 'Color', colors(c2, :), 'linewidth', 2);
            end
            c2 = c2 + 1;
            %xline(d2_times2(cc), '--', 'Color', colors(c2, :) )
        end
        for cc = 1:4
            xline(d2_times2(cc), '--', 'Color', colors(cc, :) )
        end
 
    else
        for bb = [1, 3, 7, 5]
            temp = scmatIL_1_2(scmatIL_1_2(:, 27) == bb, 1:26) ;
            c2  = c2 + 1;
            [nr, ~] = size(temp)
            hold on
            if bb == 5
                plot(0:6:144, nanmean(temp(:,1:25)), '-', 'Color', 'r', 'linewidth', 2);
            elseif bb == 7
                temp = scmat_og(scmat_og(:, 168) == 5 & scmat_og(:, 169) == 6, 1:31);
                plot(-3:6:177, nanmean(temp(:,:)), '-', 'Color', colors(c2, :), 'linewidth', 2);
            else
                plot(0:6:150, nanmean(temp(:,:)), '-', 'Color', colors(c2, :), 'linewidth', 2);
            end

        end
        legend( {'30 min', '60 min', '120 min', 'Control' }, 'AutoUpdate',  'off', 'FontSize', 12)
        for cc = 1:4
            xline(d2_times2(cc), '--', 'Color', colors(cc, :) )
        end
    end
    if aa == 1
        ylabel("Nuc/Cyto RelA", 'FontSize', 16)
    end

    %xlabel("Time (min)",'FontSize', 16)
    ax = gca;
    ax.XAxis.FontSize = 16;
    ax.YAxis.FontSize = 16;
    title(titles2{aa}, 'FontSize', 16)
    
    set(gca, 'Xlim', xlim, 'Ylim', ylim1);
    
    hold off
end
%% Fig 4E 2nd peak height
scmat = scmat_LPSTS;


scndpeakht = [scmat(:, 82), zeros(length(scmat(:, 82)), 1)];


d2_times = [10, 20, 40, 60, 90, 120, 10, 20, 40, 60, 90, 120, 120, 120, 120, 120,...
    10, 20, 40, 60, 90, 120, 10, 20, 40, 60, 90, 120, 120, 120, 120, 120];

d2_timesIL_1 = [30, 60, 120, 120];

for aa = 1:length(scmat)
        tnew = ceil( d2_times(scmat(aa, 82) )/3) + 5;
        tend = tnew+20;
        scndpeakht(aa, 2) = max(scmat(aa, tnew:tend));
        
end




xnames = { '10', '20', '40', '60', '90', '120'};

counter = 1;
figure(6)
clf
hold on
chambers = [1, 2, 3, 4, 5, 6;
            7, 8, 9, 10, 11, 13;
            17, 18, 19, 20, 21, 22;
            23, 24, 25, 26, 27, 28];
for aa = 1:6
    row = 3;
    normmax = mean(scndpeakht(ismember( scndpeakht(:, 1), [16, 18]) , 2));
    tempmat = scndpeakht(scndpeakht(:, 1) == chambers(row, aa), 2)./normmax;
    [ncells, ~] = size(tempmat);
    
    if ncells > 100
        temp_sample = datasample(tempmat, 100, 'Replace', false);
    else
        temp_sample = tempmat;
    end
    Violin(tempmat, aa, 'Bandwidth', 0.2, ...
    'ViolinColor',[0 0 0 ], 'BoxColor',[0 0 0], 'ShowData', false, 'Width', .3, 'BoxWidth', 0.1, 'ShowMean', true); % This is actually showing median
    
    counter = counter+1;

end

set(gca, 'YLim', [-0.01, 2], 'XLim', [0.25, 6.5] )
xticks([1 2 3 4 5 6 ])
xticklabels(xnames)
ylabel("IL-1 max", 'FontSize', 14)
xlabel("Time after LPS (min)", 'FontSize', 14)


%xlabel("Time (min)",'FontSize', 16)
ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;

if row == 1
    title("12.5", 'FontSize', 14)
elseif row == 2
    title("25", 'FontSize', 14)
elseif row == 3
    title("100 ng/mL LPS", 'FontSize', 14)
elseif row == 4
    title("400", 'FontSize', 14)
end
hold off

%% Fig 4E 2nd peak height
scmatIL_1_2 = scmatLOqPCR_norm;
scmatIL_1_2 = scmatIL_1_2(ismember(scmatIL_1_2(:, 27), IL_indices) , :);
scmat_og = scmatcomb_norm;

IL_indices = [1, 3, 7, 5];
scndpeakhtIL_1 = [scmatIL_1_2(:, 27), zeros(length(scmatIL_1_2(:, 27)), 1)];


overall_feat = featmat(:, maxidx_all);
for aa = 1:length(scmatIL_1_2)
        tnew = ceil( d2_timesIL_1(IL_indices == scmatIL_1_2(aa, 27) )/6);
        tend = tnew+20;
        scndpeakhtIL_1(aa, 2) = max(scmatIL_1_2(aa, tnew:26));
        
end

xnames = { '30', '60', '120'};

counter = 1;
figure(6)
clf
hold on
chambers = [1, 3, 5];
for aa = 1:3
    normmax = mean( overall_feat(ismember(scmat_og(:, 168), 6) ) );
    if aa == 3
        tempmat = overall_feat(scmat_og(:, 168) == 5 & scmat_og(:, 169) == 6, 2)./normmax;
        
    else
        tempmat = scndpeakhtIL_1(scndpeakhtIL_1(:, 1) == chambers(aa), 2)./normmax;
    end

    Violin(tempmat, aa, 'Bandwidth', 0.2, ...
    'ViolinColor',[0 0 0 ], 'BoxColor',[0 0 0], 'ShowData', false, 'Width', .3, 'BoxWidth', 0.1, 'ShowMean', true); % This is actually showing median
    
    counter = counter+1;

end

set(gca, 'YLim', [-0.01, 2], 'XLim', [0.25, 3.75] )
xticks([1 2 3])
xticklabels(xnames)
ylabel("LPS max", 'FontSize', 14)
xlabel("Time after IL-1 (min)", 'FontSize', 14)
title("0.2 ng/mL IL-1", 'FontSize', 14)



%xlabel("Time (min)",'FontSize', 16)
ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;


hold off

%% fS7A MYD peak height relative to WT peak height

WTht = [ mean( featmat_combnorm(ismember(scmatcomb_norm(:, 168), 0), 1) ), ...
    mean( featmat_combnorm(ismember(scmatcomb_norm(:, 168), 1), 1) ), ...
    mean( featmat_combnorm(ismember(scmatcomb_norm(:, 168), 2), 1) ), ...
    mean( featmat_combnorm(ismember(scmatcomb_norm(:, 168), 3), 1) )];

KOht = featmat_KO(:, 1);

scmat = scmatMydTrif_norm(any(isnan(scmatMydTrif_norm),2)==0,:);

xnames = { 'TNF', 'IL-1', 'LPS', 'PAM'};

counter = 1;
figure(2)
clf
hold on
for aa = 1:4
    tempmat = KOht(scmat(:, 87) == 32+aa)./WTht(aa);    
    [~, ncells] = size(tempmat);
    
    Violin( (tempmat), aa, 'Bandwidth', 0.2, 'ViolinColor',[0 0 0], 'BoxColor',[0 0 0], 'ShowData', false, 'Width', .32, 'BoxWidth', 0.1, 'ShowMean', true); 
    
    counter = counter+1;

end

set(gca, 'YLim', [-0.01, 3], 'XLim', [0.25, 4.75] )
xticks([1 2 3 4])
xticklabels(xnames)
xlabel("Ligand",'FontSize', 14)
ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
ylabel("MyD88 KO max response", 'FontSize', 14)
%title("Myd88 KO relative to WT", 'FontSize', 16)

hold off

%% fS7A MyD88 KO traces

interval = 6;       % Time between each frame (in min)
xlim = [0, 240];
ylim1 = [-0.2 4];
howmanytraces = 20;

scmat = scmatMydTrif_norm;

chambid = unique(scmat(:, 87));
input_names = ['TNF 30', 'IL-1 0.2', 'LPS 100', 'PAM 0.1', 'TNF 90', 'IL-1 3', 'LPS 400', 'PAM 1'];

ytit = ['TNF', "IL-1", "LPS", 'PAM'];
fig = figure(1); clf;
counter = 0;
for aa=5:8
    if ~ isempty(find(scmat(:, 87) == chambid(aa), 1))
        counter = counter+1;
        subplot(1, 4,counter);
        
        temp = scmat(scmat(:, 87) == chambid(aa), 1:41);
        [nr, ~] = size(temp);
        if nr > howmanytraces
            Sel = randperm(nr, howmanytraces);   % How many single traces to plot
        else 
            Sel = 1:length(howmanytraces);
        end
        hold on
        
        for bb=1:length(Sel)
            plot(0:6:240, temp(Sel(bb), :), '-', 'Color', [0.5 0.5 0.5, .4]);
        end
        %median
        plot(0:6:240, nanmean(temp(:,:)), '-', 'Color', 'r', 'linewidth', 2);
        %title(names(counter, 1))
        figtit = [];
        for aa = 1:2
           inputnum = unique( scmat(scmat(:, 85) == chambid(aa), 83+aa));
            if inputnum < 8
                figtit = strcat( figtit, " ", input_names(inputnum+1));
            elseif inputnum > 8
                figtit = strcat(figtit, " ", "FM");
            end
        end
        
        ylabel(ytit(counter))

        set(gca, 'Xlim', xlim, 'Ylim', ylim1);
               
        ax = gca;
        xlabel("Time (min)", 'FontSize', 12)
        ax.XAxis.FontSize = 14;
        ax.YAxis.FontSize = 14;

        hold off
    end
end

%% fS8 JNK

maxidx_all = [1, 8, 15, 22];

overall_feat = featmat(:, maxidx_all);
overall_feat(overall_feat < 0.5) = 0;
overall_feat = [overall_feat(overall_feat(:, 1) ~= 0, :), featmat(overall_feat(:, 1) ~= 0, 29:34)] ;

titles =["JNK response"];

xnames = { 'High', 'Mid'};

figure(3)
clf
plotcounter = 1;
counter = 1;
zz=1;
for aa = 1:2  %dose
    
    for bb = 2:3 %ligand 1
        tempmat = [];
        plotcounter = [1, 3, 2];
        hold on
        peak_l1 = [ mean( overall_feat(overall_feat(:, 10) == aa & ismember(overall_feat(:, 6), [0 4 8])) ), ...
            mean( overall_feat(overall_feat(:, 10) == aa & ismember(overall_feat(:, 6), [1 5 9])) ), ...
            mean( overall_feat(overall_feat(:, 10) == aa & ismember(overall_feat(:, 6), [2 6 10])) ), ...
            mean( overall_feat(overall_feat(:, 10) == aa & ismember(overall_feat(:, 6), [3 7 11])) ) ];
        for cc = 1:3 %ligand 2
            tempmat2 = overall_feat(find(overall_feat(:, 10) == aa & ismember(overall_feat(:, 6), [bb-1 bb+3 bb+7]) &...
                ismember(overall_feat(:, 7), [cc cc+4 cc+8])), 2);
            tempmat2 = tempmat2( tempmat2 < 7) ;

            tempmat2 = tempmat2 ./ peak_l1(cc+1) ;
            tempmat = [tempmat; tempmat2];
        end
        [nr, ~] = size(tempmat) ;
        if ismember(counter, [1, 3])
            Violin(tempmat, counter, 'ViolinColor',[0 0 1 ], 'BoxColor',[0 0 1], 'ViolinAlpha',0.3, 'Bandwidth', 0.15, 'ShowData', false, 'Width', .32, 'BoxWidth', 0.075, 'ShowMean', true);
        else
            Violin(tempmat, counter, 'ViolinColor',[1 0 0 ], 'BoxColor',[1 0 0], 'ViolinAlpha',0.3, 'Bandwidth', 0.15, 'ShowData', false, 'Width', .32, 'BoxWidth', 0.075, 'ShowMean', true);
        end
        allcells{aa, bb-1} = tempmat;

       

        counter = counter+1;
            
        set(gca, 'YLim', [-0.1, 2.5], 'XLim', [0.5, 4.5])
        xticks([1.5 3.5 ])
        xticklabels(xnames)
        ax = gca;
        ax.XAxis.FontSize = 14;
        ax.YAxis.FontSize = 14;      
        
        ylabel("JNK S2 max", 'FontSize', 14)
        xlabel("Dose", 'FontSize', 14)
        rectangle('Position',[.95, 1.95, 1.4, .5], 'FaceColor',[1 1 1])

        
        text(1.1, 2.38, '___', 'Color', 'blue','FontSize', 12, 'FontWeight','bold');
        text(1.6, 2.3, 'S1 IL-1','FontSize', 12, 'FontWeight','normal');

        
        text(1.1, 2.18, '___', 'Color', 'Red', 'FontSize', 12, 'FontWeight','bold');
        text(1.6, 2.1, 'S1 LPS', 'FontSize', 12, 'FontWeight','normal');

        
        hold off

    end
    plotcounter = plotcounter + 1;


end
hold off

pval_f1e = [];
for aa = 1:2
    for bb = 1
        if bb == 1
            p = ranksum(allcells{aa, bb}, allcells{aa, bb+1});
        end
        pval_f1e(aa, bb) = p*2;
        
    end
    
end

esize_f1e = [];
for aa = 1:2
    for bb = 1
        if bb == 1
            e =  mean(allcells{aa, bb})./mean(allcells{aa, bb+1});
        end
        
        esize_f1e(aa, bb) = e; 
    end
    
end

%% fS9 LPS titration raw peaks 

interval = 3;       % Time between each frame (in min)
xlim = [0, 240];
ylim1 = [-1 7];
howmanytraces = 20;

scmat = scmat_LPSTS; %scmatnorm;

condmat = [14, 1, 2, 3, 4, 5, 6, 15, 7, 8, 9, 10, 11, 13, 29, 17, 18, 19, 20, 21, 22, 30, 23, 24, 25, 26, 27, 28];

chambid = unique(scmat(:, 82));
d2_times = [-10, 10, 20, 40, 60, 90, 120, -10, 10, 20, 40, 60, 90, 120, ...
   -10, 10, 20, 40, 60, 90, 120, -10, 10, 20, 40, 60, 90, 120];

titles = {"No IL-1", "10 min", "20 min", "40 min", "60 min",  "90 min", "120 min"};
ylabs = {"12.5 ng", "25 ng", "100 ng", "400 ng"};

fig = figure(1); clf;
%set(fig,'defaultAxesColorOrder',[[1, 0, 0]; [0, 0, 1]]);
counter = 0;
zz = 1;
for aa=1:length(condmat)
%for blah=1:4
%    aa = chambers(blah);
    if ~ isempty(find(scmat(:, 82) == condmat(aa), 1))
        counter = counter+1;
        subplot(4, 7,counter);
        
        temp = scmat(scmat(:, 82) == condmat(aa), 1:81);
        [nr, ~] = size(temp);
        if nr > howmanytraces
            Sel = randperm(nr, howmanytraces);   % How many single traces to plot
        else 
            Sel = 1:length(howmanytraces);
        end
        hold on
        
        for bb=1:length(Sel)
            plot(0:3:240, temp(Sel(bb), :), '-', 'Color', [0.5 0.5 0.5, .4]);
        end
        %median
        plot(0:3:240, nanmean(temp(:,:)), '-', 'Color', 'r', 'linewidth', 2);
        xline(d2_times(aa), '--', 'linewidth', 2);
        
        if ismember( counter, 1:7)
            title(titles{counter}, 'FontSize', 14)
        end
        
        if ismember( counter, [1, 8, 15, 22 ])
            ylabel(ylabs{zz}, 'FontSize', 14)
            zz = zz+1;
        end
        
        if ismember( counter, 22:28)
            xlabel("Time (min)", 'FontSize', 14)
        end

        
        ax = gca;
        ax.XAxis.FontSize = 10;
        ax.YAxis.FontSize = 10;

        set(gca, 'Xlim', xlim, 'Ylim', ylim1);
        
        hold off
    end
end

%% f9b IL-1 titration raw peaks 

interval = 3;       % Time between each frame (in min)
howmanytraces = 20;

scmatIL_1 = scmatTS_norm;
scmatIL_1_2 = scmatLOqPCR_norm;
scmat_og = scmatcomb_norm;

chambersIL_1 = [16, 1, 17, 2, 18, 3, 4, 20];
chambersIL_1_2 = [1, 3, 5];

titles = {"No LPS", "30 min", "60 min", "90 min"};

xlim = [0, 180];
ylim1 = [-0.1 7];

figure(3)
clf
d2_times = [0, 30, 60, 120];

chambid = unique(scmatIL_1(:, 62));

counter = 0;
zz = 1;
for aa=[2, 4, 6, 7]
    if ~ isempty(find(scmatIL_1(:, 62) == chambersIL_1(aa), 1))
        counter = counter+1;
        subplot(2, 4,counter);
        
        temp = scmatIL_1(scmatIL_1(:, 62) == chambersIL_1(aa), 1:61);
        [nr, ~] = size(temp);
        if nr > howmanytraces
            Sel = randperm(nr, howmanytraces);   % How many single traces to plot
        else 
            Sel = 1:length(howmanytraces);
        end
        hold on
        
        for bb=1:length(Sel)
            plot(0:3:180, temp(Sel(bb), :), '-', 'Color', [0.5 0.5 0.5, .4]);
        end
        %median
        plot(0:3:180, nanmean(temp(:,:)), '-', 'Color', 'r', 'linewidth', 2);
        if counter > 1
            xline(d2_times(counter), '--', 'linewidth', 2);
        end
        set(gca, 'Xlim', xlim, 'Ylim', ylim1);
        ax = gca;


        title(titles{counter}, 'FontSize', 14)

        
        if counter == 1
            ylabel("3 ng/mL", 'FontSize', 14)
        end
        
        xlabel("Time", 'FontSize', 14)
        ax.XAxis.FontSize = 14;
        ax.YAxis.FontSize = 14;
        
        hold off
    end
end

figure(2)
clf
d2_times = [0, 30, 60, 120];

chambid = unique(scmatIL_1_2(:, 26));

counter = 0;
zz = 1;
for aa=[5, 1, 3, 7]
    counter = counter+1
    subplot(2, 4,counter);
    
    hold on
    c2  = c2 + 1;
    if aa == 5
        temp = scmatIL_1_2(scmatIL_1_2(:, 27) == aa, 1:26) ;
        [nr, ~] = size(temp);
        if nr > howmanytraces
            Sel = randperm(nr, howmanytraces);   % How many single traces to plot
        else
            Sel = 1:length(howmanytraces);
        end
        plot(0:6:144, nanmean(temp(:,1:25)), '-', 'Color', 'r', 'linewidth', 2);
        for bb=1:length(Sel)
            plot(0:6:144, temp(Sel(bb), 1:25), '-', 'Color', [0.5 0.5 0.5, .4]);
        end
    elseif aa == 7
        temp = scmat_og(scmat_og(:, 168) == 5 & scmat_og(:, 169) == 6, 1:31);
        [nr, ~] = size(temp);
        if nr > howmanytraces
            Sel = randperm(nr, howmanytraces);   % How many single traces to plot
        else
            Sel = 1:length(howmanytraces);
        end
        plot(-3:6:177, nanmean(temp(:,:)), '-', 'Color',  'r', 'linewidth', 2);
        for bb=1:length(Sel)
            plot(-3:6:177, temp(Sel(bb), :), '-', 'Color', [0.5 0.5 0.5, .4]);
        end
    else
        temp = scmatIL_1_2(scmatIL_1_2(:, 27) == aa, 1:26) ;
        [nr, ~] = size(temp);
        if nr > howmanytraces
            Sel = randperm(nr, howmanytraces);   % How many single traces to plot
        else
            Sel = 1:length(howmanytraces);
        end
        plot(0:6:150, nanmean(temp(:,:)), '-', 'Color',  'r', 'linewidth', 2);
        for bb=1:length(Sel)
            plot(0:6:150, temp(Sel(bb), :), '-', 'Color', [0.5 0.5 0.5, .4]);
        end
        
    end
    if counter > 1
        xline(d2_times(counter), '--', 'linewidth', 2);
    end
    
    set(gca, 'Xlim', xlim, 'Ylim', ylim1);
    ax = gca;
    
    
    title(titles{counter}, 'FontSize', 14)
    
    
    if counter == 1
        ylabel("0.2 ng/mL", 'FontSize', 14)
    end
    
    xlabel("Time", 'FontSize', 14)
    ax.XAxis.FontSize = 14;
    ax.YAxis.FontSize = 14;
    
    hold off
end