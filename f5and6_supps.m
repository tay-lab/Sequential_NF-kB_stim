%%
clear
clc
load('scmat_PS1145v2.mat')
load('scmat_LOqPCR_2.mat')
load('scmat_sequentialstim.mat')
load('featmat_sequentialstim.mat')
load('featmat_qPCR.mat')


maxidx_all = [1, 8, 15, 22];
aucidx_all = [7, 14, 21, 28];
maxtwoidx_all = [4, 11, 18, 15];
timeidx_all = [2, 9, 16, 23];
widthidx_all = [3, 10, 17, 24];

%% Fig 5a - ikki
scmat = scmat_PS1145v2;
scndpeakht = [scmat(:, 163), zeros(length(scmat(:, 163)), 1), zeros(length(scmat(:, 163)), 1)];
d2_times = [180, 300, 180, 300, 180, 300, 180, 300, 180, 300, 180, 300, 180, 300, 180, 300];
for aa = 1:length(scmat)
        tnew = ceil( d2_times(scmat(aa, 163) )/6);
        tend = tnew+20;
        scndpeakht(aa, 2) = max(scmat(aa, 1:tnew-1));
        scndpeakht(aa, 3) = max(scmat(aa, tnew+2:tend+2));
        
end

colors = hsv(8);

xnames = { '25', '100', '400'};

counter = 1;
figure(1)
clf
hold on
chambers = [8, 14, 10, 16, 12, 6;
            7, 13, 9, 15, 11, 5;
            ];
for aa = 1:6
    row = 1;
    normmax = mean(scndpeakht(ismember( scndpeakht(:, 1), [1, 3]) , 3));
    
    
    %IL-1 response
    tempmat = scndpeakht(scndpeakht(:, 1) == chambers(row, aa), 3)./normmax;
    [nr, ~] = size(tempmat)
    if ismember(aa, [1, 3, 5])
        Violin(tempmat, aa, 'ViolinColor',[0 0 1 ], 'ViolinAlpha',0.3, 'Bandwidth', 0.15, 'BoxColor',[0 0 1], 'ShowData', false, 'Width', .32, 'BoxWidth', 0.075, 'ShowMean', true);
    else
        Violin(tempmat, aa, 'ViolinColor',[1 0 0 ], 'ViolinAlpha',0.3, 'Bandwidth', 0.15, 'BoxColor',[1 0 0], 'ShowData', false, 'Width', .32, 'BoxWidth', 0.075, 'ShowMean', true);
    end
    allcells{aa} = tempmat;
    counter = counter+1;

end

set(gca, 'YLim', [-0.01, 1.5], 'XLim', [0.25, 6.5] )

xticklabels(xnames)
    
xticks([1.5 3.5 5.5 ])
rectangle('Position',[4.1, 1.12, 2, .3], 'FaceColor',[1 1 1])


text(4.8, 1.35, 'No Drug','FontSize', 12);
text(4.3, 1.4, '___', 'Color', 'blue','FontSize', 12, 'FontWeight','Bold');

text(4.8, 1.2, '+ PS1145', 'FontSize', 12);
text(4.3, 1.25, '___', 'Color', 'Red', 'FontSize', 12, 'FontWeight','Bold');


xticks([1.5 3.5 5.5 ])
xticklabels(xnames)
ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
xlabel("LPS 1st dose (ng/mL)",'FontSize', 14)
ylabel("IL-1 max", 'FontSize', 14)

hold off

xnames = { '25ng ', '100ng', '400ng'};

pval_f1e = [];
for aa = 1:3
    p = ranksum(allcells{aa*2-1}, allcells{aa*2});
    pval_f1e(aa) = p*3;
    
end

esize_f1e = [];
for aa = 1:3
    e =  mean(allcells{aa*2})./mean(allcells{aa*2-1});        
    esize_f1e(aa) = e; 
    
end


%% S12A - peak widths
scmat = scmatcomb_norm(any(isnan(scmatcomb_norm),2)==0,:);
overall_feat = featmat(:, widthidx_all(1) );
titles = {"High", "Mid", "Low"};

colors = hsv(7);

xnames = { 'TNF', 'IL-1', 'LPS', 'PAM', 'NC'};

counter = 1;
figure(2)
clf
hold on
counter = 1;
for bb = 1:3
    subplot(1, 3, bb)
    hold on
    for aa = 1:4
        if bb == 2
            cc = aa+4;
        elseif bb == 3
            cc = aa+8;
        else
            cc = aa;
        end

        tempmat = overall_feat(ismember(scmat(:, 168), [cc-1 ]) & ...
            overall_feat~= 0 )*6;
        
        Violin( (tempmat), aa, 'Bandwidth', 0.8, 'ViolinColor',[0 0 0], 'BoxColor',[0 0 0], 'ShowData', false, 'Width', .32, 'BoxWidth', 0.1, 'ShowMean', true);

        hold off
    end
    set(gca, 'YLim', [0, 50], 'XLim', [0.25, 4.5] )
    xticks([1:5])
    xticklabels(xnames)
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    xlabel("Ligand ID",'FontSize', 14)
    ylabel("Peak width(min)", 'FontSize', 14)
    title(titles(bb), 'FontSize', 14)
    counter = counter+1;

end


hold off

%% S12B - peak AUCs

LPSfeat = featmat_combnorm(ismember( scmatcomb_norm(:, 168), 6), [1, 7]);
IL_1mfeat = featmat_qPCR(ismember( scmat_LOqPCR_20201116_2(:, 42), 2), 2:3);
IL_1mhfeat = featmat_qPCR(ismember( scmat_LOqPCR_20201116_2(:, 42), 1), 2:3);

xnames = { '0.2 ng IL-1', '1 ng IL-1', 'LPS'};

colors = hsv(5);
figure(1) 
clf
hold on
counter = 1;
for aa = 1:3
    if aa == 3
        tempmat = LPSfeat(:, 1);  
    elseif aa == 1
        tempmat = IL_1mfeat(:, 1); 
    else
        tempmat = IL_1mhfeat(:, 1);
    end
    
    Violin( (tempmat), aa, 'Bandwidth', 0.3,  'ViolinColor', [0 0 0], 'ViolinAlpha', .3, 'BoxColor',[.25 .25 .25],...
        'ShowData', false, 'Width', .32, 'BoxWidth', .1, 'ShowMean', true);
    
    xticklabels(xnames)
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;  
    xlabel("Stimulus",'FontSize', 14)
    ylabel("AUC", 'FontSize', 14)
    xticks([1 2 3 ])
    set(gca, 'XLim', [0.25, 3.75], 'YLim', [0, 12])
    
    counter = counter+1;

end

hold off

