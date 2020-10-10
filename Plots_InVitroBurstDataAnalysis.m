function Plots_InVitroBurstDataAnalysis()
% 31st December 2019
load('C:\Users\Deepcutlab\Desktop\Nidhi\Code_Nidhi\In-Vitro Burst Analysis\Results\inVitroBurstDataTable.mat');
lowerOrderIndx = strcmp(inVitroBurstDataTable.cellOrder,'lowerOrder');
higherOrderIndx = strcmp(inVitroBurstDataTable.cellOrder,'higherOrder');
rowsNumSpksMoreThanOne = inVitroBurstDataTable.NumSpikes > 1;
edges = [-120;-110;-100;-90;-80;-70;-60]; % only take traces with voltages between -120mV and -60mV

%% Number of Spikes
%% Plotting NumSpikes versus hyperpolarization voltage
lowerOrderIndx = strcmp(inVitroBurstDataTable.cellOrder,'lowerOrder');
lowerOrderTable = inVitroBurstDataTable(lowerOrderIndx,[8,14,20]);
hyperVolL = lowerOrderTable.hyperpolarVoltage;
NumSpikesL = lowerOrderTable.NumSpikes;
figure(1); plot(hyperVolL,NumSpikesL,'o');
yticks(0:1:10); ylim([-1,10]);
xlabel('hyperpolarization voltage (in mV)');
ylabel('number of spikes in a burst');
title('Lower Order');

higherOrderIndx = strcmp(inVitroBurstDataTable.cellOrder,'higherOrder');
higherOrderTable = inVitroBurstDataTable(higherOrderIndx,[8,14,20]);
hyperVolH = higherOrderTable.hyperpolarVoltage;
NumSpikesH = higherOrderTable.NumSpikes;
[tempH,indx] = datasample(NumSpikesH,length(NumSpikesL),'Replace',false);
tempH2 = hyperVolH(indx);
figure(2); plot(hyperVolH,NumSpikesH,'o');
yticks(0:1:10); ylim([-1,10]);
xlabel('hyperpolarization voltage (in mV)');
ylabel('number of spikes in a burst');
title('Higher Order');


%% Plotting histograms for NumSpikes
% Subsampled to match number of data points in each subplot for higher and lower order
edgesPlot = -0.5:1:9.5;
higherOrderIndx = strcmp(inVitroBurstDataTable.cellOrder,'higherOrder');
higherOrderTable = inVitroBurstDataTable(higherOrderIndx,[8,14,19]);
f1 = figure(1); bFig1(1:12) = axes(f1); % lower order
lowerOrderIndx = strcmp(inVitroBurstDataTable.cellOrder,'lowerOrder');
lowerOrderTable = inVitroBurstDataTable(lowerOrderIndx,[8,14,19]);
f2 = figure(2); bFig2(1:12) = axes(f2); % higher order
for u = 1:6
    indxL = lowerOrderTable.HyperPolarBin == u;
    indxH = higherOrderTable.HyperPolarBin == u;
    numSpikesLower = lowerOrderTable.NumSpikes(indxL);
    numSpikesHigher = higherOrderTable.NumSpikes(indxH);
    if  length(numSpikesHigher) > length(numSpikesLower)
        numSpikesHigher = datasample(numSpikesHigher,length(numSpikesLower),'Replace',false);
    elseif length(numSpikesHigher) < length(numSpikesLower)
        numSpikesLower = datasample(numSpikesLower,length(numSpikesHigher),'Replace',false);
    end
    figure(1);
    bFig1(u) = subplot(3,2,u);
	histogram(numSpikesLower,edgesPlot);
    xlabel('num of spikes');
    xticks(0:10);
    title(strcat(num2str(edges(u)),{'mV '}, 'to', {' '},num2str(edges(u+1)),'mV'));
    figure(2);
    bFig2(u) = subplot(3,2,u);
    histogram(numSpikesHigher,edgesPlot);
    xlabel('num of spikes');
    xticks(0:10);
    title(strcat(num2str(edges(u)),{'mV '}, 'to', {' '},num2str(edges(u+1)),'mV'));
end
sgtitle(f1,'Lower Order\_Histogram of Number of Spikes');
sgtitle(f2,'Higher Order\_Histogram of Number of Spikes');
for u = 1:5
    linkaxes([bFig1(u+1),bFig1(u)],'xy');
    linkaxes([bFig2(u+1),bFig2(u)],'xy');
end

%% Time before burst after end of input pulse
%% Plotting  time before burst versus hyperpolarization voltage
lowerOrderIndx = strcmp(inVitroBurstDataTable.cellOrder,'lowerOrder');
lowerOrderTable = inVitroBurstDataTable(lowerOrderIndx,[8,14,15,20]);
rowsNumSpksNonZero = lowerOrderTable.NumSpikes ~= 0;
binsL = lowerOrderTable.HyperPolarBin;
binsL = binsL(rowsNumSpksNonZero);
hyperVolL = lowerOrderTable.hyperpolarVoltage;
hyperVolL = hyperVolL(rowsNumSpksNonZero);
timeTillFirstBurstL = lowerOrderTable.firstBurstTime;
timeTillFirstBurstL = cell2mat(timeTillFirstBurstL(rowsNumSpksNonZero));
figure(1); plot(hyperVolL,timeTillFirstBurstL,'o');
yticks(0:50:450); ylim([0,450]);
xlabel('hyperpolarization Voltage (in mV)');
ylabel('time until first burst after input stops (in ms)');
title('Lower Order');

higherOrderIndx = strcmp(inVitroBurstDataTable.cellOrder,'higherOrder');
higherOrderTable = inVitroBurstDataTable(higherOrderIndx,[8,14,15,20]);
binsH = higherOrderTable.HyperPolarBin;
binsH = binsH(rowsNumSpksNonZero);
hyperVolH = higherOrderTable.hyperpolarVoltage;
hyperVolH = hyperVolH(rowsNumSpksNonZero);
timeTillFirstBurstH = higherOrderTable.firstBurstTime;
timeTillFirstBurstH = cell2mat(timeTillFirstBurstH(rowsNumSpksNonZero));
figure(2); plot(hyperVolH,timeTillFirstBurstH,'o');
yticks(0:50:450); ylim([0,450]);
xlabel('hyperpolarization Voltage (in mV)');
ylabel('time until first burst after input stops (in ms)');
title('Higher Order');

% figure; hold on;
% plot(hyperVolL,log10(timeTillFirstBurstL),'r.','MarkerSize',12);
% plot(hyperVolH,log10(timeTillFirstBurstH),'b.','MarkerSize',12);
% xlabel('hyperpolarization Voltage (in mV)');
% ylabel('time until first burst after input stops (in ms)');
% legend('Lower Order','Higher Order');

%% Histogram of time before burst
edgesPlot = 20:40:400;
f1 = figure(1); bFig1(1:12) = axes(f1); % lower order
f2 = figure(2); bFig2(1:12) = axes(f2); % higher order
for u = 1:6
    indxL = binsL == u;
    indxH = binsH == u;
    timeLower = timeTillFirstBurstL(indxL);
    timeHigher = timeTillFirstBurstH(indxH);
    if  length(timeHigher) > length(timeLower)
        timeHigher = datasample(timeHigher,length(timeLower),'Replace',false);
    elseif length(timeHigher) < length(timeLower)
        timeLower = datasample(timeLower,length(timeHigher),'Replace',false);
    end
    figure(1);
    bFig1(u) = subplot(3,2,u);
	histogram(timeLower,edgesPlot);
    title(strcat(num2str(edges(u)),{'mV '}, 'to', {' '},num2str(edges(u+1)),'mV'));
    
    figure(2);
    bFig2(u) = subplot(3,2,u);
    histogram(timeHigher,edgesPlot);
    title(strcat(num2str(edges(u)),{'mV '}, 'to', {' '},num2str(edges(u+1)),'mV'));
end
sgtitle(f1,{'Lower Order\_Histogram of', 'time until first burst after input stops'});
sgtitle(f2,{'Higher Order\_Histogram of', 'time until first burst after input stops'});
for u = 1:5
    linkaxes([bFig1(u+1),bFig1(u)],'xy');
    linkaxes([bFig2(u+1),bFig2(u)],'xy');
end

%% First Intra-burst interval (IBI)
%% 
rowsNumSpksMoreThanOne = inVitroBurstDataTable.NumSpikes > 1;
firstIBI = inVitroBurstDataTable.firstIBI;
firstIBI = firstIBI(rowsNumSpksMoreThanOne);

%% Plot first IBI versus hyperpolarization voltage
temp = inVitroBurstDataTable(rowsNumSpksMoreThanOne,[8,16]);
plot(temp.hyperpolarVoltage,log10(temp.firstIBI),'o');
yticks(log10([2,3,4,5:5:25]));
yticklabels([2,3,4,5:5:25]);
xlabel('hyperpolarization voltage');
ylabel('first IBI');
title('first IBI vs. hyperpolarization voltage');

%% Plot histogram for all first IBI
edgesPlot = 0:0.25:25;
h = histogram(firstIBI,edgesPlot);
% E = h.BinEdges;
% y = h.BinCounts;
% xloc = E(1:end-1)+diff(E)/2;
% text(xloc, y+5, string(y));
% textStr = {'number of hyperpolarized traces',...
%     strcat('with IBI < 5ms is','{ }',string(sum(firstIBI<=5))),...
%     strcat('with IBI > 5ms is','{ }',string(sum(firstIBI>5)))};
% annotation('textbox',[.7 .5 .3 .3],'String',textStr,'FitBoxToText','on');
xlabel('first IBIs (in ms)');
xticks([2:1:25]);
title('Histogram of first Intra-burst interval');

%% Plot histogram of first IBI for lower and higher order
edgesPlot = 2:0.25:20; % 2:0.25:24.25;
lowerOrderIndx = strcmp(inVitroBurstDataTable.cellOrder,'lowerOrder');
firstIBILindx = lowerOrderIndx & rowsNumSpksMoreThanOne;
firstIBIL = inVitroBurstDataTable.firstIBI(firstIBILindx);
figure; s1 = subplot(2,1,1); histogram(firstIBIL,edgesPlot);
xticks([2:1:20]); xlabel('first IBI'); % xtickangle(90);
title('Lower Order');

higherOrderIndx = strcmp(inVitroBurstDataTable.cellOrder,'higherOrder');
firstIBIHindx = higherOrderIndx & rowsNumSpksMoreThanOne;
firstIBIH = inVitroBurstDataTable.firstIBI(firstIBIHindx);
tempfirstIBIH = datasample(firstIBIH,length(firstIBIL),'Replace',false);
s2 = subplot(2,1,2); histogram(tempfirstIBIH,edgesPlot); 
xticks([2:1:20]); xlabel('first IBI'); % xtickangle(90);
textStr = {'subsampled'};
annotation('textbox','String',textStr,'FitBoxToText','on');
title('Higher Order');
linkaxes([s1,s2],'xy');
sgtitle('Histogram of first intra-burst interval');

%% Plot histogram of first IBI for different sensory systems
visualIndx = strcmp(inVitroBurstDataTable.sensorySystem,'visual');
auditoryIndx = strcmp(inVitroBurstDataTable.sensorySystem,'auditory');
somatoIndx = strcmp(inVitroBurstDataTable.sensorySystem,'somato');
firstIBIVisual = inVitroBurstDataTable.firstIBI(visualIndx & rowsNumSpksMoreThanOne);
firstIBIAuditory = inVitroBurstDataTable.firstIBI(auditoryIndx & rowsNumSpksMoreThanOne);
firstIBISomato = inVitroBurstDataTable.firstIBI(somatoIndx & rowsNumSpksMoreThanOne);

sz = min([length(firstIBIVisual),length(firstIBIAuditory),length(firstIBISomato)]);
groups = [ones(sz,1); 2*ones(sz,1);3*ones(sz,1)];
boxPlotData = [datasample(firstIBIVisual,sz,'Replace',false)...
    ;datasample(firstIBIAuditory,sz,'Replace',false)...
    ;datasample(firstIBISomato,sz,'Replace',false)];
notBoxPlot(boxPlotData,groups,'style','sdline');
xlim([0.5,3.5]);
xticklabels({'Visual','Auditory','Somato'});
ylabel('first intra-burst interval');
title('boxplot of first IBI for 3 sensory regions');

%% Histogram of first intra-burst interval for different 
% hyperpolarization voltage bands and for lower and higher order
edgesPlot = 2:1:24;
lowerOrderIndx = strcmp(inVitroBurstDataTable.cellOrder,'lowerOrder');
binsL = inVitroBurstDataTable.HyperPolarBin(lowerOrderIndx & rowsNumSpksMoreThanOne);
higherOrderIndx = strcmp(inVitroBurstDataTable.cellOrder,'higherOrder');
binsH = inVitroBurstDataTable.HyperPolarBin(higherOrderIndx & rowsNumSpksMoreThanOne);

f1 = figure(1); bFig1(1:12) = axes(f1); % lower order
f2 = figure(2); bFig2(1:12) = axes(f2); % higher order
for u = 1:6
    indxL = binsL == u;
    indxH = binsH == u;
    ibiLower = firstIBIL(indxL);
    ibiHigher = firstIBIH(indxH);
    if  length(ibiHigher) > length(ibiLower)
        ibiHigher = datasample(ibiHigher,length(ibiLower),'Replace',false);
    elseif length(ibiHigher) < length(ibiLower)
        ibiLower = datasample(ibiLower,length(ibiHigher),'Replace',false);
    end
    figure(1);
    bFig1(u) = subplot(3,2,u);
	histogram(ibiLower,edgesPlot); ylim([0,35]); xticks(2:2:24);
    title(strcat(num2str(edges(u)),{'mV '}, 'to', {' '},num2str(edges(u+1)),'mV'));
    
    figure(2);
    bFig2(u) = subplot(3,2,u);
    histogram(ibiHigher,edgesPlot); ylim([0,35]); xticks(2:2:24);
    title(strcat(num2str(edges(u)),{'mV '}, 'to', {' '},num2str(edges(u+1)),'mV'));
end
sgtitle(f1,{'Lower Order\_Histogram of', 'first intra-burst interval'});
sgtitle(f2,{'Higher Order\_Histogram of', 'first intra-burst interval'});
for u = 1:5
    linkaxes([bFig1(u+1),bFig1(u)],'xy');
    linkaxes([bFig2(u+1),bFig2(u)],'xy');
end

%% 2nd January 2020
%% Plot boxplot of first IBI for different brain sensory regions and orders
visualIndx = strcmp(inVitroBurstDataTable.sensorySystem,'visual');
auditoryIndx = strcmp(inVitroBurstDataTable.sensorySystem,'auditory');
somatoIndx = strcmp(inVitroBurstDataTable.sensorySystem,'somato');
lowerOrderIndx = strcmp(inVitroBurstDataTable.cellOrder,'lowerOrder');
higherOrderIndx = strcmp(inVitroBurstDataTable.cellOrder,'higherOrder');
fibiVL = inVitroBurstDataTable.firstIBI(visualIndx & lowerOrderIndx & rowsNumSpksMoreThanOne);
fibiVH = inVitroBurstDataTable.firstIBI(visualIndx & higherOrderIndx & rowsNumSpksMoreThanOne);
fibiAL = inVitroBurstDataTable.firstIBI(auditoryIndx & lowerOrderIndx & rowsNumSpksMoreThanOne);
fibiAH = inVitroBurstDataTable.firstIBI(auditoryIndx & higherOrderIndx & rowsNumSpksMoreThanOne);
fibiSL = inVitroBurstDataTable.firstIBI(somatoIndx & lowerOrderIndx & rowsNumSpksMoreThanOne);
fibiSH = inVitroBurstDataTable.firstIBI(somatoIndx & higherOrderIndx & rowsNumSpksMoreThanOne);

sz = min([length(fibiVL),length(fibiVH),length(fibiAL),...
    length(fibiAH),length(fibiSL),length(fibiSH)]);
groups = [ones(sz,1); 2*ones(sz,1);3*ones(sz,1);4*ones(sz,1);5*ones(sz,1);6*ones(sz,1)];
fibiVL = datasample(fibiVL,sz,'Replace',false);
fibiVH = datasample(fibiVH,sz,'Replace',false);
fibiAL = datasample(fibiAL,sz,'Replace',false);
fibiAH = datasample(fibiAH,sz,'Replace',false);
fibiSL = datasample(fibiSL,sz,'Replace',false);
fibiSH = datasample(fibiSH,sz,'Replace',false);
boxPlotData = [fibiVL;fibiVH;fibiAL;fibiAH;fibiSL;fibiSH];
boxplot(boxPlotData,groups,'Colors',rgb('Blue'));
hold on;
xx = ones(1,sz).*(1+(rand(1,sz)-0.5)/2.2);
scatter(xx, fibiVL,'o','filled','MarkerEdgeColor','k','MarkerFaceColor', rgb('DarkGray'));
xx = ones(1,sz).*(2+(rand(1,sz)-0.5)/2.2);
scatter(xx, fibiVH,'o','filled','MarkerEdgeColor','k','MarkerFaceColor', rgb('DarkGray'));
xx = ones(1,sz).*(3+(rand(1,sz)-0.5)/2.2);
scatter(xx, fibiAL,'o','filled','MarkerEdgeColor','k','MarkerFaceColor', rgb('DarkGray'));
xx = ones(1,sz).*(4+(rand(1,sz)-0.5)/2.2);
scatter(xx, fibiAH,'o','filled','MarkerEdgeColor','k','MarkerFaceColor', rgb('DarkGray'));
xx = ones(1,sz).*(5+(rand(1,sz)-0.5)/2.2);
scatter(xx, fibiSL,'o','filled','MarkerEdgeColor','k','MarkerFaceColor', rgb('DarkGray'));
xx = ones(1,sz).*(6+(rand(1,sz)-0.5)/2.2); 
scatter(xx, fibiSH,'o','filled','MarkerEdgeColor','k','MarkerFaceColor', rgb('DarkGray'));
xticklabels({'Visual_Lower','Visual_Higher','Auditory_Lower'...
    ,'Auditory_Higher','Somato_Lower','Somato_Higher'});
ylabel('first intra-burst interval');
title('boxplot of first IBI for 3 sensory regions and 3 orders');
xlim([0.5,6.5]);

%% Plot histogram of first IBI for different brain sensory regions and orders
edgesPlot = 2:0.25:20;
f1 = figure(1); bFig1(1:12) = axes(f1); % lower order
bFig1(1) = subplot(3,1,1); histogram(fibiVL,edgesPlot); ylim([0,15]); xlabel('Visual\_Lower');
bFig1(2) = subplot(3,1,2); histogram(fibiAL,edgesPlot); xlabel('Auditory\_Lower');
bFig1(3) = subplot(3,1,3); histogram(fibiSL,edgesPlot); xlabel('Somato\_Lower');
sgtitle(f1,{'Lower Order\_Histogram of', 'first intra-burst interval'});
linkaxes([bFig1(2),bFig1(1)],'xy');
linkaxes([bFig1(3),bFig1(2)],'xy');

f2 = figure(2); bFig2(1:12) = axes(f2); % lower order
bFig2(1) = subplot(3,1,1); histogram(fibiVH,edgesPlot); xlabel('Visual\_Higher');
bFig2(2) = subplot(3,1,2); histogram(fibiAH,edgesPlot); xlabel('Auditory\_Higher');
bFig2(3) = subplot(3,1,3); histogram(fibiSH,edgesPlot); xlabel('Somato\_Higher');
sgtitle(f2,{'Higher Order\_Histogram of', 'first intra-burst interval'});
linkaxes([bFig2(2),bFig2(1)],'xy');
linkaxes([bFig2(3),bFig2(2)],'xy');

%% all IBIs
%% Average frequency of burst in lower and higher orders
ibisL = inVitroBurstDataTable.allIBIs(lowerOrderIndx & rowsNumSpksMoreThanOne);
ibisH = inVitroBurstDataTable.allIBIs(higherOrderIndx & rowsNumSpksMoreThanOne);

avgFreqL = zeros(size(ibisL));
for q = 1:size(ibisL,1)
    temp = ibisL{q,1};
    temp = temp./1000; % converting millsecs to secs
    f = 1./temp; 
    avgFreqL(q) = mean(f);
end

avgFreqH = zeros(size(ibisH));
for q = 1:size(ibisH,1)
    temp = ibisH{q,1};
	temp = temp./1000; % converting millsecs to secs
    f = 1./temp;
    avgFreqH(q) = mean(f);
end
edgesHist = 40:3:480;
avgFreqH_SubSampled = datasample(avgFreqH,length(avgFreqL),'Replace',false); 
f1 = figure(1); bFig1(1:2) = axes(f1); % lower order
bFig1(1) = subplot(2,1,1); histogram(avgFreqL,edgesHist); xlabel('Lower Order');
bFig1(2) = subplot(2,1,2); histogram(avgFreqH_SubSampled,edgesHist); xlabel('Higher Order');
linkaxes([bFig1(2),bFig1(1)],'xy');
sgtitle(f1,'Histogram of average frequency of bursts');

%% Average frequency of burst in different nuclie
nuclie = inVitroBurstDataTable.brainRegion(rowsNumSpksMoreThanOne);
ibis = inVitroBurstDataTable.allIBIs(rowsNumSpksMoreThanOne);
avgFreq = zeros(size(ibis));
for q = 1:size(ibis,1)
    temp = ibis{q,1};
	temp = temp./1000; % converting millsecs to secs
    f = 1./temp;
    avgFreq(q) = mean(f);
end
RTlgn = avgFreq(strcmp(nuclie,'LGN'));
RTdlgn = avgFreq(strcmp(nuclie,'dLGN'));
RTlp = avgFreq(strcmp(nuclie,'LP'));
RTvb = avgFreq(strcmp(nuclie,'VB'));
RTvp = avgFreq(strcmp(nuclie,'VP'));
RTpom = avgFreq(strcmp(nuclie,'POm'));
RTdmgb = avgFreq(strcmp(nuclie,'dMGB'));
RTvmgb = avgFreq(strcmp(nuclie,'vMGB'));

edgesHist = 40:3:480;
f1 = figure(1); bFig1(1:8) = axes(f1); % lower order
bFig1(1) = subplot(4,2,1); histogram(RTlgn,edgesHist); xlabel('LGN Visual Lower');
bFig1(2) = subplot(4,2,2); histogram(RTdlgn,edgesHist); xlabel('dLGN Visual Lower');
bFig1(3) = subplot(4,2,3); histogram(RTvmgb,edgesHist); xlabel('vMGB Auditory Lower');
bFig1(4) = subplot(4,2,4); histogram(RTvb,edgesHist); xlabel('VB Somato Lower');
bFig1(5) = subplot(4,2,5); histogram(RTvp,edgesHist); xlabel('VP Somato Lower');

bFig1(6) = subplot(4,2,6); histogram(RTlp,edgesHist); xlabel('LP Visual Higher');
bFig1(7) = subplot(4,2,7); histogram(RTdmgb,edgesHist); xlabel('dMGB Auditory Higher');
bFig1(8) = subplot(4,2,8); histogram(RTpom,edgesHist); xlabel('POm Somato Higher');
linkaxes([bFig1(1),bFig1(2),bFig1(3),bFig1(4),bFig1(5),bFig1(6),bFig1(7),bFig1(8)],'xy');
sgtitle(f1,{'Histogram of average frequency of bursts','for different nuclie'});

%% Average IBI in lower and higher orders
avgTimeL = zeros(size(ibisL));
for q = 1:size(ibisL,1)
    temp = ibisL{q,1};
    avgTimeL(q) = mean(temp);
end

avgTimeH = zeros(size(ibisH));
for q = 1:size(ibisH,1)
    temp = ibisH{q,1};
    avgTimeH(q) = mean(temp);
end

edgesHist = 2:0.25:20;
avgTimeH_SubSampled = datasample(avgTimeH,length(avgFreqL),'Replace',false); 
f1 = figure(1); bFig1(1:2) = axes(f1); % lower order
bFig1(1) = subplot(2,1,1); histogram(avgTimeL,edgesHist); 
xticks([2:1:30]); xlabel('Lower Order');
bFig1(2) = subplot(2,1,2); histogram(avgTimeH_SubSampled,edgesHist); 
xticks([2:1:30]); xlabel('Higher Order');
linkaxes([bFig1(2),bFig1(1)],'xy');
sgtitle(f1,'Histogram of average IBI of bursts (in ms)');

%% all IBIs versus each burst spike number
ibis = inVitroBurstDataTable.allIBIs(rowsNumSpksMoreThanOne);
figure; hold on;
for p = 1:size(ibis,1)
    temp = ibis{p,1};
    for q = 1:size(temp,1)-1
        if temp(q+1) >= temp(q)
            plot([q, q+1], [temp(q), temp(q+1)],'b');
        else
            plot([q, q+1], [temp(q), temp(q+1)],'r');
        end
    end
%     plot(temp);
end
xlabel('IBI Number in a burst');
ylabel('IBI (in ms)');

%% all IBIs versus each burst spike number separated by orders
ibisL = inVitroBurstDataTable.allIBIs(lowerOrderIndx & rowsNumSpksMoreThanOne);
ibisH = inVitroBurstDataTable.allIBIs(higherOrderIndx & rowsNumSpksMoreThanOne);
ibisH_subsampled = datasample(ibisH,length(ibisL),'Replace',false);
f1 = figure(1); bFig1(1:2) = axes(f1); % lower order
bFig1(1) = subplot(2,1,1); hold on;
for p = 1:size(ibisL,1)
    temp = ibisL{p,1};
    for q = 1:size(temp,1)-1
        if temp(q+1) >= temp(q)
            plot([q, q+1], [temp(q), temp(q+1)],'b');
        else
            plot([q, q+1], [temp(q), temp(q+1)],'r');
        end
    end
end
hold off;
bFig1(2) = subplot(2,1,2); hold on;
for p = 1:size(ibisH_subsampled,1)
    temp = ibisH_subsampled{p,1};
    for q = 1:size(temp,1)-1
        if temp(q+1) >= temp(q)
            plot([q, q+1], [temp(q), temp(q+1)],'b');
        else
            plot([q, q+1], [temp(q), temp(q+1)],'r');
        end
    end
end
linkaxes([bFig1(2),bFig1(1)],'xy');
xticks([1:1:6]);
xlabel(bFig1(1),'Lower\_Order');
xlabel(bFig1(2),'Higher\_Order');
sgtitle(f1,'Histogram of average IBI of bursts (in ms)');
xlabel('IBI Number in a burst');
ylabel('IBI (in ms)');
ylim([0,20]);

%% Duration of burst
pkTimesL = inVitroBurstDataTable.peakTimes(lowerOrderIndx &rowsNumSpksMoreThanOne);
burstDurL = zeros(size(pkTimesL));
for p = 1:size(burstDurL,1)
    temp = pkTimesL{p,1};
    burstDurL(p) = temp(end) - temp(1);
end
pkTimesH = inVitroBurstDataTable.peakTimes(higherOrderIndx &rowsNumSpksMoreThanOne);
burstDurH = zeros(size(pkTimesH));
for p = 1:size(burstDurH,1)
    temp = pkTimesH{p,1};
    burstDurH(p) = temp(end) - temp(1);
end
edgesPlot = 2:2:100;
% edgesPlot = 2:1:101;
% burstDurH_SubSampled = datasample(burstDurH,length(burstDurL),'Replace',false); 
f1 = figure(1); bFig1(1:2) = axes(f1); % lower order
bFig1(1) = subplot(2,1,1); histogram(burstDurL,edgesPlot); xlabel('Lower\_Order'); 
bFig1(2) = subplot(2,1,2); histogram(burstDurH,edgesPlot); xlabel('Higher\_Order');
sgtitle(f1,{'Histogram of duration of burst(in ms)'});
linkaxes([bFig1(2),bFig1(1)],'xy');

%% Relative change in frequency of burst of 2nd IBI relative to 1st IBI
rowsNumSpksMoreThanTwo = inVitroBurstDataTable.NumSpikes > 2;
ibisL = inVitroBurstDataTable.allIBIs(lowerOrderIndx & rowsNumSpksMoreThanTwo);
freqChangeRatioL = zeros(size(ibisL));
for p = 1:size(ibisL,1)
    temp = ibisL{p,1};
    temp = temp./1000;
    f1 = 1/temp(1);
    f2 = 1/temp(2);
    freqChangeRatioL(p) = (f2-f1)/f1;
end
ibisH = inVitroBurstDataTable.allIBIs(higherOrderIndx & rowsNumSpksMoreThanTwo);
freqChangeRatioH = zeros(size(ibisH));
for p = 1:size(ibisH,1)
    temp = ibisH{p,1};
    temp = temp./1000;
    f1 = 1/temp(1);
    f2 = 1/temp(2);
    freqChangeRatioH(p) = (f2-f1)/f1;
end
% Histogram
edgesPlot = -0.8:0.025:0.8;
freqChangeRatioH_SubSampled = datasample(freqChangeRatioH,length(freqChangeRatioL),'Replace',false); 
f1 = figure(1); bFig1(1:2) = axes(f1); % lower order
bFig1(1) = subplot(2,1,1); histogram(freqChangeRatioL,edgesPlot); xlabel('Lower\_Order');
bFig1(2) = subplot(2,1,2); histogram(freqChangeRatioH_SubSampled,edgesPlot); xlabel('Higher\_Order');
sgtitle(f1,{'Histogram of relative change in frequency between 1st and 2nd IBIs (in ms)'});
linkaxes([bFig1(2),bFig1(1)],'xy');

%% change in IBIs versus each burst spike number
ibis = inVitroBurstDataTable.allIBIs(rowsNumSpksMoreThanOne);
figure; hold on;
ibis = ibisH;
for p = 1:size(ibis,1)
    temp = ibis{p,1};
    temp = diff(temp);
    for q = 1:size(temp,1)-1
        if temp(q+1) >= temp(q)
            plot([q, q+1], [temp(q), temp(q+1)],'b');
        else
            plot([q, q+1], [temp(q), temp(q+1)],'r');
        end
    end
%     plot(temp);
end
xlabel('IBI Number in a burst');
ylabel('IBI change (in ms)');

%% Calcium rise
%% difference in time between burst calcium rise slope change and 
% hyperpolarizing input stop time for lower and higher order cells
% histogram
CaRiseTimeL = cell2mat(inVitroBurstDataTable.Ca2RiseSlopeChangeTime(lowerOrderIndx & rowsNumSpksMoreThanOne));
CaRiseTimeH = cell2mat(inVitroBurstDataTable.Ca2RiseSlopeChangeTime(higherOrderIndx & rowsNumSpksMoreThanOne));
HyperStopL = inVitroBurstDataTable.hyperpolarStopTime(lowerOrderIndx & rowsNumSpksMoreThanOne);
HyperStopH = inVitroBurstDataTable.hyperpolarStopTime(higherOrderIndx & rowsNumSpksMoreThanOne);
CaRiseChangeL =  CaRiseTimeL - HyperStopL;
CaRiseChangeH =  CaRiseTimeH - HyperStopH;

edgesPlot = 5:5:300;
CaRiseChangeH_SubSampled = datasample(CaRiseChangeH,length(avgFreqL),'Replace',false); 
f1 = figure(1); bFig1(1:2) = axes(f1); % lower order
bFig1(1) = subplot(2,1,1); histogram(CaRiseChangeL,edgesPlot); xlabel('Lower\_Order');
bFig1(2) = subplot(2,1,2); histogram(CaRiseChangeH_SubSampled,edgesPlot); xlabel('Higher\_Order');
sgtitle(f1,{'Histogram of calcium rise time (in ms)'});
linkaxes([bFig1(2),bFig1(1)],'xy');

% boxplot
sz = 191;
groups = [ones(sz,1); 2*ones(sz,1)];
boxPlotData = [CaRiseChangeL;CaRiseChangeH_SubSampled];
boxplot(boxPlotData,groups,'Colors',rgb('Blue'),'Notch','on');
hold on;
xx = ones(1,sz).*(1+(rand(1,sz)-0.5)/4);
scatter(xx, CaRiseChangeL,'o','filled','MarkerEdgeColor','k','MarkerFaceColor', rgb('DarkGray'));
xx = ones(1,sz).*(2+(rand(1,sz)-0.5)/4);
scatter(xx, CaRiseChangeH_SubSampled,'o','filled','MarkerEdgeColor','k','MarkerFaceColor', rgb('DarkGray'));
xticklabels({'Lower_Order', 'Higher_Order'});
ylabel('calcium rise time (in ms)');
title('Boxplot of calcium rise time for lower and higher order cells(in ms)');

%% difference in time between first burst spike and calcium slope change time
% histogram
CaRiseTimeL = cell2mat(inVitroBurstDataTable.Ca2RiseSlopeChangeTime(lowerOrderIndx & rowsNumSpksMoreThanOne));
CaRiseTimeH = cell2mat(inVitroBurstDataTable.Ca2RiseSlopeChangeTime(higherOrderIndx & rowsNumSpksMoreThanOne));
BSpksTimeL = inVitroBurstDataTable.peakTimes(lowerOrderIndx & rowsNumSpksMoreThanOne);
BSpksTimeH = inVitroBurstDataTable.peakTimes(higherOrderIndx & rowsNumSpksMoreThanOne);
firstBSpkTimeL = zeros(size(BSpksTimeL));
for p = 1:size(BSpksTimeL,1)
    temp = BSpksTimeL{p,1};
    firstBSpkTimeL(p) = temp(1);
end    
firstBSpkTimeH = zeros(size(BSpksTimeH));
for p = 1:size(BSpksTimeH,1)
    temp = BSpksTimeH{p,1};
    firstBSpkTimeH(p) = temp(1);
end
NaRiseTimeL = firstBSpkTimeL - CaRiseTimeL;
NaRiseTimeH = firstBSpkTimeH - CaRiseTimeH;
edgesPlot = 8:2:124;
NaRiseTimeH_SubSampled = datasample(NaRiseTimeH,length(NaRiseTimeL),'Replace',false); 
f1 = figure(1); bFig1(1:2) = axes(f1); % lower order
bFig1(1) = subplot(2,1,1); histogram(NaRiseTimeL,edgesPlot); xlabel('Lower\_Order');
bFig1(2) = subplot(2,1,2); histogram(NaRiseTimeH_SubSampled,edgesPlot); xlabel('Higher\_Order');
sgtitle(f1,{'Histogram of time between Ca slope change and first spike (in ms)'});
linkaxes([bFig1(2),bFig1(1)],'xy');

% boxplot
sz = 191;
groups = [ones(sz,1); 2*ones(sz,1)];
boxPlotData = [NaRiseTimeL;NaRiseTimeH_SubSampled];
boxplot(boxPlotData,groups,'Colors',rgb('Blue'),'Notch','on');
hold on;
xx = ones(1,sz).*(1+(rand(1,sz)-0.5)/4);
scatter(xx, NaRiseTimeL,'o','filled','MarkerEdgeColor','k','MarkerFaceColor', rgb('DarkGray'));
xx = ones(1,sz).*(2+(rand(1,sz)-0.5)/4);
scatter(xx, NaRiseTimeH_SubSampled,'o','filled','MarkerEdgeColor','k','MarkerFaceColor', rgb('DarkGray'));
xticklabels({'Lower_Order', 'Higher_Order'});
ylabel('calcium rise time (in ms)');
title('Boxplot of time from Ca slope change to first burst spike for lower and higher order cells(in ms)');

%% difference in time between first burst spike and first burst spike time
hyperVolLimitIndx = inVitroBurstDataTable.hyperpolarVoltage > -90;

HyperStopL = inVitroBurstDataTable.hyperpolarStopTime(lowerOrderIndx & rowsNumSpksMoreThanOne & hyperVolLimitIndx);
HyperStopH = inVitroBurstDataTable.hyperpolarStopTime(higherOrderIndx & rowsNumSpksMoreThanOne & hyperVolLimitIndx);
BSpksTimeL = inVitroBurstDataTable.peakTimes(lowerOrderIndx & rowsNumSpksMoreThanOne & hyperVolLimitIndx);
BSpksTimeH = inVitroBurstDataTable.peakTimes(higherOrderIndx & rowsNumSpksMoreThanOne & hyperVolLimitIndx);
firstBSpkTimeL = zeros(size(BSpksTimeL));
for p = 1:size(BSpksTimeL,1)
    temp = BSpksTimeL{p,1};
    firstBSpkTimeL(p) = temp(1);
end    
firstBSpkTimeH = zeros(size(BSpksTimeH));
for p = 1:size(BSpksTimeH,1)
    temp = BSpksTimeH{p,1};
    firstBSpkTimeH(p) = temp(1);
end
riseTimeL = firstBSpkTimeL - HyperStopL;
riseTimeH = firstBSpkTimeH - HyperStopH;
edgesPlot = 20:2.5:300;
% riseTimeH_SubSampled = datasample(riseTimeH,length(riseTimeL),'Replace',false); 
f1 = figure(1); bFig1(1:2) = axes(f1); % lower order
bFig1(1) = subplot(2,1,1); histogram(riseTimeL,edgesPlot); xlabel('Lower\_Order');
bFig1(2) = subplot(2,1,2); histogram(riseTimeH_SubSampled,edgesPlot); xlabel('Higher\_Order');
sgtitle(f1,{'Histogram of total time between hyperpolarization stop and first spike (in ms)'});
linkaxes([bFig1(2),bFig1(1)],'xy');

figure; hold on;
histogram(riseTimeL,edgesPlot,'Normalization','probability');
histogram(riseTimeH,edgesPlot,'Normalization','probability');



% boxplot
sz = 191;
groups = [ones(sz,1); 2*ones(sz,1)];
boxPlotData = [riseTimeL;riseTimeH_SubSampled];
boxplot(boxPlotData,groups,'Colors',rgb('Blue'),'Notch','on');
hold on;
xx = ones(1,sz).*(1+(rand(1,sz)-0.5)/4);
scatter(xx, riseTimeL,'o','filled','MarkerEdgeColor','k','MarkerFaceColor', rgb('DarkGray'));
xx = ones(1,sz).*(2+(rand(1,sz)-0.5)/4);
scatter(xx, riseTimeH_SubSampled,'o','filled','MarkerEdgeColor','k','MarkerFaceColor', rgb('DarkGray'));
xticklabels({'Lower_Order', 'Higher_Order'});
ylabel('calcium rise time (in ms)');
title('Boxplot of total time between hyperpolarization stop and first spike for lower and higher order cells(in ms)');

%% calcium rise time for 3 sensory regions for both lower and higher orders
visualIndx = strcmp(inVitroBurstDataTable.sensorySystem,'visual');
auditoryIndx = strcmp(inVitroBurstDataTable.sensorySystem,'auditory');
somatoIndx = strcmp(inVitroBurstDataTable.sensorySystem,'somato');
lowerOrderIndx = strcmp(inVitroBurstDataTable.cellOrder,'lowerOrder');
higherOrderIndx = strcmp(inVitroBurstDataTable.cellOrder,'higherOrder');

caRTVL = cell2mat(inVitroBurstDataTable.Ca2RiseSlopeChangeTime...
    (visualIndx & lowerOrderIndx & rowsNumSpksMoreThanOne)) - ...
    inVitroBurstDataTable.hyperpolarStopTime(visualIndx & lowerOrderIndx & rowsNumSpksMoreThanOne);
caRTVH = cell2mat(inVitroBurstDataTable.Ca2RiseSlopeChangeTime...
    (visualIndx & higherOrderIndx & rowsNumSpksMoreThanOne)) - ...
    inVitroBurstDataTable.hyperpolarStopTime(visualIndx & higherOrderIndx & rowsNumSpksMoreThanOne);
caRTAL = cell2mat(inVitroBurstDataTable.Ca2RiseSlopeChangeTime...
    (auditoryIndx & lowerOrderIndx & rowsNumSpksMoreThanOne)) - ...
    inVitroBurstDataTable.hyperpolarStopTime(auditoryIndx & lowerOrderIndx & rowsNumSpksMoreThanOne);
caRTAH = cell2mat(inVitroBurstDataTable.Ca2RiseSlopeChangeTime...
    (auditoryIndx & higherOrderIndx & rowsNumSpksMoreThanOne)) - ...
    inVitroBurstDataTable.hyperpolarStopTime(auditoryIndx & higherOrderIndx & rowsNumSpksMoreThanOne);
caRTSL = cell2mat(inVitroBurstDataTable.Ca2RiseSlopeChangeTime...
    (somatoIndx & lowerOrderIndx & rowsNumSpksMoreThanOne)) - ...
    inVitroBurstDataTable.hyperpolarStopTime(somatoIndx & lowerOrderIndx & rowsNumSpksMoreThanOne);
caRTSH = cell2mat(inVitroBurstDataTable.Ca2RiseSlopeChangeTime...
    (somatoIndx & higherOrderIndx & rowsNumSpksMoreThanOne)) - ...
    inVitroBurstDataTable.hyperpolarStopTime(somatoIndx & higherOrderIndx & rowsNumSpksMoreThanOne);

sz = min([length(caRTVL),length(caRTVH),length(caRTAL),...
    length(caRTAH),length(caRTSL),length(caRTSH)]);
groups = [ones(sz,1); 2*ones(sz,1);3*ones(sz,1);4*ones(sz,1);5*ones(sz,1);6*ones(sz,1)];
caRTVL = datasample(caRTVL,sz,'Replace',false);
caRTVH = datasample(caRTVH,sz,'Replace',false);
caRTAL = datasample(caRTAL,sz,'Replace',false);
caRTAH = datasample(caRTAH,sz,'Replace',false);
caRTSL = datasample(caRTSL,sz,'Replace',false);
caRTSH = datasample(caRTSH,sz,'Replace',false);
boxPlotData = [caRTVL;caRTVH;caRTAL;caRTAH;caRTSL;caRTSH];
boxplot(boxPlotData,groups,'Colors',rgb('Blue'),'Notch','on');
hold on;
xx = ones(1,sz).*(1+(rand(1,sz)-0.5)/2.2);
scatter(xx, caRTVL,'o','filled','MarkerEdgeColor','k','MarkerFaceColor', rgb('DarkGray'));
xx = ones(1,sz).*(2+(rand(1,sz)-0.5)/2.2);
scatter(xx, caRTVH,'o','filled','MarkerEdgeColor','k','MarkerFaceColor', rgb('DarkGray'));
xx = ones(1,sz).*(3+(rand(1,sz)-0.5)/2.2);
scatter(xx, caRTAL,'o','filled','MarkerEdgeColor','k','MarkerFaceColor', rgb('DarkGray'));
xx = ones(1,sz).*(4+(rand(1,sz)-0.5)/2.2);
scatter(xx, caRTAH,'o','filled','MarkerEdgeColor','k','MarkerFaceColor', rgb('DarkGray'));
xx = ones(1,sz).*(5+(rand(1,sz)-0.5)/2.2);
scatter(xx, caRTSL,'o','filled','MarkerEdgeColor','k','MarkerFaceColor', rgb('DarkGray'));
xx = ones(1,sz).*(6+(rand(1,sz)-0.5)/2.2); 
scatter(xx, caRTSH,'o','filled','MarkerEdgeColor','k','MarkerFaceColor', rgb('DarkGray'));
xticklabels({'Visual_Lower','Visual_Higher','Auditory_Lower'...
    ,'Auditory_Higher','Somato_Lower','Somato_Higher'});
ylabel('calcium rise time (in ms)');
title('boxplot of calcium rise time for 3 sensory regions and 3 orders');
xlim([0.5,6.5]);

%% difference in time between end of input current and first burst spike time
% for 3 sensory regions for both lower and higher orders
visualIndx = strcmp(inVitroBurstDataTable.sensorySystem,'visual');
auditoryIndx = strcmp(inVitroBurstDataTable.sensorySystem,'auditory');
somatoIndx = strcmp(inVitroBurstDataTable.sensorySystem,'somato');
lowerOrderIndx = strcmp(inVitroBurstDataTable.cellOrder,'lowerOrder');
higherOrderIndx = strcmp(inVitroBurstDataTable.cellOrder,'higherOrder');

totalRTVL = inVitroBurstDataTable.firstPeakTime(visualIndx & lowerOrderIndx & rowsNumSpksMoreThanOne) - ...
    inVitroBurstDataTable.hyperpolarStopTime(visualIndx & lowerOrderIndx & rowsNumSpksMoreThanOne);
totalRTVH =inVitroBurstDataTable.firstPeakTime(visualIndx & higherOrderIndx & rowsNumSpksMoreThanOne) - ...
    inVitroBurstDataTable.hyperpolarStopTime(visualIndx & higherOrderIndx & rowsNumSpksMoreThanOne);
totalRTAL = inVitroBurstDataTable.firstPeakTime(auditoryIndx & lowerOrderIndx & rowsNumSpksMoreThanOne) - ...
    inVitroBurstDataTable.hyperpolarStopTime(auditoryIndx & lowerOrderIndx & rowsNumSpksMoreThanOne);
totalRTAH = inVitroBurstDataTable.firstPeakTime(auditoryIndx & higherOrderIndx & rowsNumSpksMoreThanOne)- ...
    inVitroBurstDataTable.hyperpolarStopTime(auditoryIndx & higherOrderIndx & rowsNumSpksMoreThanOne);
totalRTSL = inVitroBurstDataTable.firstPeakTime(somatoIndx & lowerOrderIndx & rowsNumSpksMoreThanOne) - ...
    inVitroBurstDataTable.hyperpolarStopTime(somatoIndx & lowerOrderIndx & rowsNumSpksMoreThanOne);
totalRTSH = inVitroBurstDataTable.firstPeakTime(somatoIndx & higherOrderIndx & rowsNumSpksMoreThanOne)- ...
    inVitroBurstDataTable.hyperpolarStopTime(somatoIndx & higherOrderIndx & rowsNumSpksMoreThanOne);

sz = min([length(totalRTVL),length(totalRTVH),length(totalRTAL),...
    length(totalRTAH),length(totalRTSL),length(totalRTSH)]);
groups = [ones(sz,1); 2*ones(sz,1);3*ones(sz,1);4*ones(sz,1);5*ones(sz,1);6*ones(sz,1)];
totalRTVL = datasample(totalRTVL,sz,'Replace',false);
totalRTVH = datasample(totalRTVH,sz,'Replace',false);
totalRTAL = datasample(totalRTAL,sz,'Replace',false);
totalRTAH = datasample(totalRTAH,sz,'Replace',false);
totalRTSL = datasample(totalRTSL,sz,'Replace',false);
totalRTSH = datasample(totalRTSH,sz,'Replace',false);
boxPlotData = [totalRTVL;totalRTVH;totalRTAL;totalRTAH;totalRTSL;totalRTSH];
boxplot(boxPlotData,groups,'Colors',rgb('Blue'),'Notch','on');
hold on;
xx = ones(1,sz).*(1+(rand(1,sz)-0.5)/2.2);
scatter(xx, totalRTVL,'o','filled','MarkerEdgeColor','k','MarkerFaceColor', rgb('DarkGray'));
xx = ones(1,sz).*(2+(rand(1,sz)-0.5)/2.2);
scatter(xx, totalRTVH,'o','filled','MarkerEdgeColor','k','MarkerFaceColor', rgb('DarkGray'));
xx = ones(1,sz).*(3+(rand(1,sz)-0.5)/2.2);
scatter(xx, totalRTAL,'o','filled','MarkerEdgeColor','k','MarkerFaceColor', rgb('DarkGray'));
xx = ones(1,sz).*(4+(rand(1,sz)-0.5)/2.2);
scatter(xx, totalRTAH,'o','filled','MarkerEdgeColor','k','MarkerFaceColor', rgb('DarkGray'));
xx = ones(1,sz).*(5+(rand(1,sz)-0.5)/2.2);
scatter(xx, totalRTSL,'o','filled','MarkerEdgeColor','k','MarkerFaceColor', rgb('DarkGray'));
xx = ones(1,sz).*(6+(rand(1,sz)-0.5)/2.2); 
scatter(xx, totalRTSH,'o','filled','MarkerEdgeColor','k','MarkerFaceColor', rgb('DarkGray'));
xticklabels({'Visual_Lower','Visual_Higher','Auditory_Lower'...
    ,'Auditory_Higher','Somato_Lower','Somato_Higher'});
ylabel('total rise time (in ms)');
title({'boxplot of time from hyperpolarization stop to first burst spike',...
    'for 3 sensory regions and 3 orders'});
xlim([0.5,6.5]);

%% boxplot difference in time between end of input current and first burst spike time
% for 9 brain regions
tempT = table2cell(inVitroBurstDataTable(rowsNumSpksMoreThanOne,[4,10,13]));
tempT2 = [tempT(:,1),num2cell(cell2mat(tempT(:,3))-cell2mat(tempT(:,2)))];
RTlgn = cell2mat(tempT2(strcmp(tempT2(:,1),'LGN'),2));
RTdlgn = cell2mat(tempT2(strcmp(tempT2(:,1),'dLGN'),2));
RTlp = cell2mat(tempT2(strcmp(tempT2(:,1),'LP'),2));
RTvb = cell2mat(tempT2(strcmp(tempT2(:,1),'VB'),2));
RTvp = cell2mat(tempT2(strcmp(tempT2(:,1),'VP'),2));
RTpom = cell2mat(tempT2(strcmp(tempT2(:,1),'POm'),2));
RTdmgb = cell2mat(tempT2(strcmp(tempT2(:,1),'dMGB'),2));
RTvmgb = cell2mat(tempT2(strcmp(tempT2(:,1),'vMGB'),2));
groups = [ones(length(RTlgn),1); 2*ones(length(RTdlgn),1);...
    3*ones(length(RTlp),1); 4*ones(length(RTvb),1);...
    5*ones(length(RTvp),1);6*ones(length(RTpom),1);...
    7*ones(length(RTvmgb),1); 8*ones(length(RTdmgb),1)];
boxPlotData = [RTlgn;RTdlgn;RTlp;RTvb;RTvp;RTpom;RTvmgb;RTdmgb];
boxplot(boxPlotData,groups,'Colors',rgb('Blue'),'Notch','on');
hold on;
xx = ones(1,length(RTlgn)).*(1+(rand(1,length(RTlgn))-0.5)/2.2);
scatter(xx, RTlgn,20,'o','filled','MarkerEdgeColor','k','MarkerFaceColor', rgb('DarkGray'));
xx = ones(1,length(RTdlgn)).*(2+(rand(1,length(RTdlgn))-0.5)/2.2); 
scatter(xx, RTdlgn,20,'o','filled','MarkerEdgeColor','k','MarkerFaceColor', rgb('DarkGray'));
xx = ones(1,length(RTlp)).*(3+(rand(1,length(RTlp))-0.5)/2.2);
scatter(xx, RTlp,20,'o','filled','MarkerEdgeColor','k','MarkerFaceColor', rgb('DarkGray'));
xx = ones(1,length(RTvb)).*(4+(rand(1,length(RTvb))-0.5)/2.2);
scatter(xx, RTvb,20,'o','filled','MarkerEdgeColor','k','MarkerFaceColor', rgb('DarkGray'));
xx = ones(1,length(RTvp)).*(5+(rand(1,length(RTvp))-0.5)/2.2);
scatter(xx, RTvp,20,'o','filled','MarkerEdgeColor','k','MarkerFaceColor', rgb('DarkGray'));
xx = ones(1,length(RTpom)).*(6+(rand(1,length(RTpom))-0.5)/2.2);
scatter(xx, RTpom,20,'o','filled','MarkerEdgeColor','k','MarkerFaceColor', rgb('DarkGray'));
xx = ones(1,length(RTvmgb)).*(7+(rand(1,length(RTvmgb))-0.5)/2.2); 
scatter(xx, RTvmgb,20,'o','filled','MarkerEdgeColor','k','MarkerFaceColor', rgb('DarkGray'));
xx = ones(1,length(RTdmgb)).*(8+(rand(1,length(RTdmgb))-0.5)/2.2); 
scatter(xx, RTdmgb,20,'o','filled','MarkerEdgeColor','k','MarkerFaceColor', rgb('DarkGray'));

xticklabels({'LGN','dLGN','LP','VB','VP','POm','vMGB','dMGB'});
ylabel('total rise time (in ms)');
title({'boxplot of time from hyperpolarization stop to first burst spike',...
    'for 8 different areas of brain'});
xlim([0.5,8.5]);

%% After last burst spike, time for 90% decrease in voltage
%% After last burst spike, time for 90% decrease in voltage histogram
aftBDecTimeL = inVitroBurstDataTable.newPerct90Dec(lowerOrderIndx & rowsNumSpksMoreThanOne);
aftBDecTimeH = inVitroBurstDataTable.newPerct90Dec(higherOrderIndx & rowsNumSpksMoreThanOne);
pkTimesL = inVitroBurstDataTable.peakTimes(lowerOrderIndx & rowsNumSpksMoreThanOne);
lastPkTimeL = zeros(size(pkTimesL));
for u = 1:size(lastPkTimeL,1)
    pks = pkTimesL{u,1};
    lastPkTimeL(u) = pks(end);
end
pkTimesH = inVitroBurstDataTable.peakTimes(higherOrderIndx & rowsNumSpksMoreThanOne);
lastPkTimeH = zeros(size(pkTimesH));
for u = 1:size(lastPkTimeH,1)
    pks = pkTimesH{u,1};
    lastPkTimeH(u) = pks(end);
end

aftBDecTimeFrBL =  aftBDecTimeL - lastPkTimeL;
indx = find(aftBDecTimeL == 1);
aftBDecTimeFrBL(indx) = [];
aftBDecTimeFrBH =  aftBDecTimeH - lastPkTimeH;

edgesPlot = 20:20:1800;
aftBDecTimeFrBH_SubSampled = datasample(aftBDecTimeFrBH,length(aftBDecTimeFrBL),'Replace',false); 
f1 = figure(1); bFig1(1:2) = axes(f1); % lower order
bFig1(1) = subplot(2,1,1); histogram(aftBDecTimeFrBL,edgesPlot); xlabel('Lower\_Order');
bFig1(2) = subplot(2,1,2); histogram(aftBDecTimeFrBH_SubSampled,edgesPlot); xlabel('Higher\_Order');
sgtitle(f1,{'Histogram of voltage decrease time after last burst spike (in ms)'});
linkaxes([bFig1(2),bFig1(1)],'xy');

%% Plot some sample examples for traces different voltage decrease time
iL = inVitroBurstDataTable.i(lowerOrderIndx & rowsNumSpksMoreThanOne); iL(indx) = [];
jL = inVitroBurstDataTable.j(lowerOrderIndx & rowsNumSpksMoreThanOne); jL(indx) = [];
iH = inVitroBurstDataTable.i(higherOrderIndx & rowsNumSpksMoreThanOne); 
jH = inVitroBurstDataTable.j(higherOrderIndx & rowsNumSpksMoreThanOne); 
hyperEndTimeL = inVitroBurstDataTable.hyperpolarStopTime(lowerOrderIndx & rowsNumSpksMoreThanOne);
hyperEndTimeH = inVitroBurstDataTable.hyperpolarStopTime(higherOrderIndx & rowsNumSpksMoreThanOne);
hyperEndTimeL(indx) = [];

figure; hold on; % Lower Order
level = [50,200,400,600,800,1000,1200,1600]';
subplt = [1,3,5,7,2,4,6,8]';
for w = 1:size(level,1)
    temp = [abs(aftBDecTimeFrBL - level(w)),[1:1:size(aftBDecTimeFrBL,1)]'];
    temp = sortrows(temp,1); temp  = temp(1:3,2);
    % [~,t] = min(abs(aftBDecTimeFrBL - 50));
    p = iL(temp); q = jL(temp); hyperEndTime = hyperEndTimeL(temp);
    for r = 1:size(p,1)
        [timeaxis, Vresponse] = openingInVitroTraceFiles(p(r),q(r));
        time = timeaxis(timeaxis >= hyperEndTime(r));
        vol = Vresponse(timeaxis >= hyperEndTime(r));
        subplot(4,2,subplt(w)); title(strcat(num2str(level(w)),'ms')); plot(time,vol); hold on; 
    end
end
suptitle('90PerctDecTime\_Lower\_Order');

figure; hold on; % Higher Order
level = [50,200,400,600]';
subplt = [1,3,2,4]';
for w = 1:size(level,1)
    temp = [abs(aftBDecTimeFrBH - level(w)),[1:1:size(aftBDecTimeFrBH,1)]'];
    temp = sortrows(temp,1); temp  = temp(1:3,2);
    % [~,t] = min(abs(aftBDecTimeFrBL - 50));
    p = iH(temp); q = jH(temp); hyperEndTime = hyperEndTimeH(temp);
    for r = 1:size(p,1)
        [timeaxis, Vresponse] = openingInVitroTraceFiles(p(r),q(r));
        time = timeaxis(timeaxis >= hyperEndTime(r));
        vol = Vresponse(timeaxis >= hyperEndTime(r));
        subplot(2,2,subplt(w)); title(strcat(num2str(level(w)),'ms')); plot(time,vol); hold on; 
    end
end
suptitle('90PerctDecTime\_Higher\_Order');


%% Plotting scatter of last peak's voltages for lower and higher order cells
lastPeakVolL = inVitroBurstDataTable.lastPeakVoltage(lowerOrderIndx & rowsNumSpksMoreThanOne);
lastPeakVolH = inVitroBurstDataTable.lastPeakVoltage(higherOrderIndx & rowsNumSpksMoreThanOne);
edgesPlot = -50:1:40;
lastPeakVolH_SubSampled = datasample(lastPeakVolH,length(lastPeakVolL),'Replace',false); 
f1 = figure(1); bFig1(1:2) = axes(f1); % lower order
bFig1(1) = subplot(2,1,1); histogram(lastPeakVolL,edgesPlot); xlabel('Lower\_Order'); 
bFig1(2) = subplot(2,1,2); histogram(lastPeakVolH_SubSampled,edgesPlot); xlabel('Higher\_Order');
sgtitle(f1,{'Histogram of last burst spike voltage(in ms)'});
linkaxes([bFig1(2),bFig1(1)],'xy');

%% After last burst spike, time for 90% decrease in voltage 
% only for traces with last spike voltage lower than -20mV
% Assumption here - that flater peaks are the calcium decrease ones
lastPeakVolL = inVitroBurstDataTable.lastPeakVoltage(lowerOrderIndx & rowsNumSpksMoreThanOne);
lastPeakVolH = inVitroBurstDataTable.lastPeakVoltage(higherOrderIndx & rowsNumSpksMoreThanOne);
aftBDecTimeL = cell2mat(inVitroBurstDataTable.perct90DecInVolAfterBurst(lowerOrderIndx & rowsNumSpksMoreThanOne));
aftBDecTimeH = cell2mat(inVitroBurstDataTable.perct90DecInVolAfterBurst(higherOrderIndx & rowsNumSpksMoreThanOne));
lastPeakTimeL = inVitroBurstDataTable.lastPeakTime(lowerOrderIndx & rowsNumSpksMoreThanOne);
lastPeakTimeH = inVitroBurstDataTable.lastPeakTime(higherOrderIndx & rowsNumSpksMoreThanOne);
aftBDecTimeFrBL = aftBDecTimeL - lastPeakTimeL;
aftBDecTimeFrBH = aftBDecTimeH - lastPeakTimeH;
aftBDecTimeFrBL_short = aftBDecTimeFrBL(lastPeakVolL <= -20);
aftBDecTimeFrBH_short = aftBDecTimeFrBH(lastPeakVolH <= -20);

edgesPlot = 20:20:1800;
aftBDecTimeFrBH_short_SubSampled = datasample(aftBDecTimeFrBH_short,...
    length(aftBDecTimeFrBL_short),'Replace',false); 
f1 = figure(1); bFig1(1:2) = axes(f1); % lower order
bFig1(1) = subplot(2,1,1); histogram(aftBDecTimeFrBL_short,edgesPlot); xlabel('Lower\_Order');
bFig1(2) = subplot(2,1,2); histogram(aftBDecTimeFrBH_short_SubSampled,edgesPlot); xlabel('Higher\_Order');
sgtitle(f1,{'Histogram of 90% voltage decrease time after last burst spike (in ms)','for last peak voltages <= -20mV '});
linkaxes([bFig1(2),bFig1(1)],'xy');

%% Resting Potential
%% Histogram of resting potential
VrestL = inVitroBurstDataTable.restingPotential(lowerOrderIndx & rowsNumSpksMoreThanOne);
VrestH = inVitroBurstDataTable.restingPotential(higherOrderIndx & rowsNumSpksMoreThanOne);
edgesPlot = -70:0.25:-50;
VrestH_SubSampled = datasample(VrestH,length(VrestL),'Replace',false); 
f1 = figure(1); bFig1(1:2) = axes(f1); % lower order
bFig1(1) = subplot(2,1,1); histogram(VrestL,edgesPlot); xlabel('Lower\_Order'); 
bFig1(2) = subplot(2,1,2); histogram(VrestH_SubSampled,edgesPlot); xlabel('Higher\_Order');
sgtitle(f1,{'Histogram of resting potential (in mV)'});
linkaxes([bFig1(2),bFig1(1)],'xy');

%% Savin new information in table
%% Save firstPeakTime and lastPeakTime separately in the data table
peakTimes = inVitroBurstDataTable.peakTimes;
firstPeakTime = zeros(size(peakTimes));
lastPeakTime = zeros(size(peakTimes));
for p = 1:size(peakTimes,1)
    temp = peakTimes{p,1};
    if isempty(temp)
    elseif size(temp) == 1
        firstPeakTime(p) = temp;
    else
        firstPeakTime(p) = temp(1);
        lastPeakTime(p) = temp(end);
    end
end  
inVitroBurstDataTable = addvars(inVitroBurstDataTable,firstPeakTime,'After','peakTimes');
inVitroBurstDataTable = addvars(inVitroBurstDataTable,lastPeakTime,'After','firstPeakTime');

%% Save last peak voltage separately in the data table
peakVoltages = inVitroBurstDataTable.peakVoltages;
lastPeakVoltage = zeros(size(peakTimes));
for p = 1:size(peakVoltages,1)
    temp = peakVoltages{p,1};
    if ~isempty(temp)
        lastPeakVoltage(p) = temp(end);
    end
end  
inVitroBurstDataTable = addvars(inVitroBurstDataTable,lastPeakVoltage,'After','peakVoltages');


end