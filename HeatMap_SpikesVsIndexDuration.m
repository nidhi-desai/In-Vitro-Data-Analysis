function HeatMap_SpikesVsIndexDuration()
%%% X-axis each spindle's time, Y-axis Index Duration, 
%%% Z-axis Number of Bursts in each spindle

clear;   
load('C:\Users\Deepcutlab\Desktop\Nidhi\Code_Nidhi\PreDir_AllRats_AllSessions.mat');
cd ('C:\Users\Deepcutlab\Desktop\Nidhi\Code_Nidhi\Results\Bursts');
exportToPPTX('new','SpikesVsIndexDuration'); 
exportToPPTX('saveandclose','SpikesVsIndexDuration');

for i=1:size(PreDir,1)
    %% Loading data 
    currdir = PreDir(i,:);
    cd(currdir); 
    load('SpindleInfo.mat');
    load('SWStimespre.mat');
    thresholdDistribution = 6;

    %% For each thalamic electrode and units
    fields = fieldnames(SpindleInfo);
    indx = find(contains(fields,'SpikePerCycleThalamus_'));
    
    for j = 1:size(indx,1)
        %% Getting data to plot
        numCycles = SpindleInfo.NumOfCyclesEachSpindle;
        indexDuration = SpindleInfo.indexDuration;
        spindleTimes = SpindleInfo.SpindleStartEndTimes;
        cyclesEachSpindle = SpindleInfo.CycleNumberOccuringInSameSpindle;
        spikesPerCycle = SpindleInfo.(char(fields(indx(j))));
        shortlistIndexDuration = [];
        shortlistCyclesEachSpindle = [];
        numSpikesPerSpindle = [];
        spindleStartTime  = [];
        for p = 1:size(numCycles,1)
            if numCycles(p) >= thresholdDistribution
                shortlistIndexDuration = [shortlistIndexDuration; indexDuration(p)];
                spindleStartTime = [spindleStartTime; spindleTimes(p,1)];
                temp = cyclesEachSpindle{p};
                shortlistCyclesEachSpindle = [shortlistCyclesEachSpindle; cyclesEachSpindle(p)];
                tempSum = 0;
                for r = 1:size(temp,2)
                    tempSum = tempSum + spikesPerCycle(temp(r));
                end
                numSpikesPerSpindle = [numSpikesPerSpindle; tempSum];
            end
        end
           
        yblocks = [ceil(max(shortlistIndexDuration)):-1:floor(min(shortlistIndexDuration))];
        yblocks = yblocks';
        Y = zeros(size(shortlistIndexDuration,1),1);
        for k = 1:size(shortlistIndexDuration,1)
            Y(k) = floor(shortlistIndexDuration(k));
        end

        zeroline = find(yblocks==0); 
        mat = zeros(size(yblocks,1), size(shortlistIndexDuration,1));
        for r = 1:size(mat,2)
            if Y(r) < 0
                z = zeroline;
                counter = 0;
                for o = 1:abs(Y(r))
                    mat(z + counter, r) = numSpikesPerSpindle(r);
                    counter = counter + 1;
                end
            else
                z = zeroline-1;
                counter = 0;
                for o = 1:abs(Y(r))
                    mat(z - counter, r) = numSpikesPerSpindle(r);
                    counter = counter + 1;
                end
            end
        end
                
        %% Plotting data
        scrsz = get(0,'screensize');
        hl = figure('units','normalized','outerposition',[0 0 1 1]);
        hl = imagesc(mat);
        hold on;
        colormap(hot);
        clTicks = 0:1:max(max(mat));
        cl = colorbar('Ticks',clTicks);
        yline(zeroline-0.5,'Color','white');

        %caxis([0 0.15]);
        %clb = colorbar;
        %clb.Ticks = linspace(0, 0.15, 7);
        
        % Y ticks and Y tick labels
        yticks(0.5:1:size(mat,1)+0.5);
        yticklabels(yblocks);
        ylabel('Index Duration');
        ax1 = gca;

        % X ticks and X tick labels
        xticksGap = round(size(mat,2)/15);
        r = 1:xticksGap:size(mat,2);
        r = [r size(spindleStartTime,2)];
        r = unique(r);
        set(ax1, 'XTick', r);
        x_labels = [];
        for u = 1:size(r,2)
           x_labels = [x_labels; spindleStartTime(r(u))];
        end
        x_labels = [x_labels; spindleStartTime(end)];
        x_labels = round(x_labels);
        set(ax1, 'XTickLabel', x_labels);
        %xticklabels(x_labels);
        xlabel('Time stamp (in seconds)'); % Confirm the units of the time
        ylabel(cl, 'Number of Spikes');
        
       % Give a super title with Rat name and session number 
        ID1 = strfind(currdir, '\');
        ID2 = ID1(length(ID1)-1);
        name3 = currdir(ID2+1:end);
        name3 = strrep(name3, '\', strcat('\', '_'));
        
        elect = char(fields(indx(j)));
        ID1= strfind(elect, '_');
        ID2= ID1(1);
        name4 = elect(ID2:end);
        name4 = strrep(name4, '\', strcat('\','_'));
        name4 = strrep(name4, '_', strcat('\','_'));
        title({'HeatMap Average difference in duration of cycles vs. num of spikes for each spindle', strcat(char(name3), char(name4))});

        % Mark the different sleep periods
        spindleNumBeforeEndSleepPeriod = [];
        for p = 2:size(SWStimes,1)
               idx = find(spindleStartTime > SWStimes(p,1),1);
               spindleNumBeforeEndSleepPeriod = [spindleNumBeforeEndSleepPeriod idx];
        end
        spindleNumBeforeEndSleepPeriod = [spindleNumBeforeEndSleepPeriod size(mat,2)];
        for t = 1:size(spindleNumBeforeEndSleepPeriod,2)
           xl = xline(spindleNumBeforeEndSleepPeriod(t)+0.5,'-',{strcat('SWS period  ', num2str(t)), strcat(num2str(round(SWStimes(t,2))), ' secs')});
           xl.Color = 'cyan';
           xl.LabelHorizontalAlignment = 'left';
           xl.LabelVerticalAlignment = 'top';
           xl.LabelOrientation = 'aligned';
        end 

        
        %% Saving plots to PPTX
        cd ('C:\Users\Deepcutlab\Desktop\Nidhi\Code_Nidhi\Results\Bursts');
        exportToPPTX('open','SpikesVsIndexDuration'); 
        exportToPPTX('addslide');
        exportToPPTX('addpicture',figure(1));
        exportToPPTX('saveandclose','SpikesVsIndexDuration');
        close(gcf);
        
    end
end

end
