function [inVitroBurstData_2_12] = extractBurstInformation_inVitroData()
% Extracting information about spikes from in-vitro data to characterize
% bursts properties in lower and higher order of visual, auditory and 
% somato-sensory regions of thalamus 
% The burst information has changed with better spikes found and saved in 
% variable called peaks 2_12

clear;
cd('C:\Users\Deepcutlab\Desktop\Nidhi\DataSet_Nidhi\In_Vitro_Burst_Data');
load('C:\Users\Deepcutlab\Desktop\Nidhi\Code_Nidhi\In-Vitro Burst Analysis\Results\peaks_SavedOn_2_12_2020.mat');

%% Identify .abf files in the current folder
xx = '*.abf*';
FnameUnits = dir(xx);
extractnames = extractfield(FnameUnits,'name');
NumberCells = size(extractnames,2);
inVitroBurstData_2_12 = struct();
counter = 0;

for i = 1:NumberCells
    %% Load data 
    fileName = extractnames{i};
    [d,si] = abfload(fileName);  % si = sampling interval in microsecs
    if size(d,1) == 0 % some data are newer so cannot be opened using abfload so use abf2load
        [d,si] = abf2load(fileName);
    end
    
    %% Some files have extra information in 2nd dimension, extract only voltage matrix
    if size(d,2) == 2
        VresponseAllInput = d(:,1,:); 
    else
        VresponseAllInput = d;
    end

    %% Get voltage and time axis
    VresponseAllInput = squeeze(VresponseAllInput);
    timeaxis = linspace(0,(size(VresponseAllInput,1)*(si./1000)),size(VresponseAllInput,1)); % timeaxis is in milliseconds
    timeaxis = timeaxis';
    
    %% Find the two points of transition as response to the current applied
    % Using only the first traces as it will be the one which is most hyperpolarized
    Vresponse = VresponseAllInput(:,1);
    [TF,~] = ischange(Vresponse,'mean','MaxNumChanges',2);
    ChangeIDx = find(TF == 1)';
    
    % Find start
    StartPulseWin = Vresponse(ChangeIDx(1)-1000:ChangeIDx(1),1);
    [TFstart,~] = ischange(StartPulseWin,'linear','MaxNumChanges',1);
    InputStartLoc = find(TFstart == 1)+ ChangeIDx(1)-1000;
    
    % Find End
    EndPulseWin = Vresponse(ChangeIDx(2)-1000:ChangeIDx(2),1);
    [TFend,~] = ischange(EndPulseWin,'linear','MaxNumChanges',1);
    InputEndLoc = find(TFend == 1)+ ChangeIDx(2)-1000;
    
    for j = 1:size(VresponseAllInput,2) % Going through each output recording for different input currents for each cell
        Vresponse = VresponseAllInput(:,j);
%         figure; plot(timeaxis,Vresponse,'b');

        %% Calculate resting potential
        % average voltage between 50 ms (to remove initial variation) and 250 ms (currents are not applied before this time)
        [~,Indx50ms]  = min(abs(timeaxis - 50));
        [~,Indx250ms] = min(abs(timeaxis - 250));
        Vrest = mean(Vresponse(Indx50ms:Indx250ms));
        
        %% Finding if the trace was depolarized or hyperpolarized by comparing voltage at the end of input pulse to Vrest
        V100msBeforeInputEnd = median(Vresponse((InputEndLoc-1500):(InputEndLoc-500)));
        if (Vrest - V100msBeforeInputEnd) > 5
            hyperpolarizingInput = 1; % hyperpolarizing input was provided
        elseif (V100msBeforeInputEnd - Vrest) > 5
            hyperpolarizingInput = 2; % depoolarizing input was provided
        else
            hyperpolarizingInput = 0; % No input provided or nothing happened (default)
        end
       
       %% Calculating different parameters for hyperpolarizing traces
        if hyperpolarizingInput == 1 
            counter = counter + 1;
            %% Slope at the end of input pulse
            slopeEndHyper = atand((Vresponse(InputEndLoc) - Vresponse(InputEndLoc-1000))...
                /(timeaxis(InputEndLoc)-timeaxis(InputEndLoc-1000)));

            %% Find peaks of the spikes in a burst
            indx = cell2mat(peaks_2_12(:,1)) == i & cell2mat(peaks_2_12(:,2)) == j;
            pks = peaks_2_12{indx,3};
            locs = peaks_2_12{indx,4};
            
            %% Burst spikes properties
            if size(locs,1) >= 1
                firstSpkTimeAfterInputEnd = locs(1) - timeaxis(InputEndLoc);
            else
                firstSpkTimeAfterInputEnd = [];
            end
            
            NumSpikes = size(pks,1);
            if NumSpikes >= 2
                allIBIs = zeros(NumSpikes-1,1);
                for q = 1:NumSpikes-1
                    allIBIs(q) = locs(q+1) - locs(q);
                end
            else
                allIBIs = [];
            end
            
            %% Calculate Ca2+ build-up time
            if size(locs,1) >= 1
                tt = timeaxis(InputEndLoc:find(timeaxis == locs(1)));
                v = Vresponse(InputEndLoc:find(timeaxis == locs(1)));
                s = []; % slope in degrees
                u = 1;
                tt2 = [];
                while u <= (size(v,1)-10)
                    s = [s; atand((v(u+10)- v(u))/ (timeaxis(u+10) - timeaxis(u)))];
                    tt2 = [tt2;tt(u)];
                    u = u + 10;
                end
                f = polyfit(tt2,s,4);
                y = polyval(f,tt2);
                [~,minIndx] = min(y);
                Ca2RiseSlopeChangeTime = tt2(minIndx);
            else 
                Ca2RiseSlopeChangeTime = [];
            end
                
            %% Calculate Ca2+ 90% decrease time
            if size(locs,1) >= 1
                % The time it takes for the voltage to decrease to 10% of the height of the last peak from resting potential
                Vref = Vrest + 0.1 * (pks(end) - Vrest);
                VafterLastSpike = Vresponse(find(timeaxis == locs(end)):end);
                TimesAfterLastSpike = timeaxis(find(timeaxis == locs(end)):end);
                [~,indx] = min(abs(Vref-VafterLastSpike)); % might lead to issue
                perct90DecAfterBurst = TimesAfterLastSpike(indx);
            else 
                perct90DecAfterBurst = [];
            end
                
            %% Calculate Ca2+ decrease path, the time when it reaches 10% above Vrest
%             indx = find(timeaxis == locs(end));
%             endTimes = timeaxis(indx:end);
%             endVol = Vresponse(indx:end);
%             
%             indx = find(timeaxis ==);
%             endTimes = timeaxis(indx+164:end);
%             endVol = Vresponse(indx+164:end);
%             figure; hold on; plot(endTimes,endVol,'o');
%             f1 = fit(endTimes,endVol,'exp2');
%             plot(f1);
%             
%             
%                 tt =  timeaxis(indx:end);
%                 v = Vresponse(indx:end);
%                 s = []; % slope in degrees
%                 u = 1;
%                 tt2 = [];
%                 steps = 100;
%                 while u <= (size(v,1)-steps)
%                     s = [s; atand((v(u+steps)- v(u))/ (timeaxis(u+steps) - timeaxis(u)))];
%                     tt2 = [tt2;tt(u)];
%                     u = u + 10;
%                 end
%                 f = polyfit(tt2,s,4);
%                 y = polyval(f,tt2);
%                 [~,minIndx] = min(y);
%                 Ca2RiseSlopeChangeTime = tt2(minIndx);
%          % Plot the slope trend 
%                figure; hold on; plot(tt2,s,'o');
%                 hold on;
%                 plot(f);
%                 hold off;
%             
%             % 
%             
%                 indx = find(timeaxis == locs(end));
%                 tt = timeaxis(indx:end);
%                 v = Vresponse(indx:end);
%                 s = []; % slope in degrees
%                 u = 1;
%                 tt2 = [];
%                 steps = 10;
%                 while u <= (size(v,1)-steps)
%                     s = [s; atand((v(u+steps)- v(u))/ (timeaxis(u+steps) - timeaxis(u)))];
%                     tt2 = [tt2;tt(u)];
%                     u = u + 10;
%                 end
%                 f = fit(tt2,s,'exp2');
%                 figure; hold on; 
%                 plot(tt2,s,'o');
%                 plot(f);
%                 hold off;
% 
%                 figure; plot(tt,v);
%                 
%                 figure; plot(timeaxis,Vresponse);
            
            %% Finding brain region for a particular cell and saving them to repsective field in structure 
            k = strfind(fileName,'.abf');
            fileName2 = char(fileName(1:k-1));
            k = strfind(fileName2,'_');
            brainRegion = fileName2(k+1:end);
            if strcmp(brainRegion,'LGN') || strcmp(brainRegion,'dLGN') 
                cellOrder = 'lowerOrder';
                thalamicPart = 'visual';
            elseif strcmp(brainRegion,'LP')
                cellOrder = 'higherOrder';
                thalamicPart = 'visual';
            elseif strcmp(brainRegion,'VB') || strcmp(brainRegion,'VP')
                cellOrder = 'lowerOrder';
                thalamicPart = 'somato';
            elseif strcmp(brainRegion,'POm')
                cellOrder = 'higherOrder';
                thalamicPart = 'somato';
            elseif strcmp(brainRegion,'vMGB')
                cellOrder = 'lowerOrder';
                thalamicPart = 'auditory';
            elseif strcmp(brainRegion,'dMGB') || strcmp(brainRegion,'dMGN')
                cellOrder = 'higherOrder';
                thalamicPart = 'auditory';
            end
            cellNumber = strcat(fileName2(1:k-1)); 
            
            %% Saving paramenters for a particular cell into a structure
            % Cell information
            inVitroBurstData_2_12(counter).i = i;
            inVitroBurstData_2_12(counter).j = j;
            inVitroBurstData_2_12(counter).cellNumber = cellNumber;
            inVitroBurstData_2_12(counter).brainRegion = brainRegion;
            inVitroBurstData_2_12(counter).cellOrder = cellOrder;
            inVitroBurstData_2_12(counter).thalamicPart = thalamicPart;
            
            % Hyperpolar information
            inVitroBurstData_2_12(counter).restingPotential = Vrest;
            inVitroBurstData_2_12(counter).hyperpolarVoltage = V100msBeforeInputEnd;
            inVitroBurstData_2_12(counter).hyperpolarStartTime = timeaxis(InputStartLoc);
            inVitroBurstData_2_12(counter).hyperpolarStopTime = timeaxis(InputEndLoc);           
            inVitroBurstData_2_12(counter).slopeAtEndOfHyperpolar = slopeEndHyper;
            
            % Bursting information
            inVitroBurstData_2_12(counter).peakTimes = locs;
            inVitroBurstData_2_12(counter).peakVoltages = pks;
            inVitroBurstData_2_12(counter).NumSpikes = NumSpikes;
            inVitroBurstData_2_12(counter).firstIBI = firstSpkTimeAfterInputEnd;
            inVitroBurstData_2_12(counter).allIBIs = allIBIs;
            
            % Calcium peak during bursting information
            inVitroBurstData_2_12(counter).Ca2RiseSlopeChangeTime = Ca2RiseSlopeChangeTime;
            inVitroBurstData_2_12(counter).perct90DecInVolAfterBurst = perct90DecAfterBurst;
            
        end
    end
end
end
