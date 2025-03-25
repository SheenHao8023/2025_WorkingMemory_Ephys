%% 在单个样本文件夹下需要三个数据文件
% allData.mat、behavior.mat，以及数据日期子文件夹下setting.xml
clc
close all

% 建立所有数据的路径
all_animals = {
    'C:\Users\XinHao\Desktop\WM_DualCircuit\data\CIBRZC85\';...
};

n_animals = length(all_animals);
for i_mice = 1:n_animals
    % select recording data
    i_str = strfind(all_animals{i_mice},'\');
    data_file_tmp = dir([all_animals{i_mice},'*_behavior.mat']);
    data_file = {};
    for i_file = 1:length(data_file_tmp)
        % --------------- load solo data -------------------
        load([all_animals{i_mice},data_file_tmp(i_file).name]);
        if ~isempty(strfind(obj.sessionNameTag,'recording'))|~isempty(strfind(obj.sessionNameTag,'Recording')) 
            filename_tmp = data_file_tmp(i_file).name;
            i_str = strfind(filename_tmp,'_behavior');
            data_file{end+1,1} = [filename_tmp(1:i_str-1),'_allData.mat'];
        end
    end
    % go through recording data
    for i_file = 1   %1:length(data_file)-1
        if exist([all_animals{i_mice},data_file{i_file}],"file")
            % --------------- process solo data -------------------
            clear obj noise waveform_tmp unit
            load([all_animals{i_mice},data_file{i_file}]); % this m-file contains all the depth + laser position information
            % --------------- load recording depth + laser position data -----------------
            meta_neuronal_opto_behavior_McircuitModule
            % animal name
            i_str = strfind(all_animals{i_mice},'\');
            animal_ID_tmp = all_animals{i_mice}(i_str(end-1)+1:i_str(end)-1);
            % session date
            i_str = strfind(data_file{i_file},'_');
            session_date_tmp = data_file{i_file}(i_str(1)+1:i_str(end)-1);
            % find session
            ind_session = find(strcmp(animal_ID_tmp,mice_id_allSessions) & strcmp(session_date_tmp,date_allSessions));
            % load variables
            recording_depth_iSession = recording_depth{ind_session};
            probe_type_iSession = probe_type(ind_session,:);
            Trim=trimend{ind_session};         
            %disp([all_animals{i_mice},data_file{i_file}]);

            Path       = [all_animals{i_mice},session_date_tmp,'\']; % wherever you want to search
            searchPath = [Path ,'\**\*.xml']; % Search in folder and subfolders for  *.xml
            Files      = dir(searchPath); % Find all .xml files
            XMLfile=[Files.folder,'\',Files.name];%only fo NP
            [xcoords, ycoords]=createChannelMapFile_384_NP1(XMLfile,1,Path);

            obj.Alarm_Nums();
            obj.Pole_Time();
            obj.Cue_Time();
            
            % --------------- get data -----------------
            R_hit_tmp = ((char(obj.sides)=='r') & obj.trials.hitHistory);
            R_miss_tmp = ((char(obj.sides)=='r') & obj.trials.missHistory);
            R_ignore_tmp = ((char(obj.sides)=='r') & obj.trials.noResponseHistory);
            L_hit_tmp = ((char(obj.sides)=='l') & obj.trials.hitHistory);
            L_miss_tmp = ((char(obj.sides)=='l') & obj.trials.missHistory);
            L_ignore_tmp = ((char(obj.sides)=='l') & obj.trials.noResponseHistory);

            Pole_Time_tmp = obj.poleTime;
            Cue_Time_tmp = obj.cueTime;

            Pole_Time_tmp(:,2) = Pole_Time_tmp(:,2)+.1;  % add 0.1 s to pole offset so pole off is the end of sample period
            
            LickEarly_tmp = zeros(length(obj.eventsHistory),1);
            LickEarly_tmp(obj.trials.EarlyLicksHistory,1) = 1;
            
            delay_time=median(Cue_Time_tmp(find(LickEarly_tmp==0),1)-Pole_Time_tmp(find(LickEarly_tmp==0),2));
            sample_time=median(Pole_Time_tmp(find(LickEarly_tmp==0),2)-Pole_Time_tmp(find(LickEarly_tmp==0),1));
                      
            Stimtype=zeros(length(obj.stimProb),1);

            for i_solo_trial = 1:min([size(obj.wavesurfer.photo_input_trace1,1) length(obj.eventsHistory)])
                if ~isempty(obj.sidesTypes)
                    % get sample/delay trial
                    if ~isempty(strfind(obj.sidesTypes{i_solo_trial},'sample_period'))
                        Sample_Delay_tmp(i_solo_trial,1) = 1;
                    elseif ~isempty(strfind(obj.sidesTypes{i_solo_trial},'delay_period'))
                        Sample_Delay_tmp(i_solo_trial,1) = 2;
                    elseif ~isempty(strfind(obj.sidesTypes{i_solo_trial},'response_period'))
                        Sample_Delay_tmp(i_solo_trial,1) = 3;
                    elseif ~isempty(strfind(obj.sidesTypes{i_solo_trial},'_n'))
                        Sample_Delay_tmp(i_solo_trial,1) = 0;
                    else
                    error('unrecognized stim type selection')
                    end
                else
                Sample_Delay_tmp(i_solo_trial,1) = 0;
                end
            end

            % Stimtype(find(obj.stimProb==0))=100;
            Stimtype(obj.trim)=100;
            if length(Trim)==1
                Stimtype(Trim:end)=100;%trim the end trials when mouse stop licking
            else
                Stimtype(Trim)=100;
            end

            % choose only the stim/performing/recording periods
            if islogical(obj.trim) && (length(obj.trim) == length(obj.eventsHistory))
                i_good_trials = (~obj.trim);
            elseif (length(obj.trim) ~= length(obj.eventsHistory))
                i_good_trials = ones(length(obj.eventsHistory),1);
                i_good_trials(obj.trim) = 0;
            else
                error('');
            end

            if length(Trim)==1
                i_good_trials(Trim:end) = 0;
            else
                i_good_trials(Trim) = 0;
            end

            i_good_trials(min([size(obj.wavesurfer.photo_input_trace1,1) size(obj.eventsHistory,1)]):end)=0;
            
            % -------------- process single unit data --------------         
            for i_probe=1:length(obj.units)
                if length(obj.units{i_probe})>1
                    i_unit=0;
                    Unit={};
                    PSTH={};
                    Trials={};
                    sig_selective = [];
                    Celltype=[];
                    FR_pref=[];
                    N_trials_all=[];
                    unitCh=[];
                    Depth_all=[];
                    Spikequality=[];    
                    for i_unit_tmp = 1:length(obj.units{i_probe})
                        unit = obj.units{i_probe}{i_unit_tmp};
                
                        % if (unit.manual_quality_score=='1')
                         i_unit=i_unit+1;
                  
                        spikeamp_tmp=mean(unit.amplitudes);
                        spikeisi_tmp=diff(unit.spike_times);
                        spikeisi_tmp=length(find(abs(spikeisi_tmp)<0.0015))/length(spikeisi_tmp);
                        spikecount=length(unit.spike_times);
                        spiketrialnum=length(unique(unit.trials));

                        % if unit.cell_type == 1;  index stable trials
                        unit.stable_trials=min(unit.stable_trials):max(unit.stable_trials);
                        unit.stable_trials(unit.stable_trials>length(obj.eventsHistory))=[];
                        unit.stable_trials(unit.stable_trials<=0)=[];
                        unit_stable_trials_tmp = zeros(length(obj.eventsHistory),1);
                        unit_stable_trials_tmp(unit.stable_trials) = 1;
                        
                        %yes-posterial-right no-anterial-left
                        i_yes_trial_correct = find((R_hit_tmp==1) & Stimtype==0 & unit_stable_trials_tmp==1 & LickEarly_tmp==0 & i_good_trials==1);
                        i_no_trial_correct = find((L_hit_tmp==1) & Stimtype==0 & unit_stable_trials_tmp==1 & LickEarly_tmp==0 & i_good_trials==1);
                        i_yes_trial_error = find((R_miss_tmp==1) & Stimtype==0 & unit_stable_trials_tmp==1 & LickEarly_tmp==0 & i_good_trials==1);
                        i_no_trial_error = find((L_miss_tmp==1) & Stimtype==0 & unit_stable_trials_tmp==1 & LickEarly_tmp==0 & i_good_trials==1);
                                                                
                        % log the trial numbers
                        N_trials_all(i_unit,:) = [
                            length(i_yes_trial_correct),...
                            length(i_no_trial_correct),...
                            length(i_yes_trial_error),...
                            length(i_no_trial_error),...
                        ];
                    
                        psth1 = nan;
                        spike_times1 = {};
                        PSTH1=nan;
                        if length(i_yes_trial_correct)>=3 & length(i_no_trial_correct)>=3 
                            for trialtype_tmp=1:4
                                switch trialtype_tmp
                                    case 1
                                        trial_selected=i_yes_trial_correct';
                                    case 2
                                        trial_selected=i_no_trial_correct';
                                    case 3
                                        trial_selected=i_yes_trial_error';
                                    case 4
                                        trial_selected=i_no_trial_error';
                                    case 5
                                end
                            
                                spike_times_psth = {};
                                spk_count1 = [];
                                spk_count2 = [];
                                spk_count3 = [];
                                spk_count4 = [];
                                n_trial = 0;
                                if length(trial_selected)>=1 
                                    for i_trial = trial_selected
                                        n_trial = n_trial+1;
                                        spike_times_psth{n_trial,1} = unit.spike_times(unit.trials==i_trial) - Cue_Time_tmp(i_trial,1);
                                        spk_count1(end+1,1) = sum(unit.trials==i_trial & unit.spike_times>Pole_Time_tmp(i_trial,1) & unit.spike_times<Pole_Time_tmp(i_trial,2))/(Pole_Time_tmp(i_trial,2)-Pole_Time_tmp(i_trial,1));
                                        spk_count2(end+1,1) = sum(unit.trials==i_trial & unit.spike_times>Pole_Time_tmp(i_trial,2) & unit.spike_times<Cue_Time_tmp(i_trial,1))/(Cue_Time_tmp(i_trial,1)-Pole_Time_tmp(i_trial,2));
                                        spk_count3(end+1,1) = sum(unit.trials==i_trial & unit.spike_times>Cue_Time_tmp(i_trial,1) & unit.spike_times<Cue_Time_tmp(i_trial,1)+1.3)/1.3;
                                        spk_count4(end+1,1) = sum(unit.trials==i_trial & unit.spike_times>Pole_Time_tmp(i_trial,1) & unit.spike_times<Cue_Time_tmp(i_trial,1)+1.3)/(Cue_Time_tmp(i_trial,1)+1.3-Pole_Time_tmp(i_trial,1));
                                    end
                                    [psth1, t, PSTH1] = func_getPSTH_10msbin(spike_times_psth,-3.5,2);
                                    spike_times1 = spike_times_psth;
                                else
                                    psth1 = nan;
                                    spike_times1 = {};                              
                                    PSTH1=nan;
                                end
                            
                                Unit{i_unit,(trialtype_tmp-1)*2+1}=psth1;
                                Unit{i_unit,(trialtype_tmp-1)*2+2}=spike_times1;                            
                                PSTH{i_unit,trialtype_tmp}=PSTH1;                            
                                Trials{i_unit,trialtype_tmp}=trial_selected;
                            
                                if trialtype_tmp==1
                                    spk_count_yes1=spk_count1;
                                    spk_count_yes2=spk_count2;
                                    spk_count_yes3=spk_count3;
                                    spk_count_yes4=spk_count4;
                                end
                                if trialtype_tmp==2
                                    spk_count_no1=spk_count1;
                                    spk_count_no2=spk_count2;
                                    spk_count_no3=spk_count3;
                                    spk_count_no4=spk_count4;
                                end
                           end

                            if ~isempty(spk_count_yes1) & ~isempty(spk_count_no1)
                                 sample_selective = ttest2(spk_count_yes1,spk_count_no1);
                            else
                                 sample_selective = nan;
                            end
                            
                            if ~isempty(spk_count_yes2) & ~isempty(spk_count_no2)
                                 delay_selective = ttest2(spk_count_yes2,spk_count_no2);
                            else
                                 delay_selective = nan;
                            end
                            
                            if ~isempty(spk_count_yes3) & ~isempty(spk_count_no3)
                                 response_selective = ttest2(spk_count_yes3,spk_count_no3);
                            else
                                 response_selective = nan;
                            end
                            
                            if ~isempty(spk_count_yes4) & ~isempty(spk_count_no4)
                                 sdr_selective = ttest2(spk_count_yes4,spk_count_no4);
                            else
                                 sdr_selective = nan;
                            end
                            
                            sig_selective(i_unit,1:4) = [sample_selective delay_selective response_selective sdr_selective];
                            Celltype(i_unit)=unit.cell_type;
                            FR_pref(i_unit,:) = [mean(spk_count_yes1) mean(spk_count_yes2) mean(spk_count_yes3) mean(spk_count_yes4)]...
                                -[mean(spk_count_no1) mean(spk_count_no2) mean(spk_count_no3) mean(spk_count_no4)];
                                                     
                            unitCh(i_unit)=median(unit.channel);
                            if i_probe==2
                                 if median(unit.channel)<=32
                                    Depth_all(i_unit,1) = recording_depth_iSession(i_probe)-(probe_type_iSession(2)+1-median(unit.channel))*25;%DBC2*32
                                 else
                                    Depth_all(i_unit,1) = recording_depth_iSession(i_probe)-(probe_type_iSession(2)+1-(median(unit.channel)-32))*25;%DBC2*32
                                 end
                            else
                                 Depth_all(i_unit)=recording_depth_iSession(i_probe)-ycoords(unitCh(i_unit));
                            end
                            Spikequality(i_unit,:)=[spikeamp_tmp spikeisi_tmp spikecount spiketrialnum];
                          
                            fh=figure('position',[320         633        1312         340]);
                            for fk=1:3                                                     
                                 switch fk
                                     case 1
                                         psthA=Unit{i_unit,1};
                                         psthB=Unit{i_unit,3};
                                         spike_timesA=Unit{i_unit,2};
                                         spike_timesB=Unit{i_unit,4};
                                         titletext='correct';
                                     case 2
                                         psthA=Unit{i_unit,5};
                                         psthB=Unit{i_unit,7};
                                         spike_timesA=Unit{i_unit,6};
                                         spike_timesB=Unit{i_unit,8};
                                         titletext='error';
                                     case 3
                                         if numel(PSTH{i_unit,1})>1&numel(PSTH{i_unit,3})>1
                                             psthA=mean([PSTH{i_unit,1};PSTH{i_unit,3}]);
                                         elseif numel(PSTH{i_unit,1})>1&numel(PSTH{i_unit,3})<=1
                                             psthA=Unit{i_unit,1};
                                         elseif numel(PSTH{i_unit,1})<=1&numel(PSTH{i_unit,3})>1
                                             psthA=Unit{i_unit,5};
                                         else
                                             psthA=Unit{i_unit,1};
                                         end % case-3 if-1 end

                                        if numel(PSTH{i_unit,2})>1&numel(PSTH{i_unit,4})>1
                                            psthB=mean([PSTH{i_unit,2};PSTH{i_unit,4}]);
                                        elseif numel(PSTH{i_unit,2})>1&numel(PSTH{i_unit,4})<=1
                                            psthB=Unit{i_unit,3};
                                        elseif numel(PSTH{i_unit,2})<=1&numel(PSTH{i_unit,4})>1
                                            psthB=Unit{i_unit,7};
                                        else
                                            psthB=Unit{i_unit,3};
                                        end % case-3 if-2 end

                                        spike_timesA=Unit{i_unit,2};
                                        if size(Unit{i_unit,6},1)>=1
                                            for jt=1:size(Unit{i_unit,6},1)
                                                spike_timesA{size(Unit{i_unit,2},1)+jt,1}=Unit{i_unit,6}{jt,1};
                                            end
                                        end
                                        spike_timesB=Unit{i_unit,4};
                                        if size(Unit{i_unit,8},1)>=1
                                            for jt=1:size(Unit{i_unit,8},1)
                                                spike_timesB{size(Unit{i_unit,4},1)+jt,1}=Unit{i_unit,8}{jt,1};
                                            end
                                        end
                                        titletext='Both';
                                 end  % switch end

                                 subplot(3,1,fk); hold on
                                 if ~isempty(spike_timesA) & ~isempty(spike_timesB) 
                                    % if length(t) == length(psthA) && length(t) == length(psthB) % 判断 t 和 psthA、psthB 是否长度一致
                                    %     plot(t, psthA, 'b');
                                    %     plot(t, psthB, 'r');
                                    % else
                                    %     disp('Error: t, psthA, and psthB must have the same length.');
                                    % end
                                    plot(t, psthA, 'b');
                                    plot(t, psthB, 'r');
                                    line([0 0],[0 2.2]*max([psthA psthB]),'color','k')
                                    line([-delay_time -delay_time],[0 2.2]*max([psthA psthB]),'color','k')
                                    line([-delay_time-sample_time -delay_time-sample_time],[0 2.2]*max([psthA psthB]),'color','k')
                                    xlim([-3.5 2])
                                    y_max = max([psthA psthB])*1.2;
                                    y_scale = max([psthA psthB]);
                                    n_trials = size(spike_timesA,1)+size(spike_timesB,1);
                                
                                    for i=1:length(spike_timesB)
                                        if ~isempty(spike_timesB{i})
                                            line([spike_timesB{i} spike_timesB{i}]', [y_max+y_scale/n_trials*(i-1) y_max+y_scale/n_trials*i],'color','r')
                                        end
                                    end
                                    y_max = y_max+y_scale/n_trials*i;
                                    for i=1:length(spike_timesA)
                                        if ~isempty(spike_timesA{i})
                                            line([spike_timesA{i} spike_timesA{i}]', [y_max+y_scale/n_trials*(i-1) y_max+y_scale/n_trials*i],'color','b')
                                        end
                                    end
                                
                                    if max([psthA psthB])~=0
                                        ylim([0 2.3]*max([psthA psthB]));
                                    end
                                end
                                title(titletext);

                                if fk==1
                                    title([titletext,':','ch',num2str(median(unit.channel)),'depth(um)',num2str(Depth_all(i_unit))])
                                end
                            end % for fk end
                           
                           % keyboard
                           figpath=[all_animals{i_mice},session_date_tmp,'\SingleUnits_10msbin_Fig',num2str(i_probe)];
                           matpath=[all_animals{i_mice},session_date_tmp,'\SingleUnits_10msbin_Mat',num2str(i_probe)];
                           tfdir=isfolder(figpath);
                           if tfdir==0
                                mkdir(figpath);
                           end
                           tfdir=isfolder(matpath);
                           if tfdir==0
                                mkdir(matpath);
                           end
                           %  saveas(fh,[figpath,'\SingleUnit',num2str(i_unit),'.fig'],'fig');
                           saveas(fh,[figpath,'\SingleUnit',num2str(i_unit),'.png'],'png');
                           close;

                        else 
                           sig_selective(i_unit,1:4) = [NaN NaN NaN NaN];
                           Celltype(i_unit)=unit.cell_type;
                           FR_pref(i_unit,:) = [NaN NaN NaN NaN];              
                           unitCh(i_unit)=median(unit.channel);
                           if i_probe==2
                                if median(unit.channel)<=32
                                    Depth_all(i_unit,1) = recording_depth_iSession(i_probe)-(probe_type_iSession(2)+1-median(unit.channel))*25;%DBC2*32
                                else
                                    Depth_all(i_unit,1) = recording_depth_iSession(i_probe)-(probe_type_iSession(2)+1-(median(unit.channel)-32))*25;%DBC2*32
                                end
                           else
                                Depth_all(i_unit)=recording_depth_iSession(i_probe)/1000-ycoords(unitCh(i_unit))/1000;
                           end
                           Spikequality(i_unit,:)=[spikeamp_tmp spikeisi_tmp spikecount spiketrialnum];
                        end % length(i_yes_trial_correct)>=3 & length(i_no_trial_correct)>=3 end
                        
                    end   % i_unit_tmp end    
                    save([matpath,'\SingleUnits_10msbin_raster_psth.mat'],'Spikequality','Unit','PSTH','Trials','sig_selective','Celltype','FR_pref','N_trials_all','t','unitCh','Depth_all','Pole_Time_tmp','Cue_Time_tmp','-v7.3');
                    clear Unit PSTH Trials;
                end % if length(obj.units{i_probe})>1 end
            end % i_probe end
        
        end
    end % i_file end

end % i_mice end

%% unit stable trials plot
path='C:\Users\XinHao\Desktop\WM_DualCircuit\data\CIBRZC85\2024_01_15\SingleUnits_10msbin_Mat2\';
load([path,'SingleUnits_10msbin_raster_psth.mat']);
figure;
for i_unit=1:size(Trials,1)
    i_unit;
    trials=[];
    for i_trial=1:size(Trials,2)
        trials=[trials Trials{i_unit,i_trial}];
    end
    d=rand(1)*20;
    plot([min(trials) max(trials)],[-Depth_all(i_unit)+d -Depth_all(i_unit)+d],'b-'); hold on;
end
xlabel('trial');
ylabel('depth(um)');
saveas(gcf,[path(1:end-5),'Fig2\unit_stable_trials_plot.fig']);

%% SC SNr sig selective map
c=zeros(size(FR_pref));
c(find(FR_pref>=0))='b';
c(find(FR_pref<0))='r';
name={'sample';'delay';'response';'sdr'};
figure;
for j=1:size(FR_pref,2)
    subplot(2,2,j);
    scatter(rand(1,size(FR_pref,1))*10,-Depth_all,10*abs(FR_pref(:,j)).*(sig_selective(:,j)+0.5),c(:,j),'filled');
    xlabel('x-jitter');
    ylabel('y-depth(um)');
    title(name{j});
end
saveas(gcf,'C:\Users\XinHao\Desktop\WM_DualCircuit\data\CIBRZC85\2024_01_15\SingleUnits_10msbin_Fig1\sig_selective_chandepth_map.fig');

%% ALM SC SNr CD corr
CD_proj_all={};
trial_all={};
for jgroup=1:4
    if jgroup==1
       obj0=load('C:\Users\XinHao\Desktop\WM_DualCircuit\data\CIBRZC85\2024_01_15\SingleUnits_10msbin_Mat1\SingleUnits_10msbin_raster_psth.mat');
    elseif jgroup==4
       obj0=load('C:\Users\XinHao\Desktop\WM_DualCircuit\data\CIBRZC85\2024_01_15\SingleUnits_10msbin_Mat2\SingleUnits_10msbin_raster_psth.mat');
    end
    t=obj0.t;
    sampletime=1.3;
    delaytime=1.7;

    if jgroup==1
        depth_range=[1500 2500];%SC-like
    elseif jgroup==2
        depth_range=[2500 3500];%middle
    elseif jgroup==3
        depth_range=[3500 4500];%SNr-like
    else
        depth_range=[0 1000];%ALM
    end

    i_yes_trial_correct=[];
    i_no_trial_correct=[];
    i_yes_trial_error=[];
    i_no_trial_error=[];
    population_stable_trials = [];

    for i_unit = 1:size(obj0.Unit,1)     
        i_yes_trial_correct=[i_yes_trial_correct obj0.Trials{i_unit,1}];
        i_no_trial_correct=[i_no_trial_correct obj0.Trials{i_unit,2}];
        i_yes_trial_error=[i_yes_trial_error obj0.Trials{i_unit,3}];
        i_no_trial_error=[i_no_trial_error obj0.Trials{i_unit,4}];
        
        trial_min=min([obj0.Trials{i_unit,1} obj0.Trials{i_unit,2} obj0.Trials{i_unit,3} obj0.Trials{i_unit,4}]);
        trial_max=max([obj0.Trials{i_unit,1} obj0.Trials{i_unit,2} obj0.Trials{i_unit,3} obj0.Trials{i_unit,4}]);
        % ---------- trial selection -----------------  index stable trials
        unit.stable_trials=trial_min:trial_max;
        unit_stable_trials_tmp = zeros(1000,1);
        unit_stable_trials_tmp(unit.stable_trials) = 1;
        population_stable_trials(end+1,:) = unit_stable_trials_tmp;
    end
    trialmax=find(sum(population_stable_trials)>0);
    trialmax=max(trialmax);
    unit_list=[];

    for i_unit=1:size(obj0.Unit,1)   
        if obj0.Depth_all(i_unit)>depth_range(1)&&obj0.Depth_all(i_unit)<depth_range(2)&&...
            min(find(population_stable_trials(i_unit,:)>0))<trialmax/4&&...
            max(find(population_stable_trials(i_unit,:)>0))>trialmax*3/4
            unit_list(i_unit)=1;
        else
            unit_list(i_unit)=0;
        end
    end
    population_stable_trials(find(unit_list==0),:)=[];
    population_stable_trials_tmp = find(sum(population_stable_trials)==size(population_stable_trials,1));

    i_yes_trial_correct=intersect(population_stable_trials_tmp,i_yes_trial_correct);
    i_no_trial_correct=intersect(population_stable_trials_tmp,i_no_trial_correct);
    i_yes_trial_error=intersect(population_stable_trials_tmp,i_yes_trial_error);
    i_no_trial_error=intersect(population_stable_trials_tmp,i_no_trial_error);

    psth_R_correct=[];
    psth_L_correct=[];
    psth_R_error=[];
    psth_L_error=[];
    cdprj_R_correct=[];
    cdprj_L_correct=[];
    cdprj_R_error=[];
    cdprj_L_error=[];

    n=0;
    for j=find(unit_list)
        n=n+1;
        psth_R_correct(:,:,n)=obj0.PSTH{j,1}(find(ismember(obj0.Trials{j,1},i_yes_trial_correct)),:)';
        psth_L_correct(:,:,n)=obj0.PSTH{j,2}(find(ismember(obj0.Trials{j,2},i_no_trial_correct)),:)';
        psth_R_error(:,:,n)=obj0.PSTH{j,3}(find(ismember(obj0.Trials{j,3},i_yes_trial_error)),:)';
        psth_L_error(:,:,n)=obj0.PSTH{j,4}(find(ismember(obj0.Trials{j,4},i_no_trial_error)),:)';
    end

    CD_trialtype=squeeze(mean(psth_R_correct,2))-squeeze(mean(psth_L_correct,2));
    CD_trialtype=normc(CD_trialtype');
    CD_delay=mean(CD_trialtype(:,find(t<-0.2&t>-delaytime)),2);

    for trial=1:size(psth_R_correct,2)
        cdprj_R_correct(:,trial)=squeeze(psth_R_correct(:,trial,:))*normc(CD_delay);
    end
    for trial=1:size(psth_L_correct,2)
        cdprj_L_correct(:,trial)=squeeze(psth_L_correct(:,trial,:))*normc(CD_delay);
    end
    for trial=1:size(psth_R_error,2)
        cdprj_R_error(:,trial)=squeeze(psth_R_error(:,trial,:))*normc(CD_delay);
    end
    for trial=1:size(psth_L_error,2)
        cdprj_L_error(:,trial)=squeeze(psth_L_error(:,trial,:))*normc(CD_delay);
    end

    CD_proj_all{jgroup,1}=cdprj_R_correct;
    CD_proj_all{jgroup,2}=cdprj_L_correct;
    CD_proj_all{jgroup,3}=cdprj_R_error;
    CD_proj_all{jgroup,4}=cdprj_L_error;
    trial_all{jgroup,1}=i_yes_trial_correct;
    trial_all{jgroup,2}=i_no_trial_correct;
    trial_all{jgroup,3}=i_yes_trial_error;
    trial_all{jgroup,4}=i_no_trial_error;
    psth_all{jgroup,1}=psth_R_correct;
    psth_all{jgroup,2}=psth_L_correct;
    psth_all{jgroup,3}=psth_R_error;
    psth_all{jgroup,4}=psth_L_error;
end

%% Functional_connectivity_ XH_250323
% 首先需要 run 'ALM SC SNr CD corr' session

region = {'SC';'MRN';'SNr';'ALM'};
z = 0;
Pearson_correlation_all = cell(6,4); 
Cross_correlation_all = cell(6,4); 
Granger_causality_A_all = cell(6,4); 
Granger_causality_B_all = cell(6,4); 
% 格兰杰因果需要工具箱 MVGC: https://github.com/lcbarnett/MVGC1
% Barnett, L., & Seth, A. K. (2015). Granger causality for state-space models. Physical Review E, 91(4), 040101.
Mutual_information_all = cell(6,4); 
Transfer_entropy_A_all = cell(6,4); 
Transfer_entropy_B_all = cell(6,4); 
% 互信息和传输熵需要工具箱 JIDT：https://github.com/jlizier/jidt
% Lizier, J. T. (2014). JIDT: An information-theoretic toolkit for studying the dynamics of complex systems. Frontiers in Robotics and AI, 1, 11.
% 同时需要安装Java
javaaddpath('D:\Matlab2024a\matlab_jidt\infodynamics.jar');

for regionA = 1: 3
    for regionB = regionA+1: 4

        z = z+1;
        trial_number = cell(1, 4);
        for i_type = 1: 4
            trial_number{i_type} = intersect( trial_all{regionA, i_type}, trial_all{regionB, i_type} );
        end
        % 根据 trial number 对应获取 regionA/B 在 trial_all / CD_proj_all 的索引
        Index_A = cell(4, 1);
        Index_B = cell(4, 1);
        for i = 1:4
            Index_A{i} = find(ismember(trial_all{regionA, i}, trial_number{i}));
            Index_B{i} = find(ismember(trial_all{regionB, i}, trial_number{i}));
        end

        Pearson_correlation = cell(1, 4);
        Cross_correlation = cell(1, 4);
        Granger_causality_A = cell(1, 4);
        Granger_causality_B = cell(1, 4);
        Mutual_information = cell(1, 4);
        Transfer_entropy_A = cell(1, 4);
        Transfer_entropy_B = cell(1, 4);

        for i_condition = 1:4
            CD_proj_selectedA = CD_proj_all{regionA, i_condition}(:, Index_A{i_condition});
            CD_proj_selectedB = CD_proj_all{regionB, i_condition}(:, Index_B{i_condition});
            if size(CD_proj_selectedA) ~= size(CD_proj_selectedB)
                error('两个矩阵的大小必须相等，Timepoints × Trials');
            end
            trial_number = size(CD_proj_selectedA, 2);

            Pearson_correlation{i_condition} = NaN(trial_number, 1); 
            Cross_correlation{i_condition} = []; 
            Granger_causality_A{i_condition} = NaN(trial_number, 1);  
            Granger_causality_B{i_condition} = NaN(trial_number, 1);  
            Mutual_information{i_condition} = NaN(trial_number, 1); 
            Transfer_entropy_A{i_condition} = NaN(trial_number, 1);  
            Transfer_entropy_B{i_condition} = NaN(trial_number, 1);  

            for func_conn = 1:trial_number
                % Pearson correlation
                Pearson_correlation{i_condition} (func_conn) = corr(CD_proj_selectedA(:, func_conn), CD_proj_selectedB(:, func_conn), 'Type', 'Pearson');
                
                % Cross correlation (auto_lag)
                [cc, lags] = xcorr(CD_proj_selectedA(:, func_conn), CD_proj_selectedB(:, func_conn), 'coeff');
                Cross_correlation{i_condition}(func_conn,:) = cc; 

                % Granger Causality
                GCmatrix = [CD_proj_selectedA(:, func_conn), CD_proj_selectedB(:, func_conn)];
                % 非平稳信号采用状态空间方法；根据BIC准则选择最优滞后阶数
                [AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(GCmatrix', 20, 'LWR', false);  % Model order estimation, 'LWR' or 'OLS'
                [A, SIG] = tsdata_to_var(GCmatrix', moBIC, 'LWR'); % VAR model estimation, or 'actual', 'moAIC', numerical value)
                [F, Fsig] = var_to_pwcgc(A, SIG, GCmatrix', 'LWR', 'F'); % select GC value or p-value 
                Granger_causality_A{i_condition} (func_conn) = F(1, 2); 
                Granger_causality_B{i_condition} (func_conn) = F(2, 1); % effective connectivity
                 
                % Mutual Information
                miCalc = javaObject('infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov1');
                miCalc.initialise(1, 1);  
                miCalc.setObservations(octaveToJavaDoubleArray(CD_proj_selectedA(:, func_conn)), octaveToJavaDoubleArray(CD_proj_selectedB(:, func_conn)));
                Mutual_information{i_condition} (func_conn) = miCalc.computeAverageLocalOfObservations();
                 
                % Transfer entropy
                teCalc = javaObject('infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorKraskov');
                teCalc.initialise(1);  % 设置历史长度为 1 (Schreiber k=1)
                teCalc.setProperty('k', '4');  % 设置 Kraskov 方法的参数 K = 4（使用 4 个最近邻点）
                teCalc.setProperty("NOISE_LEVEL_TO_ADD", "0");
                teCalc.setObservations(CD_proj_selectedA(:, func_conn), CD_proj_selectedB(:, func_conn));
                Transfer_entropy_A{i_condition} (func_conn) = teCalc.computeAverageLocalOfObservations();
                teCalc.setObservations(CD_proj_selectedB(:, func_conn), CD_proj_selectedA(:, func_conn));
                Transfer_entropy_B{i_condition} (func_conn) = teCalc.computeAverageLocalOfObservations();

            end
        end
        Pearson_correlation_all{z, 1} = Pearson_correlation{1};  
        Pearson_correlation_all{z, 2} = Pearson_correlation{2};  
        Pearson_correlation_all{z, 3} = Pearson_correlation{3};  
        Pearson_correlation_all{z, 4} = Pearson_correlation{4}; 
        Cross_correlation_all{z, 1} = Cross_correlation{1};  
        Cross_correlation_all{z, 2} = Cross_correlation{2};  
        Cross_correlation_all{z, 3} = Cross_correlation{3};  
        Cross_correlation_all{z, 4} = Cross_correlation{4}; 
        Granger_causality_A_all{z, 1} = Granger_causality_A{1};  
        Granger_causality_A_all{z, 2} = Granger_causality_A{2};  
        Granger_causality_A_all{z, 3} = Granger_causality_A{3};  
        Granger_causality_A_all{z, 4} = Granger_causality_A{4}; 
        Granger_causality_B_all{z, 1} = Granger_causality_B{1};  
        Granger_causality_B_all{z, 2} = Granger_causality_B{2};  
        Granger_causality_B_all{z, 3} = Granger_causality_B{3};  
        Granger_causality_B_all{z, 4} = Granger_causality_B{4}; 
        Mutual_information_all{z, 1} = Mutual_information{1};  
        Mutual_information_all{z, 2} = Mutual_information{2};  
        Mutual_information_all{z, 3} = Mutual_information{3};  
        Mutual_information_all{z, 4} = Mutual_information{4}; 
        Transfer_entropy_A_all{z, 1} = Transfer_entropy_A{1};  
        Transfer_entropy_A_all{z, 2} = Transfer_entropy_A{2};  
        Transfer_entropy_A_all{z, 3} = Transfer_entropy_A{3};  
        Transfer_entropy_A_all{z, 4} = Transfer_entropy_A{4}; 
        Transfer_entropy_B_all{z, 1} = Transfer_entropy_B{1};  
        Transfer_entropy_B_all{z, 2} = Transfer_entropy_B{2};  
        Transfer_entropy_B_all{z, 3} = Transfer_entropy_B{3};  
        Transfer_entropy_B_all{z, 4} = Transfer_entropy_B{4}; 
        
    end
end

%% FC visualization
data_sets = {'Pearson_correlation_all', 'Granger_causality_A_all', 'Granger_causality_B_all', 'Mutual_information_all','Transfer_entropy_A_all', 'Transfer_entropy_B_all'};
titles = {'Pearson Correlation', 'Granger Causality A', 'Granger Causality B', 'Mutual Information', 'Transfer Entropy A', 'Transfer Entropy B'};
subplots_titles1 = {'SC--MRN', 'SC--SNr', 'SC--ALM', 'MRN--SNr', 'MRN--ALM', 'SNr--ALM'};
subplots_titles2 = {'SC-->MRN', 'SC-->SNr', 'SC-->ALM', 'MRN-->SNr', 'MRN-->ALM', 'SNr-->ALM'};
subplots_titles3 = {'SC<--MRN', 'SC<--SNr', 'SC<--ALM', 'MRN<--SNr', 'MRN<--ALM', 'SNr<--ALM'};
subplots_titles = {subplots_titles1, subplots_titles2, subplots_titles3, subplots_titles1, subplots_titles2, subplots_titles3};
ylimitations = {[-0.2 0.2], [0 0.01], [0 0.01], [0 0.55], [-0.02 1.2] , [-0.02 1.2]};
save_path = 'C:\Users\XinHao\Desktop\WM_DualCircuit\data\CIBRZC85\2024_01_15\FCfigures2\';
% subgroups
for d = 1:length(data_sets)
    figure;
    data = eval(data_sets{d});  % 获取当前数据集
    for z = 1:6
        subplot(2, 3, z);
        hold on;
        plotdata = data(z, :);  % 取出当前区域对应的4个条件数据
        all_data = [];
        group_labels = [];
        for i = 1:4
            all_data = [all_data; plotdata{i}(:)];
            group_labels = [group_labels; repmat(i, length(plotdata{i}), 1)];
        end
        x = 1:4;  
        violinplot(all_data, group_labels, 'ShowMean', true);
        ylim(ylimitations{d});
        title(subplots_titles{d}{z});
        xticks(x);
        xticklabels({'YesCorrect', 'NoCorrect', 'YesError', 'NoError'});
        ylabel(titles{d});
        grid on;
        hold off;
    end
    sgtitle(strrep(data_sets{d}, '_', ' '));  
    saveas(gcf, fullfile(save_path, [data_sets{d}, '.png']));  
    close(gcf);  
end

% total trials
for d = 1:length(data_sets)
    figure;
    currentdata = eval(data_sets{d});
    for i = 1:6
        subplot(3, 2, i); hold on;
        all_data = []; 
        for j = 1:4
            data = currentdata(i, j);
            if isempty(data)
                continue;
            end
            all_data = [all_data; data];
        end
        all_data = vertcat(all_data{:});
        violinplot(all_data);
        xticklabels({'Total Trials'});
        ylim(ylimitations{d});
        ylabel(titles{d});
        title(subplots_titles{d}{z});
    end
    sgtitle(strrep(data_sets{d}, '_', ' '));  
    saveas(gcf, fullfile(save_path, [data_sets{d}, '_total.png']));  
    close(gcf);  
end

% cross correlation plots
conditions = {'YesCorrect', 'NoCorrect', 'YesError', 'NoError'};
colors = {'#018A67', '#1868B2', '#F3A332', '#DE582B'}; 
figure;
for i = 1:6
    subplot(3, 2, i); hold on;
    for j = 1:4
        data = Cross_correlation_all{i, j};
        mean_data = mean(data, 1);  
        std_data = std(data, 0, 1);
        plot(mean_data, 'Color', colors{j}, 'LineWidth', 2);
        [max_val, max_idx] = max(mean_data);
        plot(max_idx, max_val, 'o', 'Color', colors{j}, 'MarkerFaceColor', colors{j}, 'MarkerSize', 6); 
        text(max_idx, max_val, sprintf('(%d, %.2f)', max_idx, max_val), ...
             'FontSize', 10, 'Color', colors{j}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
        %fill([1:length(mean_data), fliplr(1:length(mean_data))],[mean_data + std_data, fliplr(mean_data - std_data)],colors{j}, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end
    title(subplots_titles1{i});
    xlabel('Lags (point)');
    ylabel('correlation');
    ylim([-0.2, 0.2]);  
    legend(conditions, 'Location', 'bestoutside', 'FontSize', 10, 'Box', 'off');
end
sgtitle('Cross Correlation');
saveas(gcf, fullfile(save_path, 'crosscorr_allpoint.png'));  
close(gcf);  