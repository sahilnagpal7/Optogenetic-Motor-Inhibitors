% Velocity Plotting code - Part 2 - Optogenetic control of kinesins -1, -2, -3 and dynein reveals their specific roles in vesicular transport

% Positive, Negative and "Neutral" velocity categorization:

clear all; clc; close all;
set(0,'DefaultFigureWindowStyle','docked');
addpath('/Users/sahilnagpal/Desktop/Abdullah-codes-20Nov22/General codes');

% Name = 'Velocity Analysis [AVG Velo - 30 seconds before condition and (+/-) 0.2 um/s Neutral velocity]';
cd('/Volumes/AGHlab-data/NEWWW/23-08-09_kif16b_ovr_opto/Dish1_opto/Cell5/2/');
File = num2str(input("Enter the name of the velocity file you want to analyse?: ",'s'))
load(File);

velocity_2 =   velocity_1(~cellfun('isempty',velocity_1)); % removing all the empty values
Sliding_time =   Slid_time(~cellfun('isempty',Slid_time)); % removing all the empty values
velocity_final  = [];
L = 0;
LastValueTrackID = numel(velocity_2) % getting the last cell number to do for loop

for tn=1:LastValueTrackID
    count = 0;
    if isempty(velocity_2{tn})
        count = count+1; % finding the number of conditions which they have 0 velocities
        velocity_2{tn} = [];
        Sliding_time{tn} = [];
    end
    if count == numel(velocity_2{tn}) % if the number of conditions with zero velocity match with the number of condition itself, marking that velocity cell = 0.
        velocity_2{tn} = [];
        Sliding_time{tn} = [];
    end

    if numel(velocity_2{tn})<=400 % NOTE: Considering only those trajectories having datapoints greater than 400.
        velocity_2{tn}=[];
        Sliding_time{tn} = [];
    end
end


velocity_final =  velocity_2(~cellfun('isempty',velocity_2))'; % removing those empty cells
Sliding_time_processed = Sliding_time(~cellfun('isempty',Sliding_time))';

for tn = 1:numel(velocity_final)
    % This is for the first transition condition (from dark 1 to lit state):
    if numel(find(Sliding_time_processed{tn} == 120.0450)) == 1
        condition_indx = find(Sliding_time_processed{tn} > 85 & Sliding_time_processed{tn} <= 115); %% 120 105
        velocity_mean(tn) = mean(velocity_final{tn}(condition_indx));
        velocity_mean_test1 = velocity_mean;

    end
    clear velocity_mean_2 Sliding_time_processed_1 condition_indx
end

velocity_mean = velocity_mean';
Sliding_time_final = Sliding_time_processed;


%% Categorizing the velocities into positive, negative and "Neutral" velocity
for pk = 1:numel(velocity_mean)
    if velocity_mean(pk) > 0.01 % Positive velocity - change the values here with 0.01 increments
        velocity_values = velocity_final{pk};
        Slid_time_values = Sliding_time_final{pk};
        velocity_positive{pk}=velocity_values;
        Sliding_time_positive{pk} = Slid_time_values; 
        velocity_negative{pk} = {};
        Sliding_time_negative{pk} = {};
        velocity_neut{pk} = {};
        Sliding_time_neut{pk} = {};
        velo_pos(pk) = velocity_mean(pk);
        velo_neg(pk)=nan;
        velo_neut(pk)=nan;

    elseif velocity_mean(pk) < -0.01 % Negative velocity - change the values here with 0.01 increments
        velocity_values = velocity_final{pk};
        Slid_time_values = Sliding_time_final{pk};
        velocity_negative{pk}=velocity_values;
        Sliding_time_negative{pk} = Slid_time_values; 
        velocity_positive{pk} = {};
        Sliding_time_positive{pk} = {};
        velocity_neut{pk} = {};
        Sliding_time_neut{pk} = {};
        velo_pos(pk) = nan;
        velo_neg(pk)=velocity_mean(pk);
        velo_neut(pk)=nan;

    else % Neutral velocity
        velocity_values = velocity_final{pk};
        Slid_time_values = Sliding_time_final{pk};
        velocity_neut{pk}=velocity_values;
        Sliding_time_neut{pk} = Slid_time_values; 
        velocity_positive{pk} = {};
        Sliding_time_positive{pk} = {};
        velocity_negative{pk} = {};
        Sliding_time_negative{pk} = {};
        velo_pos(pk) = nan;
        velo_neut(pk)=velocity_mean(pk);
        velo_neg(pk)=nan;
    end
end


%% Removing the empty cells
velocity_positive =   velocity_positive(~cellfun('isempty',velocity_positive));
velocity_negative =   velocity_negative(~cellfun('isempty',velocity_negative));
velocity_neut =   velocity_neut(~cellfun('isempty',velocity_neut));
Sliding_time_positive =   Sliding_time_positive(~cellfun('isempty',Sliding_time_positive));
Sliding_time_negative =   Sliding_time_negative(~cellfun('isempty',Sliding_time_negative));
Sliding_time_neut =   Sliding_time_neut(~cellfun('isempty',Sliding_time_neut));

Slid_time_cat = cat(1,Sliding_time_processed{:});
min_slid_time = min(Slid_time_cat);
max_slid_time = max(Slid_time_cat);
slid_time_ref = min_slid_time:0.302:max_slid_time;
Velo_average =[];

%% Positive average velocity calculation:
for ind = 1:numel(slid_time_ref) % Reference time - 1:1262
    check_1 = slid_time_ref(ind);
    ind;
    for kar = 1:numel(Sliding_time_positive)% running through each trajectory 
        kar;
        for pol = 1:numel(Sliding_time_positive{kar}); % running through each values of that trajectory (Varies)
            pol;
            if check_1 == Sliding_time_positive{kar}(pol);
                Velo_average = velocity_positive{kar}(pol);
                velo = Velo_average;
                Velo_average1(kar)=velo;
            else
                Velo_average(kar) = nan;
            end
        end
    end
    Velo_average2_pos{ind}=Velo_average1;
    Velo_average3_pos(ind) = mean(Velo_average1);
    Velo_average3_pos_STD(ind)=std(Velo_average1); % adding STDev to the average plot
end

clear Velo_average1 Velo_average ind kar pol;

%% Negative average velocity calculation:
for ind = 1:numel(slid_time_ref) 
    check_1 = slid_time_ref(ind);
    ind;
    for kar = 1:numel(Sliding_time_negative)% running through each trajectory 
        kar;
        for pol = 1:numel(Sliding_time_negative{kar}) % running through each values of that trajectory (Varies)
            pol;
            if check_1 == Sliding_time_negative{kar}(pol)
                Velo_average = velocity_negative{kar}(pol);
                velo = Velo_average;
                Velo_average1(kar)=velo;
            else
                Velo_average(kar) = nan;
            end
        end
    end
    Velo_average2_neg{ind}=Velo_average1;
    Velo_average3_neg(ind) = mean(Velo_average1);
    Velo_average3_neg_STD(ind)=std(Velo_average1); % adding STDev to the average plot
end


%% Neutral average velocity calculation
for ind = 1:numel(slid_time_ref) % Reference time 
    check_1 = slid_time_ref(ind);
    ind;
    for kar = 1:numel(Sliding_time_neut)% running through each trajectory - 1:21
        kar;
        for pol = 1:numel(Sliding_time_neut{kar}) % running through each values of that trajectory (Varies)
            pol;
            if check_1 == Sliding_time_neut{kar}(pol)
                Velo_average = velocity_neut{kar}(pol);
                velo = Velo_average;
                Velo_average1(kar)=velo;
            else
                Velo_average(kar) = nan;
            end
        end
    end
    Velo_average2_neut{ind}=Velo_average1;
    Velo_average3_neut(ind) = mean(Velo_average1);
    Velo_average3_neut_STD(ind)=std(Velo_average1); % adding STDev to the average plot
end

%% Percentage
velo_percent_total = numel(velocity_neut) + numel(velocity_negative) + numel(velocity_positive);
velo_perc_pos = numel(velocity_positive)/velo_percent_total;
velo_perc_neg = numel(velocity_negative)/velo_percent_total;
velo_perc_neut = numel(velocity_neut)/velo_percent_total;

%% FIGURE 1 - All trajectories + average:
graylight = [.7 .7 .7];
light_red = [211 94 96]./255;
light_blue = [114 147 203]./255;
light_green = [132 186 91]./255;
blue = [57 106 177]./255;
red = [204 37 41]./255;
green = [62 150 81]./255;

figure(1); hold on; xlabel('Sliding time (s)'); ylabel('Velocity (\mum/s)'); %title(Name);
for tn=1:numel(velocity_positive) % Positive
    tn
    plot(Sliding_time_positive{tn}-120, velocity_positive{tn},'Color', 'g', 'LineWidth',0.2);
end
for tn=1:numel(velocity_negative) % Negative
    tn
    plot(Sliding_time_negative{tn}-120, velocity_negative{tn}, 'Color', 'm', 'LineWidth',0.2);
end
for tn=1:numel(velocity_neut) % Inbetween
    tn
    plot(Sliding_time_neut{tn}-120, velocity_neut{tn}, 'Color', graylight, 'LineWidth',0.2);
end
legend('on')
xline([115.47 - 120], '-', {'Dark 1 to Lit transition'}); % These numbers need to be flexible!
xline([250.47 - 120], '-', {'Lit to Dark 2 transition'});
plot(slid_time_ref-120, Velo_average3_pos, 'Color', 'g', 'LineWidth',2.5, 'DisplayName','Positive Velocity');
plot(slid_time_ref-120, Velo_average3_neg, 'Color', 'm', 'LineWidth',2.5, 'DisplayName','Negative Velocity');
plot(slid_time_ref-120, Velo_average3_neut, 'Color', graylight, 'LineWidth',2.5, 'DisplayName','Neutral Velocity');
publication_fig(0,0,1);
hold off;


%% FIGURE 2 - Average velocity plot
figure(2); hold on; xlabel('Sliding time (s)'); ylabel('Velocity (\mum/s)'); % title(Name);
xticks(-120:50:270);
publication_fig(0,0,1);

% if Quest == 1 % For this program, put the lower time series as 0, 135 and 270 and upper time series as 120, 255 and 390
disp('How many conditions do you want to consider?');
Transition_no = input('Enter a number (either 1, 2 or 3): ');
for tn = 1:Transition_no
    Time_S_Trans1 = input('Enter lower time series: ');
    Time_S_Trans2 = input('Enter upper time series: ');
    Time_Pul_Trans{tn} = [Time_S_Trans1; Time_S_Trans2];
    times_pulse{tn} = Time_Pul_Trans{tn};
end

kpr = 0;
for tn = 1:numel(times_pulse) % Checking number of columns
    kpr=kpr+1;
    tp1=times_pulse{tn}; % i has 2 row values (as it is in cell format) --> 1st row: Time in which the condition started (eg: 0, 135, 270) and 2nd Row: Time in which the condition ended (eg: 120, 255, 390)
    Velo_average3_pos_seperate{kpr}=Velo_average3_pos(slid_time_ref>=tp1(1) & slid_time_ref<tp1(2)); % Finding positions that are between the condition starting state and condition ending state
    Velo_average3_pos_seperate_STD{kpr}=Velo_average3_pos_STD(slid_time_ref>=tp1(1) & slid_time_ref<tp1(2));
    slid_time_ref_seperate{kpr}=slid_time_ref(slid_time_ref>=tp1(1) & slid_time_ref<tp1(2)); % Finding time frames that are between the condition starting state and condition ending state
    Velo_average3_neg_seperate{kpr}=Velo_average3_neg(slid_time_ref>=tp1(1) & slid_time_ref<tp1(2)); % Finding positions that are between the condition starting state and condition ending state
    Velo_average3_neg_seperate_STD{kpr}=Velo_average3_neg_STD(slid_time_ref>=tp1(1) & slid_time_ref<tp1(2));
    Velo_average3_neut_seperate{kpr}=Velo_average3_neut(slid_time_ref>=tp1(1) & slid_time_ref<tp1(2)); % Finding positions that are between the condition starting state and condition ending state
    Velo_average3_neut_seperate_STD{kpr}=Velo_average3_neut_STD(slid_time_ref>=tp1(1) & slid_time_ref<tp1(2));
end


graydark  = [0.4 0.4 0.4];
graylight = [.7 .7 .7];
blacklight  = [0.2 0.2 0.2];
light_green = [132 186 91]./255;
light_blue = [114 147 203]./255;

for kfr = 1:numel(Velo_average3_pos_seperate)
    if kfr/2 ~= 1
        plot(slid_time_ref_seperate{kfr}-120, Velo_average3_pos_seperate{kfr}, 'g', 'LineStyle', '-', 'LineWidth', 4, 'DisplayName','Positive Velocity');
        plot(slid_time_ref_seperate{kfr}-120, Velo_average3_neg_seperate{kfr}, 'Color', 'm','LineStyle', '-', 'LineWidth', 4, 'DisplayName','Negative Velocity');
        plot(slid_time_ref_seperate{kfr}-120, Velo_average3_neut_seperate{kfr}, 'Color', graylight,'LineStyle', '-', 'LineWidth', 4, 'DisplayName','Neutral Velocity');
    elseif kfr/2 == 1
        plot(slid_time_ref_seperate{kfr}-120, Velo_average3_pos_seperate{kfr}, 'g', 'LineStyle', '-', 'LineWidth', 4);
        plot(slid_time_ref_seperate{kfr}-120, Velo_average3_neg_seperate{kfr}, 'Color', 'm','LineStyle', '-', 'LineWidth', 4);
        plot(slid_time_ref_seperate{kfr}-120, Velo_average3_neut_seperate{kfr}, 'Color', graylight,'LineStyle', '-', 'LineWidth', 4);
        set(findall(gca, 'Type', 'Line'),'LineWidth',5);
    end

end
ylim([-0.32, 0.32]);
hold on
plot([0 0], get(gca, 'ylim'), 'k', 'LineWidth', 4);
hold off;



%% FIGURE 3 - Average velocity plot with STD
figure(3); hold on; xlabel('Sliding time (s)'); ylabel('Velocity (\mum/s)'); %title(Name);
xticks(-130:10:270);


kpr = 0;
for tn = 1:numel(times_pulse) % Checking number of columns
    kpr=kpr+1;
    tp1=times_pulse{tn}; % i has 2 row values (as it is in cell format) --> 1st row: Time in which the condition started (eg: 0, 135, 270) and 2nd Row: Time in which the condition ended (eg: 120, 255, 390)
    Velo_average3_pos_seperate{kpr}=Velo_average3_pos(slid_time_ref>=tp1(1) & slid_time_ref<tp1(2)); % Finding positions that are between the condition starting state and condition ending state
    Velo_average3_pos_seperate_STD{kpr}=Velo_average3_pos_STD(slid_time_ref>=tp1(1) & slid_time_ref<tp1(2));
    slid_time_ref_seperate{kpr}=slid_time_ref(slid_time_ref>=tp1(1) & slid_time_ref<tp1(2)); % Finding time frames that are between the condition starting state and condition ending state
    Velo_average3_neg_seperate{kpr}=Velo_average3_neg(slid_time_ref>=tp1(1) & slid_time_ref<tp1(2)); % Finding positions that are between the condition starting state and condition ending state
    Velo_average3_neg_seperate_STD{kpr}=Velo_average3_neg_STD(slid_time_ref>=tp1(1) & slid_time_ref<tp1(2));
    Velo_average3_neut_seperate{kpr}=Velo_average3_neut(slid_time_ref>=tp1(1) & slid_time_ref<tp1(2)); % Finding positions that are between the condition starting state and condition ending state
    Velo_average3_neut_seperate_STD{kpr}=Velo_average3_neut_STD(slid_time_ref>=tp1(1) & slid_time_ref<tp1(2));
end


graydark  = [0.4 0.4 0.4];
graylight = [.7 .7 .7];
blacklight  = [0.2 0.2 0.2];
light_green = [132 186 91]./255;
light_blue = [114 147 203]./255;

for kfr = 1:numel(Velo_average3_pos_seperate)
    if kfr/2 ~= 1
        legend('on');
        errorbar(slid_time_ref_seperate{kfr}-120, Velo_average3_pos_seperate{kfr}, Velo_average3_pos_seperate_STD{kfr}, 'color', 'g', 'LineWidth', 0.02); % adding STDev to the average plot
        errorbar(slid_time_ref_seperate{kfr}-120, Velo_average3_neg_seperate{kfr}, Velo_average3_neg_seperate_STD{kfr}, 'color', 'm', 'LineWidth', 0.02); % adding STDev to the average plot
        errorbar(slid_time_ref_seperate{kfr}-120, Velo_average3_neut_seperate{kfr}, Velo_average3_neut_seperate_STD{kfr}, 'color', graylight, 'LineWidth', 0.02);
        
        plot(slid_time_ref_seperate{kfr}-120, Velo_average3_pos_seperate{kfr}, '-g', 'LineWidth',2, 'DisplayName','Positive Velocity');
        plot(slid_time_ref_seperate{kfr}-120, Velo_average3_neg_seperate{kfr}, 'Color', 'm', 'LineWidth',2, 'DisplayName','Negative Velocity');
        plot(slid_time_ref_seperate{kfr}-120, Velo_average3_neut_seperate{kfr}, 'Color', graylight, 'LineWidth',2, 'DisplayName','Neutral Velocity');
    elseif kfr/2 == 1
        errorbar(slid_time_ref_seperate{kfr}-120, Velo_average3_pos_seperate{kfr}, Velo_average3_pos_seperate_STD{kfr}, 'color', 'g', 'LineWidth', 0.02); % adding STDev to the average plot
        errorbar(slid_time_ref_seperate{kfr}-120, Velo_average3_neg_seperate{kfr}, Velo_average3_neg_seperate_STD{kfr}, 'color', 'm', 'LineWidth', 0.02); % adding STDev to the average plot
        errorbar(slid_time_ref_seperate{kfr}-120, Velo_average3_neut_seperate{kfr}, Velo_average3_neut_seperate_STD{kfr}, 'color', graylight, 'LineWidth', 0.02);
        plot(slid_time_ref_seperate{kfr}-120, Velo_average3_pos_seperate{kfr}, '-g', 'LineWidth',2);
        plot(slid_time_ref_seperate{kfr}-120, Velo_average3_neg_seperate{kfr}, '-m', 'LineWidth',2);
        plot(slid_time_ref_seperate{kfr}-120, Velo_average3_neut_seperate{kfr}, 'Color', graylight, 'LineWidth',2);
    end

end
xline([115.47 - 120], '-', {'Dark 1 to Lit transition'});
xline([250.47 - 120], '-', {'Lit to Dark 2 transition'});

publication_fig(0,0,1);
hold off;
