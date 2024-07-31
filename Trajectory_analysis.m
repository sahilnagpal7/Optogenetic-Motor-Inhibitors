% Cargo Trajectory Analysis - Optogenetic control of kinesins -1, -2, -3 and dynein reveals their specific roles in vesicular transport -Adapted from code originally written by former graduate student in our lab Linda Balabanian

% Cumulates cell means within condition. 
% 4 conditions to select for tracks. 
%Takes xlm file generated from ImageJ plugin - TrackMate - Modify the file
%name with x and y cordinates of cell center (followed by 'x', 'y') and 
% area of cell (followed by 'a')

clear all, warning off, close all
set(0,'DefaultFigureWindowStyle','docked')
addpath('/Users/sahilnagpal/Desktop/Desktop/Fiji/Fiji.app/scripts/'); %Modify this

% Call on all files with the desired condition
cells = dir('*_Tracks*.xml');

for d=1:length(cells)
    cell=cells(d, :);  % calling initially on 1st cell of this condition, and then so on.
    
file_path_tracks = cell.name;
clipZ = true;     % Remove Z coordinates, if you know you can.
scaleT = true;       % Use physical time for T.
tracks = importTrackMateTracks(file_path_tracks, clipZ, scaleT);
n_tracks = numel( tracks );

% (1) FIRST, calculate Rg of tracks like I normally would, but importantly,
% doing this before having cell center as 0, 0.
new_tracks = [];
for s = 1 : n_tracks
    if (numel(tracks{s}(:, 1)) > 5) && (tracks{s}(1, 1) == 0) && (tracks{s}(1, 2) ~= tracks{s}(2, 2)) && (tracks{s}(1, 3) ~= tracks{s}(2, 3))
        % 1st condition: Select only tracks with minimum number of points, now 5
        % 2nd condition: Select only trajectories of the cell that commence in the 1st frame of movie
        % 3rd and 4th condition: Getting rid of tracks that are not true trajectories, not moving at all in x and y
        new_tracks1 = tracks{s}(:, 1);
        i = find(new_tracks1 < 120);  % To cut movie at certain frame.
        new_tracks1 = new_tracks1(i);
        
        new_tracks2 = tracks{s}(:, 2);
        new_tracks2 = new_tracks2(i);
        
        new_tracks3 = tracks{s}(:, 3);
        new_tracks3 = new_tracks3(i);
        
        new_tracks_final =[new_tracks1,new_tracks2,new_tracks3];
        new_tracks{s} = new_tracks_final;
    end
end
new_tracks = new_tracks';
tracks = new_tracks;
tracks(cellfun('isempty', tracks)) = [];   % getting rid of all cells without values
n_tracks = numel( tracks );        % number of tracks


% Radius of gyration Rg for each track
all_x = [];
all_y = [];
all_a = [];
all_b = [];
all_c = [];
mean_x_cum = zeros(n_tracks, 1);
mean_y_cum = zeros(n_tracks, 1);
Rg = zeros(n_tracks, 1);
for s = 1 : n_tracks
    % x and y coordinates of every spot in a track
    all_x{s} = tracks{s}(:, 2);
    all_y{s} = tracks{s}(:, 3);
    mean_x = mean(all_x{s}(:, 1));
    mean_y = mean(all_y{s}(:, 1));
    % mean of x and y coordinates for each track
    mean_x_cum(s) = mean_x;
    mean_y_cum(s) = mean_y;
    % calculation of a = (x-X)^2 and b = (y-Y)^2 for each x and y
    for x_value = 1 : numel(all_x{s}(:, 1))
        a = ((all_x{s}(x_value, 1))-mean_x).^2;
        b = ((all_y{s}(x_value, 1))-mean_y).^2;
        c = a+b;
        all_a{s}(x_value) = a;
        all_b{s}(x_value) = b;
        all_c{s}(x_value) = c;   % c values then are all summed for a track
    end
    Rg(s) = sqrt((1/numel(all_x{s}(:, 1)))*sum(all_c{s}(1, :))); % Rg for a track
end
% You can make a histogram of hist(Rg) to see the distribution for Rg's of a cell


% (2) SECOND, re-taking the tracks but rendering their position relative to the
% cell center as coordinates 0, 0.

% Extracting the x and y coordinate values of the MTOC center from the file name. 
x_center = str2double(regexp(cell.name, '\d+\.?\d+(?=x)', 'match', 'once')); % This means that if I put an 'x' after the x coordinate value in the filename, it will be selected.
y_center = str2double(regexp(cell.name, '\d+\.?\d+(?=y)', 'match', 'once'));
% Extracting the cell area value from the file name.
cell_area = str2double(regexp(cell.name, '\d+\.?\d+(?=a)', 'match', 'once'));
% From the cell area value extracted above, we calculate the radius of the
% circle (approximate as the cells are not really circular) using the
% formula Area = pi*radius^2
cell_radius = sqrt(cell_area/pi);

new_tracks = [];
for s = 1 : n_tracks
    if (numel(tracks{s}(:, 1)) > 5) && (tracks{s}(1, 1) == 0) && (tracks{s}(1, 2) ~= tracks{s}(2, 2)) && (tracks{s}(1, 3) ~= tracks{s}(2, 3))
        % 1st condition: Select only tracks with minimum number of points, now 5
        % 2nd condition: Select only trajectories of the cell that commence in the 1st frame of movie
        % 3rd and 4th condition: Getting rid of tracks that are not true trajectories, not moving at all in x and y
        new_tracks1 = tracks{s}(:, 1);
        i = find(new_tracks1 < 120); %To cut movie at certain frame.
        new_tracks1 = new_tracks1(i);
        
        new_tracks2 = tracks{s}(:, 2) - x_center;
        new_tracks2 = new_tracks2(i);
        
        new_tracks3 = tracks{s}(:, 3) - y_center;
        new_tracks3 = new_tracks3(i);
        
        new_tracks_final =[new_tracks1,new_tracks2,new_tracks3];
        new_tracks{s} = new_tracks_final;
    end
end
new_tracks = new_tracks';
tracks = new_tracks;
tracks(cellfun('isempty', tracks)) = [];   % getting rid of all cells without values
n_tracks = numel( tracks );        % number of tracks


% Converts cartesian coordinates into polar coordinates
% Using this to measure radius from cell center (MTOC)
% And plotting it out (not yet here)
theta = [];
rho = [];
for s = 1: n_tracks
[theta{s},rho{s}] = cart2pol(tracks{s}(:, 2), tracks{s}(:, 3));  % theta and rho of each spot of track
end
rho = rho';
% And plotting it out 
polar_plot=figure
polaraxes
hold on
for s = 1: n_tracks
    polarplot(theta{s}, rho{s})
end

FigName = strcat('PolarPlot_cell',num2str(d));
savefig(polar_plot, FigName);

% Taking the mean of rhos for each track
mean_rho = zeros(n_tracks, 1);
% tau_int = zeros(n_tracks, 1);
for s = 1 : n_tracks
mean_rho(s) = mean(rho{s}); % mean of the rhos along a track (rho of each spot within track averaged)
end

% Normalizing the mean rho to the cell radius
mean_rho_norm = mean_rho/cell_radius;


%Rg Categorisation

Low_Rg= [];
Medium_Rg= [];
High_Rg= [];

for s = 1 : n_tracks
    if Rg(s) < 0.3
        Low_Rg = vertcat(Low_Rg, Rg(s));
    elseif ((Rg(s) > 0.3) && (Rg(s) < 0.5))
        Medium_Rg = vertcat(Medium_Rg, Rg(s));
    elseif Rg(s) > 0.5
        High_Rg = vertcat(High_Rg, Rg(s));
    end
end

%  Low_Rg= find(Rg(s) < 0.3);
n_Low_Rg= numel(Low_Rg);
% 
% % Low_Rg= find((Rg(s) > 0.3) && (Rg(s) < 0.6));
n_Medium_Rg= numel(Medium_Rg);
% 
% % High_Rg= find(Rg(s) > 0.6);
n_High_Rg= numel(High_Rg);

perc_Low_Rg(d) = (n_Low_Rg/n_tracks)*100;
perc_Medium_Rg(d) = (n_Medium_Rg/n_tracks)*100;
perc_High_Rg(d) = (n_High_Rg/n_tracks)*100;

% Plotting colour-coded tracks based on Rg
% Low Rg: Black
% Medium Rg: Yellow
% High Rg: Purple
Colored_Rg=figure     
hold on
for s = 1 : n_tracks
    if Rg(s) < 0.3
        x = tracks{s}(:, 2);
        y = tracks{s}(:, 3);
        plot(x, y, '.-', 'Color', [0 0 0], 'LineWidth', 3, 'MarkerSize', 10)
    elseif ((Rg(s) > 0.3) && (Rg(s) < 0.5))
        x = tracks{s}(:, 2);
        y = tracks{s}(:, 3);
        plot(x, y, '.-', 'Color', [1 1 0], 'LineWidth', 3, 'MarkerSize', 10)
    elseif Rg(s) > 0.5
        x = tracks{s}(:, 2);
        y = tracks{s}(:, 3);
        plot(x, y, '.-', 'Color', [0.6 0 0.6], 'LineWidth', 3, 'MarkerSize', 10)
    end
end
axis equal
FigName = strcat('ColoredRg_cell',num2str(d));
savefig(Colored_Rg, FigName);

% (3) THIRD, bringing the Rg of trajectory (1) and mean rho of trajectory
% (2) together in a table for each trajectory of the cell.
%table_cell = [Rg mean_rho];   % Take this one if you don't want to normalize rho mean
% OR
table_cell = [Rg mean_rho_norm]; % Take this one for normalized rho mean

% (4) FINALLY, cumulating the table (3) across all cells within the same
% condition.
number = numel(mean_rho_norm); % number of tracks per cell

% (5) Also, for plotting the mean of rho means of tracks of cell (by mean, here, I mean the mean of rhos of each
% spot within track and then the mean of that for all cell's tracks).
% Collecting all mean of rho means of cells within condition.
mean_rhoMeans = mean(mean_rho_norm);


% To look at DIRECTIONALITY of each trajectory (inward vs outward net motion), 
% use the polar plot rho values of the trajectory 
% start points rho r1 and trajectorys end points rho r2.
% If r2 - r1 is positive, its outward motion. If r2 - r1 is negative, its inward motion.

outward = [];
outward2 = [];
inward = [];
inward2 = [];
stationary = []; 
stationary2 = [];
for s = 1 : n_tracks
 if (Rg(s) >= 0.2) && (rho{s}(1, 1) < rho{s}(end, 1))
  outward = vertcat(outward, Rg(s));   % saving the Rg
  outward2 = vertcat(outward2, mean_rho_norm(s));   % saving the mean_rho_norm
 elseif (Rg(s) >= 0.2) && (rho{s}(1, 1) > rho{s}(end, 1))
  inward = vertcat(inward, Rg(s));
  inward2 = vertcat(inward2, mean_rho_norm(s));
 elseif Rg(s) < 0.2
  stationary = vertcat(stationary, Rg(s));
  stationary2 = vertcat(stationary2, mean_rho_norm(s));
 end
end

% Plotting colour-coded tracks based on inward vs outward motion vs stationary
% Outward motion: BLUE 
% Inward motion: GREEN
% Stationary: GREY
directionality_fig=figure     
hold on
for s = 1 : n_tracks
    if (Rg(s) >= 0.2) && (rho{s}(1, 1) < rho{s}(end, 1))
        x = tracks{s}(:, 2);
        y = tracks{s}(:, 3);
        plot(x, y, '.-', 'Color', [0 0.45 0.70], 'LineWidth', 3, 'MarkerSize', 10)
    elseif (Rg(s) >= 0.2) && (rho{s}(1, 1) > rho{s}(end, 1))
        x = tracks{s}(:, 2);
        y = tracks{s}(:, 3);
        plot(x, y, '.-', 'Color', [0 0.8 0.2], 'LineWidth', 3, 'MarkerSize', 10)
    elseif Rg(s) < 0.2
        x = tracks{s}(:, 2);
        y = tracks{s}(:, 3);
        plot(x, y, '.-', 'Color', [0.8 0.8 0.8], 'LineWidth', 3, 'MarkerSize', 10)
         
    end
   
end
axis equal
FigName = strcat('Directionality_cell',num2str(d));
savefig(directionality_fig, FigName);



%% Calculating and saving the percentage for peripheral, juxtanuclear and perinuclear for
% each type of motion, i.e. outward vs inward vs stationary, for each cell

% 1) outward
out_peri = find(outward2 >= 0.85);
n_out_peri = numel(out_peri); % number of outward runs that are peripheral for a cell
out_jux = find(outward2 >= 0.5);
outward2_new = outward2(out_jux);
out_jux_new = find(outward2_new < 0.85);
n_out_jux = numel(out_jux_new);
out_perinuc = find(outward2 < 0.5);
n_out_perinuc = numel(out_perinuc);

% % 2) inward
in_peri = find(inward2 >= 0.85);
n_in_peri = numel(in_peri);
in_jux = find(inward2 >= 0.5);
inward2_new = inward2(in_jux);
in_jux_new = find(inward2_new < 0.85);
n_in_jux = numel(in_jux_new);
in_perinuc = find(inward2 < 0.5);
n_in_perinuc = numel(in_perinuc);

% % 3) stationary
sta_peri = find(stationary2 >= 0.85);
n_sta_peri = numel(sta_peri);
sta_jux = find(stationary2 >= 0.5);
stationary2_new = stationary2(sta_jux);
sta_jux_new = find(stationary2_new < 0.85);
n_sta_jux = numel(sta_jux_new);
sta_perinuc = find(stationary2 < 0.5);
n_sta_perinuc = numel(sta_perinuc);
%% 
peripheral= n_out_peri + n_in_peri + n_sta_peri;
juxtaneuclear= n_out_jux + n_in_jux + n_sta_jux;
perinuclear= n_out_perinuc + n_in_perinuc + n_sta_perinuc;

percent_perinuclear(d) = (perinuclear/n_tracks)*100;
percent_juxtanuclear(d) = (juxtaneuclear/n_tracks)*100;
percent_peripheral(d) = (peripheral/n_tracks)*100;

% Looking at overall percentages, within cell, of outward vs inward vs stationary
% Saving output for each cell
percent_outward(d) = (numel(outward)/n_tracks)*100;
percent_inward(d) = (numel(inward)/n_tracks)*100;
percent_stationary(d) = (numel(stationary)/n_tracks)*100;


% Cumulating from all cells and saving output
n1 = numel(outward);
n2 = numel(inward);
n3 = numel(stationary);

end







