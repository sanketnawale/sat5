clear all; close all; clc;

%% 🌍 Simulation Parameters
tic
Simulation_T = 110 * 60;  % Total simulation time (110 minutes in seconds)
Time_Step = 60;           % Time step (1 min = 60 sec)
MonteCarlo = 5000;        % Monte Carlo simulation (number of packets per node)
Nodes = 3;                % Number of ground nodes (Rome, Milan, NodeRM)
Pkct_Per_Hour = 10;      % Packets per hour for each node

%% 🌍 Ground Nodes (Rome, Milan, NodeRM)
Node_Coordinates = [ 
    41.9028, 12.4964;  
    45.4642, 9.1900;   
    41.9, 12.5         
];

%% 🛰️ Satellite Constellation (Walker)
Sat_Per_Plane = 18;
Num_Planes = 6;
Total_Sats = Sat_Per_Plane * Num_Planes;
Orbital_Inclination = deg2rad(87);
H = 1200e3;
Earth_Radius = 6378e3;
Time_Vector = 0:Time_Step:Simulation_T;

% 🛰️ Generate Walker Delta Constellation
oev = walker_delta(Sat_Per_Plane, Num_Planes, 1, pi, Earth_Radius + H, Orbital_Inclination);
Num_Satellites = size(oev, 1);
num_steps = length(Time_Vector);

%% 📡 Obtain Satellite Geometry
[Distances, Elevation_Angles, Ground_Distances, Visibility, Num_Visible_Sats, Sat_IDs, Latitudes, Longitudes,Sat_To_NodeRM_Delay] = ...
    Satellite_Geometry(H, Node_Coordinates, oev, Earth_Radius, Time_Vector);
Visible_Sat_Matrix = zeros(length(Time_Vector), 4);

%% Additional storage for packet reception times info
% (For Rome and Milan, we build a string for each time step indicating
% the packet arrival (reception) times per visible satellite.
% For NodeRM, we store a simple message.)
Rome_PktReception = cell(length(Time_Vector), 1);
Milan_PktReception = cell(length(Time_Vector), 1);
NodeRM_PktReception = cell(length(Time_Vector), 1);
% Initialize global Packet ID counter
PacketIDCounter = 0; % This will generate unique IDs for all packets
%% 🌍 LoRa Duty Cycle (Sleep Mode)
Duty_Cycle_Percentage = 1; % Example: Node transmits only 1% of the time
LoRa_Sleep_Time = (100 / Duty_Cycle_Percentage - 1) * 60; % Convert to seconds



%% 🌍 LR-FHSS Communication Parameters
Payload = 100;            
Header_N_DR8 = 3;        
Code_Rate = 1/3;         
Header_duration = 0.233; 
F_duration = 0.05;       
Header_ToA_DR8 = Header_N_DR8 * Header_duration;

% Time on Air Calculation
[ToA_DR8, ToA_DR8_WH] = ToA_Packets_DR8(Payload, Header_ToA_DR8, 2);
ToA_DR8(1) = ToA_DR8(1) + (6.472 / 1000);

% Fragmentation Details
fragment_duration = 50 / 1000;
fragment_50_ms = floor(ToA_DR8_WH(1) / fragment_duration);
Last_fragment_duration = ((ToA_DR8_WH(1) / fragment_duration) - fragment_50_ms) * fragment_duration;
fragment_PHY_length = fragment_50_ms + 1;
fragment_length = Header_N_DR8 + 1 + fragment_PHY_length;

% Define Monte Carlo Simulations
OBW_channels = 280;
Collision_Threshold = 2 * 50e-3;

% 🚀 Define Link Budget Parameters
Tx_Power = 14;  
Antenna_Gain = 2;  
Noise_Floor = -174 + 10 * log10(137e6);  

%% 📊 Initialize Matrices
SuccessRate = zeros(Nodes, length(Time_Vector));
Collisions = zeros(Nodes, length(Time_Vector));
Received_Packets_NodeRM = zeros(1, length(Time_Vector));
Signal_Delay = zeros(Nodes, Num_Satellites, num_steps);


Duty_Cycle = 5;  % 10% duty cycle (adjustable)
sleep_time = (1 - Duty_Cycle) * 60;  % Convert to seconds
wakeup_time = Duty_Cycle * 60;  % Convert to seconds

%% 📡 Main Simulation Loop
for t = 1:length(Time_Vector)
     current_time_min = Time_Vector(t) / 60;
    fprintf('\n⏳ Time %.2f min: \n', current_time_min);
    Visible_Sat_Matrix(t, :) = [current_time_min, Num_Visible_Sats(1, t), Num_Visible_Sats(2, t), Num_Visible_Sats(3, t)];
    NodeRM_Packet_Times = [];
    for n = 1:Nodes-1
        if Num_Visible_Sats(n, t) == 0
            continue;
        end
    % 🚀 **Check if the Node is Sleeping**
        if mod(current_time_min, (sleep_time + wakeup_time) / 60) >= wakeup_time / 60
            fprintf('💤 Node %d is sleeping at %.2f min\n', n, current_time_min);
            continue;  % Skip packet transmission for sleeping node
        end
        

        Num_Packets = 10;
        lambda = 1 / (60 / Pkct_Per_Hour);
        Inter_Arrival_Times = exprnd(1/lambda, 1, Num_Packets);
        Tx_Timestamps = cumsum(Inter_Arrival_Times);
        Tx_Timestamps(Tx_Timestamps > wakeup_time) = [];

        Visible_Sats = Sat_IDs{n, t};
        Sat_Receive_Times = cell(Total_Sats, 1);

        % Packet Transmission Logic for Rome and Milan Nodes
        for pkt = 1:length(Tx_Timestamps)
            % Increment Packet ID Counter
            PacketIDCounter = PacketIDCounter + 1;
            UniquePktID = PacketIDCounter; % Assign a unique packet ID

        for chosen_sat = Visible_Sats
        % Calculate the arrival time at the satellite
            arrival_time = Tx_Timestamps(pkt) + Signal_Delay(n, chosen_sat, t);

        % Store the packet with its unique ID and arrival time
            Sat_Receive_Times{chosen_sat} = [Sat_Receive_Times{chosen_sat}; UniquePktID, arrival_time];
        end

    % Log the packet transmission (optional, for debugging)
        fprintf('📦 Node %d Transmitted Packet ID %d at %.10f seconds\n', n, UniquePktID, Tx_Timestamps(pkt));
    end
%ine Probabilities for Packet Loss
        Header_Loss_Prob = 0.1;   % 10% chance to lose a header
        Fragment_Loss_Prob = 0.15; % 15% chance to lose a fragment

        % 🚀 **Collision Detection & Fragment-Based Tracking**
        target_collided = zeros(1, fragment_length);  % Track collisions per fragment
        target_discarded = zeros(1, fragment_length); % Track discarded packets due to capture effect

        % Simulate random header loss
        for h = 1:Header_N_DR8
            if rand() < Header_Loss_Prob
                target_collided(h) = 1;  % Mark header as lost
            end
        end

% Simulate random fragment loss
        for frag = (Header_N_DR8+2):fragment_length
            if rand() < Fragment_Loss_Prob
                target_collided(frag) = 1;  % Mark fragment as lost
            end
        end

       for s = Visible_Sats
  if ~isempty(Sat_Receive_Times{s})
    % Extract Packet IDs and Arrival Times
    sat_packet_data = Sat_Receive_Times{s}; % Matrix: [PacketID, ArrivalTime]
    sat_packet_data = sortrows(sat_packet_data, 2); % Sort by arrival time
    sat_arrival_times = sat_packet_data(:, 2); % Extract arrival times
    sat_packet_ids = sat_packet_data(:, 1); % Extract Packet IDs

    % 📡 Compute SNR & Apply Rician Fading
    SNR = Tx_Power + Antenna_Gain - Noise_Floor - (20*log10(Distances(n, s, t)/1e3));
    K_factor = 5;
    sigma = sqrt(SNR / (2 * (K_factor + 1))); 
    Fading_SNR = SNR + sigma * randn;
    Decoding_Threshold = 30;  

    % ✅ Detect Collisions
    collisions = sum(diff(sat_arrival_times) < Collision_Threshold);
    total_packets = length(sat_arrival_times);

    % ✅ Display the exact arrival timestamps of packets at this satellite
    formatted_arr = ['[', strtrim(num2str(sat_arrival_times', '%.2f ')), ']'];
    fprintf('⏰ Node %d, Satellite %d arrival packet timings (within %.2f sec): %s\n', ...
            n, s, Time_Step, formatted_arr);

    % ✅ Save arrival times in a formatted string for logging
    pkt_str = sprintf('Sat %d: %s; ', s, formatted_arr);

    if Fading_SNR > Decoding_Threshold
        SuccessRate(n, t) = total_packets - collisions;

        % ✅ If Rome successfully transmits, relay to NodeRM
       if n == 1  % If Rome successfully transmits, relay to NodeRM
    % Extract non-collided packets (packet ID + arrival time)
    non_collided_packets = sat_packet_data(collisions+1:end, :);

    % ✅ Iterate over all non-collided packets
    for pkt_idx = 1:size(non_collided_packets, 1)
        packet_id = non_collided_packets(pkt_idx, 1);
        original_arrival_time = non_collided_packets(pkt_idx, 2);

        % ✅ Iterate over all satellites that received the packet
        for sat_id = Visible_Sats
            % ✅ Compute Total Delay (Node → Satellite + Satellite → NodeRM)
            total_propagation_delay = (Signal_Delay(n, sat_id, t) + Sat_To_NodeRM_Delay(sat_id, t)) / 60; % Convert sec -> min

            % ✅ Compute Final Corrected Arrival Time at NodeRM
            nodeRM_reception_time = original_arrival_time + total_propagation_delay;

            % ✅ Store the packet with the correct reception time
            NodeRM_Packet_Times = [NodeRM_Packet_Times; packet_id, nodeRM_reception_time];

            % ✅ Debugging Output
            fprintf('📡 Node %d → Sat %d → NodeRM | Packet %d | Tx Time: %.6f min | Arrival at NodeRM: %.6f min\n', ...
                n, sat_id, packet_id, original_arrival_time, nodeRM_reception_time);
        end
    end


end

    else
        Collisions(n, t) = collisions;
    end
else
    fprintf('⏰ Node %d, Satellite %d: No packet arrivals during this time step.\n', n, s);
    pkt_str = sprintf('Sat %d: []; ', s);
end

% ✅ Save packet reception logs for Rome & Milan
% Save packet reception logs for Rome & Milan
if n == 1
    Rome_PktReception{t} = [Rome_PktReception{t}; {UniquePktID, arrival_time}];
elseif n == 2
    Milan_PktReception{t} = [Milan_PktReception{t}; {UniquePktID, arrival_time}];
end



        % ✅ Debugging Output
        fprintf('📊 Node %d transmitted %d packets, %d collisions\n', n, Num_Packets, Collisions(n, t));
    end

    % 🚀 **Decoding at NodeRM with Fragment & Header Validation**
   % 🚀 **Decoding at NodeRM with Detailed Failure Analysis**
% 🚀 **Decoding at NodeRM with Detailed Failure Analysis**
if ~isempty(NodeRM_Packet_Times)
    % ✅ Ensure sorting by arrival time
    NodeRM_Packet_Times = sortrows(NodeRM_Packet_Times, 2); % Sort by arrival time
    
    % ✅ Identify visible satellites at NodeRM
   % ✅ Identify visible satellites at NodeRM
NodeRM_Visible_Sats = Sat_IDs{3, t};
common_sats = intersect(NodeRM_Visible_Sats, Sat_IDs{1, t}); % Find common satellites

if ~isempty(common_sats)
    % ✅ Compute NodeRM's SNR
    NodeRM_SNR = Tx_Power + Antenna_Gain - Noise_Floor - (20*log10(Distances(3, common_sats(1), t)/1e3));
    sigma_rm = sqrt(NodeRM_SNR / (2 * (K_factor + 1)));  
    Fading_SNR_RM = NodeRM_SNR + sigma_rm * randn;

    % ✅ Check Headers and Payload Fragments
    Success_header = Header_N_DR8 - length(nonzeros(target_collided(1:Header_N_DR8)));  
    Threshold = fragment_length - round(fragment_PHY_length * (1 - Code_Rate)) - Header_N_DR8 - 1;
    Success_fragment = fragment_length - length(nonzeros(target_collided((Header_N_DR8 + 2):end))) - Header_N_DR8 - 1;

    % ✅ Packet Reception Decision
    if Success_header < 1
        fprintf('❌ NodeRM failed at %.2f min: MISSING HEADERS (Received %d/%d)\n', current_time_min, Success_header, Header_N_DR8);
    elseif Success_fragment < Threshold
        fprintf('❌ NodeRM failed at %.2f min: INSUFFICIENT FRAGMENTS (Received %d/%d)\n', current_time_min, Success_fragment, Threshold);
    elseif Fading_SNR_RM <= Decoding_Threshold
        fprintf('❌ NodeRM failed at %.2f min: LOW SNR (%.2f dB, Threshold = %.2f dB)\n', current_time_min, Fading_SNR_RM, Decoding_Threshold);
    else
        % ✅ Successful Packet Reception at NodeRM
        Received_Packets_NodeRM(t) = size(NodeRM_Packet_Times, 1);
        fprintf('✅ NodeRM successfully received packets at %.2f min:\n', current_time_min);

        % ✅ Adjust the packet arrival time to include propagation delay **before logging**
        sat_to_NodeRM_delay = mean(Signal_Delay(3, common_sats, t)); % Compute avg delay

        for pkt_idx = 1:size(NodeRM_Packet_Times, 1)
            packet_id = NodeRM_Packet_Times(pkt_idx, 1);
            original_arrival_time = NodeRM_Packet_Times(pkt_idx, 2);
            
            % ✅ Compute final arrival time at NodeRM (adding delay)
               % ✅ Corrected arrival time at NodeRM (adding propagation delay + 1 min)
    nodeRM_reception_time = original_arrival_time + total_propagation_delay ; % Convert sec -> min and add 1 min buffer


            % ✅ Store the corrected reception time
            NodeRM_Packet_Times(pkt_idx, 2) = nodeRM_reception_time;

            fprintf('   🛰️ Packet ID: %d, Corrected Arrival Time: %.10f min\n', packet_id, nodeRM_reception_time);
        end

        % ✅ Store received packets for logging
        formatted_packets = strjoin(arrayfun(@(id, time) ...
            sprintf('ID:%d Time:%.2f', id, time), ...
            NodeRM_Packet_Times(:, 1), NodeRM_Packet_Times(:, 2), 'UniformOutput', false), '; ');

        NodeRM_PktReception{t} = sprintf('Received: %s', formatted_packets);
    end
else
    fprintf('❌ NodeRM sees no common satellites with Rome at %.2f min\n', current_time_min);
end

end


    end
end
   

disp(array2table(Visible_Sat_Matrix, 'VariableNames', {'Time_Min','Rome_Sats','Milan_Sats','NodeRM_Sats'}));
% ✅ Define Excel File Path
outputFile = fullfile(pwd, 'Visible_Satellites_Log.xlsx');  % Save in current directory

% ✅ Convert Data to Table
Visible_Sat_Table = array2table(Visible_Sat_Matrix, ...
    'VariableNames', {'Time_Min', 'Rome_Sats', 'Milan_Sats', 'NodeRM_Sats'});

% ✅ Write to Excel File
writetable(Visible_Sat_Table, outputFile, 'Sheet', 'Visibility Data');

% ✅ Confirm Save
fprintf('📄 Visibility data saved to: %s\n', outputFile);
%% Build a detailed table containing visible satellite IDs and packet reception times
Time_Min = round((Time_Vector)' / 60, 1);

% ✅ Convert cell arrays to string format for Excel compatibility
Rome_PktReception = cellfun(@(x) strjoin(string(x), ', '), Rome_PktReception, 'UniformOutput', false);
Milan_PktReception = cellfun(@(x) strjoin(string(x), ', '), Milan_PktReception, 'UniformOutput', false);
NodeRM_PktReception = cellfun(@(x) strjoin(string(x), ', '), NodeRM_PktReception, 'UniformOutput', false);

% ✅ Convert Satellite IDs into readable strings
Rome_Sat_IDs_str = cellfun(@(x) mat2str(x), Sat_IDs(1,:)', 'UniformOutput', false);
Milan_Sat_IDs_str = cellfun(@(x) mat2str(x), Sat_IDs(2,:)', 'UniformOutput', false);
NodeRM_Sat_IDs_str = cellfun(@(x) mat2str(x), Sat_IDs(3,:)', 'UniformOutput', false);

% ✅ Build the final table for Excel export
DetailedTable = table(Time_Min, Rome_Sat_IDs_str, Milan_Sat_IDs_str, NodeRM_Sat_IDs_str, ...
    Rome_PktReception, Milan_PktReception, NodeRM_PktReception, ...
    'VariableNames', {'Time_Min','Rome_Sat_IDs','Milan_Sat_IDs','NodeRM_Sat_IDs', ...
    'Rome_PktReception','Milan_PktReception','NodeRM_PktReception'});

% ✅ Write the modified table to Excel
%writetable(DetailedTable, filename, 'Sheet', 'DetailedResults');
%fprintf('✅ Detailed results saved to: %s\n', filename);


%% Save DetailedTable to Excel with a specified path
% Specify the folder path (change this to your desired folder)
folderPath = 'D:\thesis\walker\Analysis-and-Simulation-of-LoRaWAN-LR-FHSS-main (1)\Analysis-and-Simulation-of-LoRaWAN-LR-FHSS-main';
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end
filename = fullfile(folderPath, 'DetailedResults.xlsx');
writetable(DetailedTable, filename, 'Sheet', 'DetailedResults');
fprintf('Detailed results saved to: %s\n', filename);

%% 📊 Final Results
fprintf('✅ Overall Success Rate for Rome: %.2f%%\n', mean(SuccessRate(1, :)) / Num_Packets * 100);
fprintf('✅ Overall Success Rate for Milan: %.2f%%\n', mean(SuccessRate(2, :)) / Num_Packets * 100);
fprintf('📡 Total Packets Successfully Received by NodeRM: %d\n', sum(Received_Packets_NodeRM));
toc;

%% GRAPH PLOTTING

% Create a common time axis (in minutes)
time_minutes = Visible_Sat_Matrix(:, 1);

% Graph 1: Collisions Over Time for Rome (Node 1) and Milan (Node 2)
figure;
plot(time_minutes, Collisions(1, :), 'r-o', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
plot(time_minutes, Collisions(2, :), 'b-o', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('Time (min)'); ylabel('Number of Collisions');
title('Collisions over Time for Rome and Milan');
legend('Rome (Node 1)', 'Milan (Node 2)'); grid on;

% Graph 2: Successful Transmissions Over Time for Rome and Milan
figure;
plot(time_minutes, SuccessRate(1, :), 'r-o', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
plot(time_minutes, SuccessRate(2, :), 'b-o', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('Time (min)'); ylabel('Successful Transmissions');
title('Successful Transmissions over Time for Rome and Milan');
legend('Rome (Node 1)', 'Milan (Node 2)'); grid on;

% Graph 3: NodeRM Packet Reception over Time
figure;
stem(time_minutes, Received_Packets_NodeRM, 'g', 'LineWidth', 1.5, 'Marker', 'o');
xlabel('Time (min)'); ylabel('NodeRM Reception (1 = Received)');
title('NodeRM Packet Reception over Time'); grid on;



% ✅ Check if Latitudes and Longitudes exist before plotting
if exist('Latitudes', 'var') && exist('Longitudes', 'var')
    figure;
    hold on;
    grid on;
    
    % 🌍 **Plot Ground Stations**
    scatter(Node_Coordinates(:, 2), Node_Coordinates(:, 1), 100, 'ro', 'filled');  % Red markers for nodes
    text(Node_Coordinates(:, 2) + 1, Node_Coordinates(:, 1), {'Rome', 'Milan', 'NodeRM'}, 'FontSize', 12);

    xlabel('Longitude (°)');
    ylabel('Latitude (°)');
    title('2D Animated Ground Tracks of Satellites');
    xlim([-180, 180]);
    ylim([-90, 90]);

    % 🌟 **Initialize Satellite Plot Objects**
    sat_plots = gobjects(Num_Satellites, 1);
    ground_tracks = gobjects(Num_Satellites, 1);

    for s = 1:Num_Satellites
        % **Plot empty ground track (will update over time)**
        ground_tracks(s) = plot(Longitudes(s, 1), Latitudes(s, 1), 'b--', 'LineWidth', 1); % Dashed line for ground track
        sat_plots(s) = plot(Longitudes(s, 1), Latitudes(s, 1), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b'); % Blue circles for satellites
    end

    % 🛰️ **Animate the Satellite Movement**
    for t = 1:num_steps
        for s = 1:Num_Satellites
            % Update satellite positions in the plot
            set(sat_plots(s), 'XData', Longitudes(s, t), 'YData', Latitudes(s, t));

            % **Update ground track by plotting past positions**
            set(ground_tracks(s), 'XData', Longitudes(s, 1:t), 'YData', Latitudes(s, 1:t));
        end

        % 📌 **Update Plot Title with Time**
        title(sprintf('2D Animated Ground Tracks of Satellites (Time: %.2f min)', Time_Vector(t) / 60));

        pause(0.1);  % Small pause for animation effect
    end

    hold off;
else
    fprintf('⚠️ Warning: Latitudes and Longitudes are not available for plotting.\n');
end