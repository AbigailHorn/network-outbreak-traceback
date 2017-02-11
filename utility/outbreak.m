% outbreak6: gives times consumed, which can be used to calculate number of illnesses averted
% Abby's outbreak with f=.0077 but 1/500 reporting (vs. 1/50)
% dispersion = 'max' or 'min'
% distances = adjacency matrix with distances as edge weights
% locations = cell in which each row is a vector of the (x,y) tuples (as output from assign_distances)
% flows = adjacency matrix with outgoing flows as edge weights
% init_vols = square diagonal matrix of initial volumes
% contam_retailers = 2 x num_contam_retailers matrix. first row is the contam node IDs, 2nd is the reporting times at them
% transport_dev_frac = variance of transportation time will be set to (.25*cum_time)^2 (transport_dev_frac specifically is the fraction * cum_time)
% storage_is_determ = binary switch for stochastisity of all post-transport variables, including shelf life, etc.
function [contam_farm, contam_reports] = outbreak(dispersion, src_farm, distances, locations, flows, init_vols, num_stages, return_all_reports, show_plots, transport_dev_frac, storage_is_determ)
    'Propagating outbreak on input network...'
    % calculate volumes and acreage at each farm node
    farm_vol_props = init_vols.'; % column vector of proportions of initial volumes at each farm
    num_farms = length(farm_vol_props);
    farm_vols = 35737*num_farms*farm_vol_props; % daily volume per farm (in heads of lettuce)
    farm_acrg = 228*num_farms*farm_vol_props; % farm acerage
    
    % generate contamination event
    f = .01; % fraction of farm's producing land that's contaminated
    if src_farm <= 0
        contam_farm = randsample(num_farms, 1, true, farm_vol_props); % weight selection of contaminated farm by intial volume fractions
    else
        contam_farm = src_farm;
    end 
    
    contam_days = ceil(f*130); % number of days the outbreak continues
    if f < 1/130
        ppd = ceil(f*farm_acrg*20400/960);  % column vector of contaminated pallets per day
        lpd = ceil(f*farm_acrg*20400/23040); % loads per day
    else
        ppd = ceil(farm_vols/960); % vector of each farm's pallets per day
        lpd = ceil(farm_vols/23040);
    end
    
    contam_paths = cell([contam_days, 1]); % column cell; each entry will store paths of all contaminated pallets loosed that day (as array)
    
    % Based on user-input dispersion ('min' or 'max'), sample the appropriate number of distinct contamination paths. 
    % Then generate contam_times_till_retail (and contam_dists_till_retail) by palette.
    switch dispersion
        case 'min' % outbreak propagates by load
            for day = 1:contam_days % for each day in the outbreak
                % for each contaminated load, sample a path. 
                contam_paths{day} = sample_paths(contam_farm, lpd(contam_farm), flows, num_stages); % each col is a path; lpd cols total
            end % end for
                
            % convert distinct paths to times for each pallet
            contam_times_till_retail = zeros(contam_days*lpd(contam_farm)*24, 2); % col1: time contamination reaches retailer; col2: ID of retailer reached
            contam_dists_till_retail = zeros(contam_days*ppd(contam_farm), 2); % col1: dist contamination travels to retailer; col2: ID of retailer reached
            l_counter = 1; % counts which pallet we're currently tracking
            for day = 1:length(contam_paths) % for each day 
                for path_ind = 1:lpd(contam_farm) % for all paths in that day
                    path = contam_paths{day}(:, path_ind);
                    contam_ret = path(num_stages);
                    contam_times_till_retail((24*(l_counter-1)+1):(24*l_counter), :) = repmat([time(path, distances, transport_dev_frac) contam_ret], 24, 1); % 24 pallets per load
                    contam_dists_till_retail((24*(l_counter-1)+1):(24*l_counter), :) = repmat([dist(path, distances) contam_ret], 24, 1); % 24 palettes per load
                    l_counter = l_counter + 1;
                end % end for
            end % end for
            
        otherwise % default: max dispersion, outbreak propagates by pallet
            for day = 1:contam_days % for each day in the outbreak
                % for each contaminated pallet, sample a path.
                contam_paths{day} = sample_paths(contam_farm, ppd(contam_farm), flows, num_stages); % each col is a path; ppd cols total
            end % end for
                
            % convert each pallet's path to time and distance
            contam_times_till_retail = zeros(contam_days*ppd(contam_farm), 2); % col1: time contamination reaches retailer; col2: ID of retailer reached
            contam_dists_till_retail = zeros(contam_days*ppd(contam_farm), 2); % col1: dist contamination travels to retailer; col2: ID of retailer reached
            p_counter = 1; % counts which pallet we're currently tracking
            for day = 1:length(contam_paths) % for each day 
                for path_ind = 1:ppd(contam_farm) % for all paths in that day
                    path = contam_paths{day}(:, path_ind);
                    contam_ret = path(num_stages);
                    contam_times_till_retail(p_counter,:) = [time(path, distances, transport_dev_frac) contam_ret];
                    contam_dists_till_retail(p_counter,:) = [dist(path, distances) contam_ret];
                    p_counter = p_counter + 1;
                end % end for
            end % end for   
    end % end switch
    
    % histogram showing contaminant times till retail, BY PALLET
    if show_plots
        figure();
        % histogram with one bin per day of outbreak
        firstDay = floor(min(contam_times_till_retail(:,1)));
        lastDay = ceil(max(contam_times_till_retail(:,1)));
        histogram(contam_times_till_retail(:,1),  lastDay-firstDay + 1);
        title(strcat('Contaminated Pallet Times Till Retail from Origin at Farm', {' '}, int2str(contam_farm))); 
        xlabel('Time (days)'); 
        ylabel('Counts');
    
        % scatterplot of time contamination reaches retailer vs. dist from location of origin
        figure()
        scatter(contam_times_till_retail(:,1), contam_dists_till_retail(:,1))
        title('Contamination Times Till Retail vs. Distance Traveled'); 
        xlabel('Time (days)'); 
        ylabel('Distance (miles)');
    end % end if show_plots
    
    % disaggregate pallets down to heads 
    times_by_head = zeros(960*length(contam_times_till_retail), 2); % preallocate
    for pal = 1:length(contam_times_till_retail)
        times_by_head(((pal-1)*960+1):(pal*960), :) = repmat(contam_times_till_retail(pal, :), 960, 1);
    end
    
    strcat(int2str(length(times_by_head)), {' '}, 'heads of lettuce contaminated.')
    
    % add time spent at retail + time spent inside house
    times_till_consumed = []; % initialize empty column vector
    for i = 1:length(times_by_head)
        % add storage/incubation stochasticity if desired
        if storage_is_determ
            store_to_plate = 1.5 + 3; 
        else
            store_to_plate = exprnd(1.5) + wblrnd(3, 1);
        end
        % discard all heads past their shelf lives
        if store_to_plate <= 25 % if total time from farm to table is less than the shelf life 
            times_till_consumed = [times_till_consumed; (times_by_head(i, 1) + store_to_plate) times_by_head(i, 2)]; % 2nd col remains retailer ID     
        end % end if
    end % end for
    
    tic
    % deal with the heads that *are* consumed! construct epidemic curve, keep track of useful data
    reporting_times = [];
    consuming_times = []; 
    for i = 1:length(times_till_consumed)
        for ppl = 1:unidrnd(4) % for each person that consumes the contaminant
            if unidrnd(10) == 1 % if the person develops symptoms
                if unidrnd(100) == 1 % if the person reports their symptoms
                    % add storage/incubation stochasticity if desired
                    if storage_is_determ
                        reporting_times = [reporting_times; (times_till_consumed(i,1) + 4) times_till_consumed(i, 2)]; % deterministic
                    else
                        reporting_times = [reporting_times; (times_till_consumed(i,1) + (normrnd(4,1))) times_till_consumed(i, 2)]; % time of report = consumption time + incubation time
                    end
                    consuming_times = [consuming_times; (times_till_consumed(i,1)) times_till_consumed(i, 2)]; % Edit for outbreak6 : time of consumption for those reporting
                end % end if
            end % end if
        end % end for 
    end % end for
    strcat('Spent', {' '}, num2str(toc), {' '}, 'seconds calculating reporting times.')
    
    if show_plots
        figure();
        firstDay = floor(min(reporting_times(:,1)));
        lastDay = ceil(max(reporting_times(:,1)));
        histogram(reporting_times(:,1), lastDay-firstDay);
        title('Times of illness reports'); 
        xlabel('Time (days)'); 
        ylabel('Counts');
    
        % plot geographic spread of outbreak
        visualize_spread(locations, contam_farm, contam_paths, flows, dispersion);
    end 
    
    % arrange output information into standard format
    if return_all_reports % return every report, along with reporting and consuming times
        reporting_times = reporting_times.';
        consuming_times = consuming_times.'; 
        contam_reports = zeros(size(reporting_times));  
        contam_reports(1, :) = reporting_times(2, :); % first row = nodes contaminated
        contam_reports(2, :) = reporting_times(1, :); % second row = times contam reported at
        contam_reports(3, :) = consuming_times(1, :); % third row = times contaminated food eaten  
    
    else % else just return each contam_node along with the time of its first contamination report
        contam_reports = [];
        [num_stages, ~] = size(contam_paths{day});
        for day = 1:length(contam_paths) % for each day in outbreak
            contam_reports = union(contam_reports, contam_paths{day}(num_stages, :));
        end % end for
        
        mapped_indices = zeros(1, max(contam_reports)); % will map node_ID to its index in contam_retailers
        for i = 1:length(contam_reports)
            mapped_indices(contam_reports(i)) = i; % so, mapped_indices(node_ID) = the index in contam_retailers that describes node_ID
        end % end for
        
        contam_reports = contam_reports.'; % transpose
        contam_reports(2, 1:length(contam_reports)) = Inf*ones(1, length(contam_reports)); % initialize to find first incidence times
        for report = 1:length(reporting_times) % find first times of incidence
            % if the current report happened earlier than the earliest logged one at same stage, update
            if reporting_times(report, 1) < contam_reports(2, mapped_indices(reporting_times(report, 2)))
                contam_reports(2, mapped_indices(reporting_times(report, 2))) = reporting_times(report, 1);
            end % end if
        end % end for
    end % end if
end % end outbreak

% display plot of network on map, such that
% red = contaminated retailer nodes
% purple = noncontaminated retailer nodes
% blue = origin farm node
% black = all other nodes
% ---
% dists = adjacency matrix with distances as input

% locations = cell with num_nodes_in_stage_s-by-2 matrix in entry for stage {s},
% containing x and y coordinate pairs for nodes, by ID order
function[] = visualize_spread(locations, start, contam_paths, flows, dispersion)
    figure(); % run to make sure opens on next figure
    total_num_retailers = 0; % counter for later
    
    % recall: contam_paths is matrix in which each column is a contaminant path
    contam_retailers = [];
    [stages, ~] = size(contam_paths{1});
    for day = 1:length(contam_paths)
        contam_retailers = [contam_retailers contam_paths{day}(stages, :)]; % nonunique entries (this changes on line 235)
    end
    
    viz_data = []; % col 1 = x val, col 2 = y val, col 3 = color
    
    % 1 = contam retail nodes
    % 2 = noncontam retail nodes
    % 3 = start farm node
    % 4 = all other nodes
    node_id = 1;
    for s = 1:length(locations) % for each stage
        for loc = locations{s}.' % iterate through each (x,y) pair in the stage
            if s == length(locations) % if in retail stage
                viz_data(node_id, 3) = 2; % initially mark retail nodes as uncontam retail nodes
                total_num_retailers = total_num_retailers + 1;
            else
                viz_data(node_id, 3) = 4; % else mark as 'other node'
            end % end if
            viz_data(node_id,1) = loc(1); % map stage_of_node(n,1) to its x
            viz_data(node_id,2) = loc(2); % map stage_of_node(n,2) to its y
        node_id = node_id + 1;
        end% end for
    end % end for
    
    % set color of origin farm
    viz_data(start, 3) = 3;
    
    for contam_ret = unique(contam_retailers)
        viz_data(contam_ret, 3) = 1;
    end % end for
    
    % scatter plot the network according to the coloring scheme
    % this part is an almost exact duplicate of scatter_plot_network in assign_dists. potentially clean up to address redundancy?
    
    % edge_out is vector of nodes edge leaves from
    % edge_out is vector of nodes edge goes into
    % heuristic formatting things
   headLength = 8;
   headWidth = 5;
   scaleFactor = .96;
   % edge_out is vector of nodes edge leaves from
   % edge_out is vector of nodes edge goes into
   edge_ends = find(flows);
   [edge_out, edge_in] = ind2sub(size(flows), edge_ends);
   [ne, ~] = size(edge_out); % number of edges in network
   for e = 1:ne % for each edge in the graph
       % '-ok' = solid black line. Look into making directed
       xOut = viz_data(edge_out(e),1);
       yOut = viz_data(edge_out(e), 2);
       xIn = viz_data(edge_in(e),1);
       yIn = viz_data(edge_in(e), 2);
       xComp = xIn-xOut;
       yComp = yIn-yOut;
       %xComponent = n_d(edge_in(e),3) - n_d(edge_out(e),3);
       %yComponent = n_d(edge_in(e), 4) - n_d(edge_out(e), 4);
       %plot([n_d(edge_out(e),3), n_d(edge_in(e),3)], [n_d(edge_out(e), 4), n_d(edge_in(e), 4)], '-ok')
       % hq = quiver(n_d(edge_out(e),3), n_d(edge_out(e),4), xComponent*(1+extendBy), yComponent*(1+extendBy));
            ah = annotation('arrow',...
                'headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth, 'Color', [.35, .35, .55]);
            set(ah,'parent',gca);
            set(ah,'position',[xOut, yOut, xComp, yComp*scaleFactor]);
       hold on
   end % end for
   
   % contam retail, uncontam retail, origin farm, else
   if length(unique(viz_data(:,3))) < 4
       colorSpec = [1 .184 .012; .235 .906 1; 0.0 0.0 0.0]; % if no uncontaminated retailers left, need to shift colorspec
       colorSpec = [1 .6 .6; 1 .184 .012; 0.3 0.3 0.3];
   else
       colorSpec = [1 .6 .6; .235 .906 1; 1 .184 .012; 0.3 0.3 0.3];
   end
   
    % x, y, groups, colrs, symbol, size
    gscatter(viz_data(:, 1), viz_data(:,2), viz_data(:,3), colorSpec, [], 40, 'off');
    % hold 
    titleStr = strcat('Geography of Contamination Spread', {' ('},dispersion, ' Dispersion)');
   contamSpreadInfoStr = strcat(int2str(length(unique(contam_retailers))), {' out of '}, int2str(total_num_retailers), {' retailers contaminated'});
   title(titleStr,'FontSize',16);
   xlabel(contamSpreadInfoStr, 'FontSize', 14);
   hold off
    
end


% Given the column vector of a path, the network's distance-weighted adjacency matrix, and 
% the fraction of a path's mean travel-time spanned by one standard deviation of the travel time distribution,
% samples (from a binomial distribution) a possible transit time for the path, in days. Returns this time.
% (Setting transport_dev_frac = 0 results in deterministic travel times.)
function [time] = time(path, dists, transport_dev_frac)
    time = 0; 
    for i = 2:length(path) % sum up mean travel time for each path leg
        determ = dists(path(i-1), path(i))/630; % 45mph * 14 driving hours / day  %200; % 45?? % deterministic time component
        time = time + determ; 
    end % end for 
    time = normrnd(time, transport_dev_frac*time); % used to have default transport_dev_frac = 0.25
end % end time


% Given the column vector of a path and the network's distance-weighted adjacency matrix, 
% computes and returns the total distance along the path.
function [cum_dist] = dist(path, dists)
    cum_dist = 0;
    for i = 2:length(path) % for each edge in the path
        cum_dist = cum_dist + dists(path(i-1), path(i));
    end
end % end dists


% Given a starting node ID, a number of paths to sample, an adjacency matrix (with edges weighted by
% transition probability), and the total number of nodes in the path, samples num_paths paths along
% the graph (beginning at the start node) and returns them as a path_length-by-num_paths matrix.
function [paths] = sample_paths(start_node, num_paths, flows, path_length)
    paths = zeros([path_length, num_paths]); % each column is a distinct path
    
    for p = 1:num_paths % sample num_paths paths
        paths(:, p) = sample_path(start_node, flows, path_length); % each column is a path
    end % end for
end % end function


% Given a starting node ID, an adjacency matrix (with edges weighted by
% transition probability), and the total number of nodes in the path, samples a
% path along the graph (beginning at the start node) and returns a column vector of the node IDs along the path. 
function [path] = sample_path(start, flows, path_length)
    path = zeros([path_length, 1]); % initialize
    path(1) = start;
    
    for s = 1:(path_length-1) % sample all nodes after the start node
        path(s+1) = randsample(length(flows), 1, true, flows(path(s), :));
    end % end for 
end % end sample_path