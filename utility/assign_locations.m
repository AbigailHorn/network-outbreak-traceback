% January 2017 implementation of location assignment

% --- INPUTS ---
% node_assignments = cell in which each row is a vector of the node ids in its stage (as output from random_layered_graph)
% adj = adjacency matrix of the network, with transition probabilities as edge weights.
% show_plots = a boolean. If true, will display informative pop-up graphs.

% --- OUTPUTS ---
% dist_mat = adjacency matrix with distances as edge weights.
% locations = cell in which each row is a column vector of the (x,y) tuples.

function [dist_mat, locations] = assign_locations(node_assignments, adj, show_plots) 
    w = 2500; % grid width
    h = 1500; % grid height
    node_data = []; % purely an internal variable: node_data(n, :) = [stage of n, n's index within stage, n's x-coord, n's y-coord]
    
    % maps node to stage it's in, as well as its index within the stage (for housekeeping)
    % meaning: node_data(n, 1) = the stage of node n; node_data(n, 2) = index of node n within its stage
    end_of_prev_stage = 0; % temporary variable
    for s = 1:length(node_assignments)
        for n = node_assignments{s}
            node_data(n,1) = s; % map stage_of_node(n,1) to its stage
            node_data(n,2) = n - end_of_prev_stage; % map stage_of_node(n,2) to its index within stage
        end
        end_of_prev_stage = length(node_data);
    end
    
    [num_stages, ~] = size(node_assignments);
    locations = cell(num_stages, 1); % initialize 1 row for each node layer
    
    % uniformly distribute locations of start and end stages
    locations{1, 1} = unif_locs(w, h, length(node_assignments{1}));
    
    locations{num_stages, 1} = unif_locs(w, h, length(node_assignments{num_stages}));
    
    % handle networks with more than two stages
    switch num_stages
        case 3
            locations{2, 1} = grid_midpoint_locs(adj, node_assignments{2}, 'in', node_data, locations);
        case 4
            locations{2, 1} = grid_midpoint_locs(adj, node_assignments{2}, 'in', node_data, locations);
            locations{3, 1} = grid_midpoint_locs(adj, node_assignments{3}, 'out', node_data, locations);
        case 5
            locations{2, 1} = grid_midpoint_locs(adj, node_assignments{2}, 'in', node_data, locations);
            locations{4, 1} = grid_midpoint_locs(adj, node_assignments{4}, 'out', node_data, locations);
            locations{3, 1} = grid_midpoint_locs(adj, node_assignments{3}, 'out', node_data, locations);
        otherwise
            if num_stages > 5
                locations{2, 1} = grid_midpoint_locs(adj, node_assignments{2}, 'in', node_data, locations);
                locations{num_stages-1, 1} = grid_midpoint_locs(adj, node_assignments{num_stages-1}, 'out', node_data, locations);
                locations{num_stages-2, 1} = grid_midpoint_locs(adj, node_assignments{num_stages-2}, 'out', node_data, locations);
                for s = 2:(num_stages-3)
                    locations{s, 1} = grid_midpoint_locs(adj, node_assignments{s}, 'in', node_data, locations);
                end % end for
            end % end if
    end % end switch
    
    % create "adjacency matrix" with distances as edge weights
    dist_mat = adj;
    for r = 1:length(dist_mat) % for every row
        for c = find(dist_mat(r, :)) % for every column in the row with a nonzero element
            % find start and end loc of the corresponding edge
            start_loc = locations{node_data(r,1),1}(node_data(r,2),:); 
            end_loc = locations{node_data(c,1),1}(node_data(c,2),:);
            dist_mat(r,c) = sqrt(sum((start_loc - end_loc).^2)); % cartesian distance
        end % end for
    end % end for
    
    
    % maintain node_data to include location x-coord in 3rd col, y-coord in 4th
    for s = 1:length(locations)
        start_n = node_assignments{s}(1); % first node index in current stage
        end_n = node_assignments{s}(length(node_assignments{s})); % last node index in current stage
        node_data(start_n:end_n, 3) = locations{s, 1}(:, 1);
        node_data(start_n:end_n, 4) = locations{s, 1}(:, 2);
    end
    
    if show_plots
        figure();
        scatter_plot_network(node_data, adj);
    end % end if
    
end % end function
   
% h,w = height and width of potential location territory
% num = number of locations to assign within it
% --
% locations = a num-row, two-column matrix in which each row represents the
% (x, y) coordinates of a node. (loications are distributed uniformly at random)
function [locations] = unif_locs(w, h, num)
    locations = [unidrnd(w, 1, num).' unidrnd(h, 1, num).'];
end

% stage = list of all node labels in current stage
% adj = the network's adjacency matrix
% direxn = either 'in' or 'out' depending on whether location assignment
% is based on either incoming or outgoing edges
% defining_stage = index of stage that defines distances of this stage (usually prev or next)
% locs = current list of locations
% --
% locations = a num-row, two-column matrix in which each row represents the
% (x, y) coordinates of a node. (locations are distributed uniformly at random)
function [locations] = grid_midpoint_locs(adj, stage, direxn, stage_of_node, locs)
    locations = []; % initialize
    current_stage = stage_of_node(stage(1), 1);
    % if we want to assign based on incoming nodes, transform adj matrx so
    % assignment algorithm can proceed in same manner as for outgoing
    if strcmp(direxn, 'in')
        adj = adj.';
        defining_stage = current_stage - 1;
    else
        defining_stage = current_stage + 1;
    end
    
    % for each node in the stage
    for n = stage
        out_locs = zeros(0,2); % initialize to empty 2-col vector
        % assign it outgoing grid midpoint distance
        for out = find(adj(n, :)) % for each outgoing edge end node
            out_locs = [out_locs; locs{defining_stage}(stage_of_node(out,2), :)];
        end
        
        % set avg_loc to correct grid midpoint location
        [r, ~] = size(out_locs);
        if r == 1
            avg_loc = out_locs;
        else
            avg_loc = sum(out_locs)/length(out_locs);
        end
        locations = [locations; avg_loc];
    end % end for
end % end grid_midpoing_locs

% n_d = node_data array as upkept in main function
% displays a scatter plot of locations of all nodes, with each stage a different color.
function [] = scatter_plot_nodes(n_d) 
   % x, y, groups, colrs, symbol, size
   gscatter(n_d(:, 3), n_d(:,4), n_d(:,1), [], [], 40);
end % end scatter_plot_nodes

function [] = scatter_plot_network(n_d, adj)
   % heuristic formatting things
   headLength = 8;
   headWidth = 5;
   scaleFactor = .96;
   % edge_out is vector of nodes edge leaves from
   % edge_out is vector of nodes edge goes into
   edge_ends = find(adj);
   [edge_out, edge_in] = ind2sub(size(adj), edge_ends);
   [ne, ~] = size(edge_out); % number of edges in network
   for e = 1:ne % for each edge in the graph
       % '-ok' = solid black line. Look into making directed
       xOut = n_d(edge_out(e),3);
       yOut = n_d(edge_out(e), 4);
       xIn = n_d(edge_in(e),3);
       yIn = n_d(edge_in(e), 4);
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
   
   % x, y, groups, colrs, symbol, size
   gscatter(n_d(:, 3), n_d(:,4), n_d(:,1), [], [], 40);
   %hold on
   hold off
   
   title('Network Geography ','FontSize',16);
    
end % end scatter_plot_networks