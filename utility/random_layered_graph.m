% Build a layered graph with some element of randomness.
% node_layers = cell in which each row is a vector of the node labels in its layer
% i_vol_dist = 'bin', 'geo', or 'unif'
% ----
% adj = adjacency matrix but with flows as edge weights

% IAP Updated Params
% node_counts = array with linear indices such that node_counts(s) = the number of nodes in stage s

% Notes to user:
% * Every <distrib>_assign_degrees function only assigns degrees *on top* of
% the guaranteed-at-least-one per node.
% * Every <distrib>_assign_degrees function only returns a number of edges that is somewhat close to num_edges (the target)
function [adj, node_assignments, init_vols, stage_ends] = random_layered_graph(node_counts,i_vol_dist,layer_params, show_plots, display_network)

    [num_layers, ~] = size(layer_params);
    num_stages = length(node_counts);
    
    % explicitly assign the nodes to stages
    node_assignments = cell(num_stages, 1); % 1 row for each stage
    node_assignments{1} = 1:node_counts(1); % initialize
    for s = 2:num_stages
        last_in_prev_stage = node_assignments{s-1}(end); 
        node_assignments{s} = (last_in_prev_stage + 1):(last_in_prev_stage + node_counts(s)); % allocate nodes to next stage
    end
    
    num_nodes = node_assignments{end}(end); % total number of nodes in network
    adj = zeros(num_nodes); % initialize adjacency matrix to size of network
    
    
    % Assign initial volumes to the first stage, based on i_vol_dist
    n1 = node_counts(1); % number of nodes in first stage
    switch i_vol_dist
        case 'bin'
            init_vols = bin_vols(n1);
        case 'geo'
            init_vols = geo_vols(n1);
        case 'unif'
            init_vols = unif_vols(n1);
        case 'determ'
            init_vols = determ_vols(n1);
    end % end switch
    
    % plot i_vol_dist
    if show_plots
        figure()
        histogram(diag(init_vols));
        title(strcat('Initial Volumes in Stage 1,', {' '}, i_vol_dist)); % title histogram
        xlabel('Proportion of Total Volume');
        ylabel('Counts');
    end % end if
    
    
    % ASSIGN EDGES WITHIN LAYERS 
    for l=1:num_layers
        prev_nodes = node_assignments{l}; % for readability
        next_nodes = node_assignments{l+1}; % vector of integer node IDs in right stage
        
        % from input params, compute num edges in each layer
        switch layer_params{l, 2}
            case 'out'
                num_edges = round(length(prev_nodes)*layer_params{l, 1});
            case 'in'
                num_edges = round(length(next_nodes)*layer_params{l, 1});
        end % end switch
        
        % ensure the layer has at least as many edges as nodes on either side
        num_edges = max([num_edges, length(prev_nodes), length(next_nodes)]);
        
        % based on user-specified distributions, assign initial degrees to nodes
        switch layer_params{l, 3} % assign degs to 'prev' stage
            case 'identical'
                prev_degrees = identical_assign_degrees(num_edges, length(prev_nodes));
            case 'bin'
                prev_degrees = binomial_assign_degrees(num_edges, length(prev_nodes), length(next_nodes));
            case 'unif'
                prev_degrees = unif_assign_degrees(num_edges, length(prev_nodes), length(next_nodes));
            case 'geo'
                prev_degrees = geo_assign_degrees(num_edges, length(prev_nodes), length(next_nodes));
        end % end 'prev' stage switch
        
        switch layer_params{l, 4} % assign degs to 'next' stage
            case 'identical'
                next_degrees = identical_assign_degrees(num_edges, length(next_nodes));
            case 'bin'
                next_degrees = binomial_assign_degrees(num_edges, length(next_nodes), length(prev_nodes));
            case 'unif'
                next_degrees = unif_assign_degrees(num_edges, length(next_nodes), length(prev_nodes));
            case 'geo'
                next_degrees = geo_assign_degrees(num_edges, length(next_nodes), length(prev_nodes));
                
        end % end 'next' stage switch

        % add edges to the adjacency matrix, guided by the generated degree distributions
        if sum(prev_degrees) <= sum(next_degrees)
            [a, prev_degs, next_degs] = stitch(adj, prev_nodes, next_nodes, prev_degrees, next_degrees);
            adj = a;
        else 
            [at, next_degs, prev_degs] = stitch(adj.', next_nodes, prev_nodes, next_degrees, prev_degrees);
            adj = at.';       
        end
        
        % ASSIGN VOLUME FLOWS TO EDGES
        first = node_assignments{l, 1}(1); % label of first node on prev side of layer
        last = node_assignments{l,1}(length(node_assignments{l})); % label of last node on next side of layer
        
        switch i_vol_dist % assign flows by same distribution as initial volumes
        case 'bin'
            for node = first:last % iterate over all nodes on left
                edges_to = find(adj(node, :));
                adj(node, adj(node,:) > 0) = bin_vols(length(edges_to));
            end % end for
        case 'geo'
            for node = first:last % iterate over all nodes on left
                edges_to = find(adj(node, :));
                adj(node, adj(node,:) > 0) = geo_vols(length(edges_to));
            end % end for
        case 'unif'
            for node = first:last % iterate over all nodes on left
                edges_to = find(adj(node, :));
                adj(node, adj(node,:) > 0) = unif_vols(length(edges_to));
            end % end for
        case 'determ'
            for node = first:last % iterate over all nodes on left
                edges_to = find(adj(node, :));
                adj(node, adj(node,:) > 0) = determ_vols(length(edges_to));
            end % end for
        end % end switch
        
        if show_plots
            % generate and display degree histograms
            figure()
            histogram(prev_degs);
            hold on
            histogram(next_degs);
            title(strcat('In/Out Degrees of Layer ', {' '}, num2str(l))); % title histogram
            legend('Out', 'In')
            xlabel('Degree')
            ylabel('Counts')
            hold off
        end % end if
    end % end for that iterates through layers
    
    % create and show network visualization
    if display_network
        bg = biograph(adj);
        nodes = bg.Nodes;
        for i = 1:length(nodes)
            if i < 10
                nodes(i).Label = [' ', int2str(i), ' ']; % for uniformity of circle display size
            else
                nodes(i).Label = int2str(i);
            end
            nodes(i).Shape = 'circle';
        end
        set(bg, 'ShowTextInNodes', 'Label');
        view(bg)       
    end
    
    % return this variable for user later
    stage_ends = zeros(1,num_stages);
    for i = 1:num_stages
        stage_ends(i) = node_assignments{i}(end);
    end

end % end function



% volume-assignment functions

% Given a number of nodes to distribute volume among, samples and returns a
% plausible diagonal matrix of fractional volumes from a binomial distribution.
function [vols] = bin_vols(num_nodes)
    Vm = 20; % default value for max volume ratio of highest magnitude to lowest magniture
    p = .5; % default value
    
    vols = ones([1, num_nodes]) + binornd(Vm - 1, p, [1, num_nodes]);
    vols = round(vols/sum(vols), 4); % normalize and round to 4 digits after decimal
end 

% Given a number of nodes to distribute volume among, samples and returns a
% plausible diagonal matrix of fractional volumes from a geometric distribution.
function [vols] = geo_vols(num_nodes)
    p = 0.5; % default value
    
    vols = ones([1, num_nodes]) + geornd(p, [1, num_nodes]);
    vols = round(vols/sum(vols), 4);
end 

% Given a number of nodes to distribute volume among, samples and returns a
% diagonal matrix of fractional volumes from a uniform distribution.
function [vols] = unif_vols(num_nodes)
    Vm = 50; % default value of largest possible ratio of max to min vol
    
    vols = ones([1, num_nodes]) + unidrnd(Vm, [1, num_nodes]);
    vols = round(vols/sum(vols), 4);
end 

% Given a number of nodes to distribute volume among, returns a diagonal
% in which each of num_nodes has identical fractional volume.
function [vols] = determ_vols(num_nodes) % assign same amount of volume to every node
    vols = round(ones([1, num_nodes])/num_nodes, 4);
end 

% degree-assignment functions

% Given a target number of edges, returns a uniform distribution of them across num_nodes nodes.
% (Will be a uniform distribution if num_nodes divides num_edges, or two different integers otherwise.)
% (For example: num_edges = 4 with num_nodes = 3 might result in a distribution [1 1 2].)
function [degrees] = identical_assign_degrees(num_edges, num_nodes)
    target_mean_deg = num_edges / num_nodes; 
    int_mean_deg = floor(target_mean_deg); % smallest int that's <= target_mean_deg
    extras = num_edges - int_mean_deg*num_nodes;
    
    degrees = repmat(int_mean_deg, 1, num_nodes); % 1xnum_nodes matrx filled with int_avg_deg
    
    % randomly select indices of nodes that'll receive the extra edges (no replacement, bc no node should get more than one extra edge)
    extra_receivers = datasample(1:num_nodes, extras, 'Replace' , false);
    for i = 1:length(extra_receivers)
        degrees(extra_receivers(i)) = int_mean_deg + 1;
    end % end for
end % end function

% Given a target number of edges, a number of nodes to distribute them
% among, and the maximum degree allowed for any node, samples and returns a
% vector describing an approximately-binomial degree distribution.
function [degrees] = binomial_assign_degrees(num_edges, num_nodes, max_deg)
    target_mean_deg = num_edges / num_nodes;
    n = max_deg - 1; % n in binom distrib, shifted by one b/c every node already has one edge
    p = (target_mean_deg - 1)/(n); % p in binom distrib. Derived from E[x] = np
    
    degrees = ones(1, num_nodes) + binornd(n, p, [1, num_nodes]); % ensure deg >=1 
end % end function

% Given a target number of edges and a number of nodes to distribute them
% among,  samples and returns a vector describing a uniform degree distribution.
function [degrees] = unif_assign_degrees(num_edges, num_nodes, max_deg)
    target_mean_deg = num_edges / num_nodes;
    max_deg = min(ceil(target_mean_deg)*2-1, max_deg); % ensure we assign roughly num_edges, but no more than max_deg per node
    
    degrees = unidrnd(max_deg, [1, num_nodes]); 
end

% Given a target number of edges, a number of nodes to distribute them
% among, and the maximum degree allowed for any node, samples and returns a
% vector describing a truncated geometric degree distribution (as described in 'Network Generator Variables V5')
function [degrees] = geo_assign_degrees(num_edges, num_nodes, max_deg) 
    p0 = 1 / (num_edges/num_nodes); % 1/(the expected degree)
    p0 = min([1 p0]); % do not allow to be > 1
    p0 = max([1/max_deg p0]); % do not allow to generate an expected value greater than max deg
    
    C = 0;
    for m = 1:max_deg
        C = C + p0*(1-p0)^(m-1); 
    end % end for
    
    degrees = rand([1, num_nodes]);
    degrees = ceil(log(1-C*degrees)/log(1-p0)); % ensure that every node has deg >=1
end

% functions that assemble the network from parameters

% Given the adjacency matrix so far, as well as lists of nodes in the previous and
% next stages (with corresponding target degree distributions), returns the
% original adjacency matrix updated to contain edges that (as much as possible) satisfy the newly
% input degree distributions. (Also returns the tweaked degree distributions represented in adj.) 
function [adj, prev_degs, next_degs] = stitch(adjacency, prev_stage, next_stage, prev_degs, next_degs)
    adj = adjacency;
    
    % make sure right and left side have same total degree by randomly adding edges to smaller side
    if sum(next_degs) ~= sum(prev_degs) % only fix them if they're unequal
        if sum(next_degs) < sum(prev_degs) % if the next stage has lower degree
            next_degs = add_to_equalize(next_degs, sum(prev_degs), length(prev_degs));
        else % else do same thing but in other direction
            prev_degs = add_to_equalize(prev_degs, sum(next_degs), length(next_degs));
        end % end if
    end % end if

    % connect nodes according to generated degree distribution
    remaining_prev_nodes = [prev_stage; prev_degs]; % node ID in first row, corresponding degree in 2nd row
    remaining_next_nodes = [next_stage; next_degs];

    changed_degree_in_last_iteration = true; % double-checked convergence condition
    for node = remaining_prev_nodes(1, :) % iterate through prev_nodes
        if changed_degree_in_last_iteration % make sure something happened in last round
            node_deg = remaining_prev_nodes(2, remaining_prev_nodes(1,:) == node);
            [~, perm_to_sort_cols] = sort(remaining_next_nodes(2,:), 'descend'); 
            remaining_next_nodes = remaining_next_nodes(:, perm_to_sort_cols); % sort nodes on other side by descending degree
            connect_to_nodes = remaining_next_nodes(1, 1:node_deg); % IDs of nodes to attach to
            
            for next_node = connect_to_nodes % iterate through elements in linkto
                col_ind_of_next_node = find(remaining_next_nodes(1,:) == next_node, 1);
                deg_of_next_node = remaining_next_nodes(2, col_ind_of_next_node);
                if deg_of_next_node > 0 % shouldn't happen, but just to be safe
                    adj(node, next_node) = 1; % just one way, since it's a directed graph
                    remaining_next_nodes(2, col_ind_of_next_node) = deg_of_next_node - 1; % decrement degree 
                    changed_degree_in_last_iteration = true;
                end % end if
            end % end for  
        end % end if
    end % end for
end % end stitch

% Given a vector v of integers and an integer new_mag that must be greater
% than the current magnitude of v, adds 1 incrementally to random positions
% of v until it has magnitude new_mag, subject to the constraint that no entry of
% v can exceed max_entry. Returns the resulting v.
function [v] = add_to_equalize(v, new_mag, max_entry)
    assert(sum(v) <= new_mag, 'new_mag must be greater than the current magnitude of v.');
    assert(new_mag <= max_entry * length(v), 'impossible to reach desired new_mag of v when constrained by max_entry.');
    
    diff = new_mag - sum(v); % number of edges remaining to add
    strcat('Adding', {' '}, num2str(diff), ' edges to smaller-degree side. (It originally had', {' '}, num2str(sum(v)), '.)')
    while diff ~= 0
        i = floor(rand()*length(v)) + 1; % random integer between 1 and length(v), inclusive
        % linear search until reach valid index to add to
        while v(i) >= max_entry
            i = mod(i + 1, length(v)) + 1;
        end
        v(i) = v(i) + 1;
        diff = diff - 1;  
    end % end while
end
