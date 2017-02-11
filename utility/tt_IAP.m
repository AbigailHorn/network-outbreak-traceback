% Notes:
% ts_info weighted by path probability before ts optimization

% --- PRECONDITION ---
% flows = an adjacency matrix for the network with transition probabilities as edge weights
% stage_starts_ends(stage) = last node in stage
% prior_pmf = row vect of prior probability. prior_pmf(s) = probability of source node s if contam_observed_nodes is empty
% contam_observed_nodes also contains time in 2nd row
% distances = adjacency matrix with distances as edge weights
% tau = contamination time uncertainty
% P = percentage of contaminated nodes that must be connected to a source
% before it's deemed "feasible"

% --- POSTCONDITION ---
% pmf = the likelihood of each node being the source


% The new one! January 2017 time traceback formulation
% ~INPUTS~ 
% time_estimator = 'exact', 'bfs', or 'MaxP'
% use_volumes = a boolean indicting whether or not the returned estimates will be influenced by network volume information
% contam_reports = first row is the node ID report observed at, second row is the time
% heuristic = 'exact', 'bfs', or 'mean'
% include_volumes = true or false. pairs with heuristic time_estimator components
% ~OUTPUTS~
% pmf = a vector such that pmf(node_ID) the likelihood of node_ID being the source
% t_hat_star = the MLE outbreak start time
function [pmf, t_hat_star] = tt_IAP(time_estimator, use_volumes, flows, stage_ends, prior_pmf, contam_reports, distances, P)

pmf = zeros([1, stage_ends(1)]); % pmf will eventually be row vect of a prob for each potential source
num_stages = length(stage_ends); % keep generalizable to number of stages other than 4
flows = sparse(flows); % for runtime


% identify feasible sources as those which reach a fraction > P of all contaminated nodes
feasible_sources = 1:stage_ends(1); % initialize to list of all farm nodes
reaches = zeros(size(feasible_sources)); % contains num of contaminated nodes reached by farm of same index
for c_node = unique(contam_reports(1,:)) % compute all potential sources
    last_parents = c_node;
    while ~isempty(parents(last_parents, flows))
        last_parents = parents(last_parents, flows);
    end % end while
    reaches(last_parents) = reaches(last_parents) + 1; % update reach counts
end % end for
feasible_sources = find(reaches >= P*length(unique(contam_reports(1,:)))); % feasible sources IDs are those which reach a fraction-greater-than-P of contaminated nodes


% Necessary setup for the heuristic we're using
characteristic_contam_reports = unique(contam_reports(1, :));
characteristic_contam_reports(2, :) = ones(size(characteristic_contam_reports)) * Inf;
% If 'bfs', second row of characteristic_contam_times will contain time of first report
% If 'mean', second row of characteristic_contam_times will contain mean report time
switch time_estimator
    case 'exact_estimator'
        
    case 'bfs'
        % use first arrival time at each node
        for time = contam_reports % for every (retailer_ID, time) report tuple
            col_of_node = characteristic_contam_reports(1, :) == time(1); % logical array
            if time(2) < characteristic_contam_reports(2, col_of_node)
                characteristic_contam_reports(2, col_of_node) = time(2);
            end
        end % end for 
        
    case 'MaxP'
        for col = 1:length(characteristic_contam_reports(1,:))
            node = characteristic_contam_reports(1, col);
            reports_of_node = contam_reports(:, contam_reports(1, :) == node);
            characteristic_contam_reports(2, col) = mean(reports_of_node(2, :)); % set 'characteristic time' as mean reported time
        end % end for
        
end % end switch

characteristic_contam_reports

% start keeping edge_data, such that edge_data(k) = [from to length]
edge_data = zeros([nnz(flows), 3]); % one row for each edge, of [from to length]
edge_lin_inds = find(flows);
for ind_ind = 1:length(edge_lin_inds)
    [from, to] = ind2sub(size(flows), edge_lin_inds(ind_ind));
    edge_data(ind_ind, :) = [from to distances(edge_lin_inds(ind_ind))]; 
end % end for

edge_data

% Calulate volume-based path likelihood term, based on [Path Likelihood] section
I_t = eye(stage_ends(num_stages-1));
Q = flows(1:stage_ends(end-1), 1:stage_ends(end-1));
A = inv(I_t - Q)*flows(1:stage_ends(end-1), stage_ends(end-1)+1:stage_ends(end));
A

% Put it all together and make a pmf
for s = feasible_sources
    
    prob_s = 1; % Will be the probability of s being source. Multiply all prior, time, volume components onto this.
    
    % spatial/volume-based term
    path_likelihood_of_observations = 1; % multiply all path probabilities together onto this variable
    for node_ID = contam_reports(1, :)
        path_likelihood = A(s, node_ID-stage_ends(end-1)); % likelihood that the current observation started at s, summed over all possible paths
        path_likelihood_of_observations = path_likelihood_of_observations * path_likelihood; % probability of source s generating the whole collection of reports
    end % end for

    % ok do the time term now
    switch time_estimator
        case 'bfs'
            c_s = c_max_prob(s, characteristic_contam_reports(1,:), flows, edge_data); % c matrix representing MaxP tree
        case 'MaxP'
            c_s = c_shortest(s, characteristic_contam_reports(1,:), flows, edge_data); % c matrix representing shortest path tree
    end % end switch
    mu_theta = edge_data(:, 3) / 630; % mean travel time of each edge, in days
    Sigma_theta = diag(.1 * mu_theta); % CHANGE THIS COEFF
    t_vec = characteristic_contam_reports(2, :).'
    mu_d = ones(size(t_vec)) * 8.5; % mean delay from 2 rounds of storage + incubation times
    
    mu_s_star = c_s * mu_theta + mu_d
    Sigma_s = c_s * Sigma_theta * c_s.' + eye(length(characteristic_contam_reports(1,:))) * 12.25;
    
    t_hat_star = (t_vec - mu_s_star).' * inv(Sigma_s) * ones(size(mu_d)); % eq 12 numerator
    t_hat_star = t_hat_star / (ones(size(mu_d.')) * inv(Sigma_s) * ones(size(mu_d))) % divide by eq 2 denominator
    
    % now use t_hat_star for other stuff
    num_obs = length(t_vec); % retroactively define italic O
    size(ones([1, num_obs])*t_hat_star)
    size(c_s*mu_theta)
    size(mu_d)
    mu_s_hat = ones([num_obs, 1])*t_hat_star + c_s*mu_theta + mu_d
    
    max_temporal_likelihood = exp(-0.5*(t_vec - mu_s_hat).' * inv(Sigma_s) * (t_vec - mu_s_hat));
    max_temporal_likelihood = 1/sqrt(norm(Sigma_s)) * max_temporal_likelihood % add the coeff
    
    
    
    pmf(s) = path_likelihood_of_observations * prior_pmf(s) * max_temporal_likelihood % add time term onto this, eventually
end % end for over all feasible sources
pmf
pmf = pmf / sum(pmf)
feasible_sources
end % end function

% directly return the relevant c matrix data structure
% this one is based on distance! 
function c_s = c_shortest(s, contam_nodes, dists, edge_data)
    % prepare for shortest-path search over distances
    dists_bg = biograph(dists);
    
    paths = [];
    % compute shortest cascade
    for node_ID = contam_nodes % create matrix data structure. first col is probability of path; remaining cols are nodes in path
        [~, path, ~] = shortestpath(dists_bg, s, node_ID, 'Method', 'Acyclic');
        paths = [paths; path];
    end % end for 
    % paths now stores the path that would have been the shortest to each observation

    c_s = c_matrix(paths, edge_data);  
end

% directly return the relevant c matrix data structure
% this is for the 'mean' heuristic
% edge_data(s) = [from to dist] (note that edge_data(:, 3)) = theta
function c_s = c_max_prob(s, contam_nodes, flows, edge_data)
    paths = [];
    % compute max cascade
    for node_ID = contam_nodes % create matrix data structure. first col is probability of path; remaining cols are nodes in path
        [~, path] = max_prob_path(s, node_ID, flows) % each row is a different path     
        paths = [paths; path];
    end % end for 
  
    % convert max_prob_diff_traj to format of c (an O-by-K matrix)
    % 1 row for each observation, and cols are binary: 1 if in shortest path from s to the obs, else 0
    c_s = c_matrix(paths, edge_data);
end

% Given a list of paths (as lists of nodes) (one path per row) representing a diffusion trajectory,
% and edge_data as maintained internally,
% Returns the num_paths-by-K c_s matrix describing the diffusion trajectory
function c_s = c_matrix(paths, edge_data)
    [num_paths, ~] = size(paths);
    [K, ~] = size(edge_data); % get the number of edges in the graph
    c_s = sparse(zeros([num_paths, K])); % initialize
    
    for ind_of_path = 1:num_paths % iterate over each path
        
        for ind_in_path = 2:length(paths(ind_of_path, :)) % iterate through path, starting right after start node
            from_ID = paths(ind_of_path, ind_in_path-1); 
            to_ID = paths(ind_of_path, ind_in_path);
            edge_ind = intersect(find(edge_data(:,1) == from_ID), find(edge_data(:,2) == to_ID)); % get the index of the edge that connects the nodes
            c_s(ind_of_path, edge_ind) = 1; 
        end % end for 
    end % end for
end
       




% corresponds to eq. 5 in Time and Location source identification
% means = matrix describing means of each edge
% variances = matrix describing variances of each edge
% ti = time observed at contam node
% ts = guessed time observed at source node
% tau = uncertainty; 'window' taken from gaussian
% p = output probability from equation 5
function p = arrival_time_prob(ts, ti, path, dists, tau)
    travel_time = ti - ts;
    
    [mu_si, std_dev_si] = path_prob_stats(path, dists);
    
%     mu_si = 0;
%     var_si = 0;
%     % compute total mean, variance for normal representing entire path
%     for i = 1:(length(path)-1)
%         mu_si = mu_si + means(path(i), path(i+1)); % add mean, var for each edge
%         var_si = var_si + double(variances(path(i), path(i+1)));
%     end % end for 
    
    % equation 5
    p = normcdf(travel_time + tau, mu_si, std_dev_si) - normcdf(travel_time - tau, mu_si, std_dev_si);
end % end function

function [mean, std_dev] = path_prob_stats(path, dists)
    % get the path time
    total_dist = 0; % initialize total path dist
    for i = 2:length(path) % for each edge in the path
        total_dist = total_dist + dists(path(i-1), path(i)); % deterministic time component 
    end % end for
    
    mean = total_dist/630; % deterministic distance travel component
    determ_variance = (0.25*mean)*(0.25*mean); % according to "propagatio ratio"
    mean = mean + 8.5; % add offset mu_ab as given in travel time densities specs
    % truncate mean appropriately
    %mean = max(mean, total_dist/20);
    %mean = min(mean, total_dist/80);
    std_dev = sqrt(determ_variance + 1.5*1.5 + 3*3 + 1*1); % + 15.25 % likewise add offset, also convert from variance to std_dev
end % end function


% start_node = ID of path start node
% end_node = ID of path end node
% flows = the adjacency matrix of the network we're operating over
% path = [start_node ... ... end_node]; where path(i) is the i-th node in
% the max probability path from start to end
function [prob, path] = max_prob_path(start_node, end_node, flows)
    % transform to reduce the problem to one of shortest-path searching
    transf_flows = 1 - log(flows);
    transf_flows(transf_flows == Inf) = 0; % replace every resulting in Inf with 0
    transf_flows_bg = biograph(transf_flows);

    [prob, path, ~] = shortestpath(transf_flows_bg, start_node, end_node, 'Method', 'Acyclic');
    prob = exp(-(prob - length(path)+1)); 
end 


% node_set = array of node IDs that we are querying for the overall set of ancestors of. Assumed to all be in same stage.
% flows = the adjacency matrix
% prents = an array containing all the IDs of parents of node_set
function prents = parents(node_set, flows)
    prents = [];
    for j = 1:length(node_set)
        prents = union(prents, find(flows(:, node_set(j)))); % find each node that leads into nodes in j
        [r, c] = size(prents);
        if r > c % make sure it comes out as a row vector
            prents = prents.';
        end % end if
    end % end for
end % end ancestors
