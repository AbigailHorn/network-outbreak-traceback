# Locating the Sources of Large-Scale Outbreaks on Networks

![alt tag](https://raw.githubusercontent.com/elenapolozova/network-outbreak-traceback/master/utility/regional_network_pic.png)

One Paragraph of project description goes here

## Getting Started

Download all files and check out the ‘examples’ folder! heuristic_tests.m walks through the whole process of network generation, location assignment, outbreak simulation, and traceback from start to finish. heuristic_vs_exact.m contains a tractable example outbreak for the exact estimators.

Documentation for each utility function can be found in comments at the top of the files.

### What’s Included

* Network connectivity generation

```
[flows, node_layers, init_vols, stage_ends] = random_layered_graph(node_counts, 
    init_vol_dist, network_params, show_plots, display_network);
```

* Assignment of network nodes to physical locations

```
[dists, node_locs] = assign_locations(node_layers, flows, show_loc_plots);
```

* Simulated outbreak propagation on generated networks

```
[contam_farm, contam_retailers] = outbreak(dispersion, src_farm, dists,
    node_locs, flows, init_vols, num_stages, return_all_reports, 
    plot_or_not, transport_dev_frac, storage_is_determ);
```

* Traceback algorithm over resulting outbreak information

```
[pmf, t_s_stars] = time_traceback(method_id, flows, stage_ends, init_vols,
    contam_retailers, distances, P, transport_dev_frac);

```


## Authors

* **Elena Polozova** - *MATLAB Implementation* - 
* **Abigail Horn**  - *Theoretical Derivations* -


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details



