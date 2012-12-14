%% Prior!
clear mdp
clear misc_data
misc_data = struct();
[misc_data.segments, misc_data.slic_inds, misc_data.edge_G] = ...
    make_diff_graph(f1, 10, 0.1);
misc_data
mdp = make_mdp_prior(PHI1, misc_data);

%% Run!
run_mdp(@neal3_iter, [], [], PHI1, f1, mdp, 1000, 0, 1);
