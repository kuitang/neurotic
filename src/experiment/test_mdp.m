%% Prior!
clear mdp

[misc_data.segments, misc_data.slic_inds, misc_data.edge_G] = ...
    make_diff_graph(f1, 10, 0.1);

mdp = make_mdp_prior(PHI, misc_data);

%% Run!
run_mdp(@neal3_iter, [], [], PHI, f1, mdp, 1000, 0, 1);
