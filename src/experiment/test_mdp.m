%% Prior!
clear mdp

[misc_data.segments, misc_data.slic_inds, misc_data.edge_G] = make_diff_graph(f1, 10, 0.1);

[mdp, init_state, params]  = make_mdplight_prior(PHI, misc_data);

%% Run!
run_mdp(@neal8_iter, init_state, params, PHI, f1, mdp, 1000, 0, 1);
