%% Prior!
clear mdp
mdp = make_mdp_prior(PHI1, misc_data);

%% Run!
run_mdp(@neal3_iter, [], [], PHI1, f1, mdp, 1000, 0, 1);
