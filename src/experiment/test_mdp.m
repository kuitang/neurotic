%% Prior!
clear mdp
mdp = make_mdp_prior(PHI2);

%% Run!
run_mdp(@neal3_iter, [], [], PHI2, f1, mdp, 1000, 0, 1);
