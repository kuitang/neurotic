%% Generate data
prior_params = [ 0.1 0.1 0.0005 ];
true_m = 0.8;
true_s = 0.9;

[a b] = beta_sm_to_ab(true_s, true_m);
a = true_s * true_m;
b = true_s * (1 - true_m);

SN = 1000;

xs = betarnd(a * ones(SN,1), b * ones(SN, 1));

%% Initialize with MLE
mle_params = betafit(xs);
[s m] = beta_ab_to_sm(mle_params(1), mle_params(2))

m_star = interval_to_positive(m);

%% Sample!
Nsamps = 3000;
Nburn_in = 1000;
samps = zeros(Nsamps, 2);

for n = 1:Nburn_in
    [s m_star] = beta_posterior_rnd(prior_params, s, m_star, xs, 0.01);
end

samps(1,:) = [s m_star];

Naccept = 0;
Nreject = 0;
for n = 2:Nsamps
    s_old      = samps(n - 1,1);
    m_star_old = samps(n - 1,2);
    [s_new m_star_new acc] = beta_posterior_rnd(prior_params, s_old, m_star_old, xs, 0.01);
    samps(n,:) = [s_new m_star_new];
    if acc        
        Naccept = Naccept + 1;
    end
    if ~acc
        Nreject = Nreject + 1;
    end    
end
disp(['Rejections: ' num2str(Nreject)]);

accept_rate = Naccept / Nsamps;
if accept_rate > 0.9
    warning(['Unusually high accept rate: ' num2str(accept_rate)]);
end

disp(['Accept rate: ', num2str(Naccept / Nsamps)]);

figure(1);
hist(samps(:,1));
title('s samples');

samps(:,2) = positive_to_interval(samps(:,2));

figure(2);
hist(samps(:,2));
title('m samples');