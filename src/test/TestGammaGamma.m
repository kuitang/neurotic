classdef TestGammaGamma < TestCase
    
    properties       
       % A default gg, for the more mechnically tests when we don't care
       % so much about inferential accuracy
       gg;
       
       % Some people believe in testing randomized algorithms with a fixed
       % seed. I am not one of those people. Run the fitting tests for
       % multiple iterations instead.
       N_fit_iters = 1000;
       N_post_samples = 50;
       
       % Absolute tolerance (because there are probabilities, and can get
       % close to zero, where relative tolerance becomes problematic)
       pred_tol = 0.1;
       
       % This is a RELATIVE tolerance. So the two values being compared
       % must be within 10%.
       param_tol = 0.1;
    end

    % Tests are methods
    methods
        
        function self = TestGammaGamma(name)
            % This is the super() syntax in MATLAB.
            self = self@TestCase(name);
        end

        % Now the tests
        function setUp(self)
            self.gg = GammaGamma(1, 2, 1);
        end

        function tearDown(self)
            
        end                
        
        function helpTestFits(self, gg, Npoints, shape, true_rate)
            for i = 1:self.N_fit_iters
                X = gamrnd(shape, 1/true_rate, Npoints, 1);
                gg.fit(X);

                % Draw some samples from the posterior distribution. Their
                % mean should be close to the true means.
                mean_post_rate = mean(gamrnd(gg.post_shape, 1/gg.post_rate, self.N_post_samples, 1));
                assertElementsAlmostEqual(mean_post_rate, true_rate, 'relative', self.param_tol);
                
                % The predictive likelihood should be similar to the true
                % likelihoods                
                data_likes = gampdf(X, shape, 1/true_rate);
                pred_likes = gg.pred_like(X);                                
                
                assertElementsAlmostEqual(data_likes, pred_likes, 'absolute', self.pred_tol); 
            end    
        end
        
        function testFitsCorrectPrior(self)
            % The prior rate has mean 1 and variance 1. Generate data from
            % Gamma(shape=2, rate=1). We expect that gg.pred_like(X) will
            % have absolute difference at most tol with respect to
            % gampdf(X, 2, 1).
                        
            gg  = GammaGamma(2, 1, 1);            
            self.helpTestFits(gg, 1000, 2, 1);
        end
        
        function testFitsWrongPrior(self)
            % Now we will set a prior rate with mean 5 variance 10, which
            % comes to shape = 2.5, rate = 0.5 (since E[X] = shape/rate and
            % V[X] = shape/rate^2; a little algebra does the rest).
            %
            % Our data rate will actually be 10. (This is 0.1 in MATLAB).
            %
            % In addition, we will set our model shape to 3 but set the
            % data shape to 2.9. The Gamma-Gamma model is not robust to shape
            % changes. (Experimentally, even 2.9 fails over 1000 trials).
            %
            % Note that the Gamma model is not particularly robust to shape
            % changes.
            %
            % We expect
            % a large amount of data (1000 points) to overwhelm the prior,
            % particularly since it is quite broad (the variance is 10).
            
            gg  = GammaGamma(3, 2.5, 0.5);
            self.helpTestFits(gg, 10000, 3, 10);            
        end                

        function assertAlmostEqualTo(self, other)
            % Should really pull this out into a foreach...
            fns = properties(self.gg);
            for i=1:length(fns)
                assertElementsAlmostEqual(self.gg.(fns{i}), other.(fns{i}));                
            end
        end
        
        function testEdgeCases(self)
            % The initialized class should be prior_* = post_*
            assertEqual(self.gg.prior_shape, self.gg.post_shape);
            assertEqual(self.gg.prior_rate,  self.gg.prior_rate);
            
            % Fitting empty data should revert to the prior
            clean_copy = copy(self.gg);
            self.gg.fit([]);
            self.assertAlmostEqualTo(clean_copy);
            
            % Removing from empty data should revert to the prior
            self.gg.remove_point(1);
            self.assertAlmostEqualTo(clean_copy);
            
            % Adding a point and then removing it should revert to the
            % prior
            self.gg.remove_point(1);
            self.assertAlmostEqualTo(clean_copy);
        end        
        
        function testOnlineAlmostEqual(self)
            % Make sure the online updates and batch updates return the
            % same results
            N = 100;
            X = gamrnd(5, 5, N, 1); % doesn't really matter
            
            online = copy(self.gg);
            
            self.gg.fit(X);
                        
            for n = 1:N                
                online.add_point(X(n));
            end
                        
            self.assertAlmostEqualTo(online);
            
            % Make sure both predictors return the same results
            assertElementsAlmostEqual(self.gg.pred_like(X), online.pred_like(X));                        
            
            % Now test removals
            x1 = X(1,:);
            self.gg.fit(x1);
            
            for n = 2:N
                online.remove_point(X(n));
            end
            
            self.assertAlmostEqualTo(online);
            assertElementsAlmostEqual(self.gg.pred_like(X), online.pred_like(X));
            
        end

    end
end
