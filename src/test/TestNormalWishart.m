classdef TestNormalWishart < TestCase

    % Initialized by setup and used by methods (instead of passing around parameters)
    % Subclassing is useful when you have a lot of state you need to test.
    properties
        nw x D;        
    end

    % Tests are methods
    methods

        % Constructor. Takes name of test method and passes to the base constructor
        function self = TestNormalWishart(name)
            % This is the super() syntax in MATLAB.
            self = self@TestCase(name);
        end

        % Now the tests
        function setUp(self)
            self.nw = NormalWishart([2 -1 1 -2], eye(4), 5, 1);  
            self.x  = [10 -5 2 6];
            self.D  = 4;
        end

        function tearDown(self)
            
        end

        function assertAlmostEqualTo(self, other)
            assertElementsAlmostEqual(self.nw.post_mean, other.post_mean);
            assertElementsAlmostEqual(self.nw.post_dof,  other.post_dof);
            
            assertElementsAlmostEqual(self.nw.pred_mean, other.pred_mean);
            assertElementsAlmostEqual(self.nw.pred_dof,  other.pred_dof);

            assertElementsAlmostEqual(self.nw.post_chol, other.post_chol);
            assertElementsAlmostEqual(self.nw.pred_chol, other.pred_chol);
        end
        
        function testEdgeCases(self)
            self.nw.fit([]);            
            % Should revert to the prior
            assertEqual(self.nw.post_mean, self.nw.prior_mean);
            assertEqual(self.nw.post_n, self.nw.prior_n);
            
            post_cov = self.nw.post_chol' * self.nw.post_chol;
            assertElementsAlmostEqual(post_cov, self.nw.prior_cov);
                        
            self.nw.fit(self.x);              
            % TODO: Do some tests!
        end
        
        function testOnlineEdgeCases(self)
            self.nw.add_point(self.x);
            
            % TODO: Update tests!            
            
            self.nw.remove_point(self.x);
            
            % TODO: Update tests!

        end
        
        function testOnlineAlmostEqual(self)
            % Make sure the online updates and offline updates return the
            % same results
            N = 100;
            X = mvnrnd([2 0 5 -1], 2*eye(4), N);
            
            other = copy(self.nw);
            
            self.nw.fit(X);
                        
            for n = 1:N                
                other.add_point(X(n,:));
            end
                        
            self.assertAlmostEqualTo(other);            
            % Make sure both predictors return the same results
            assertElementsAlmostEqual(self.nw.pred_like(X), other.pred_like(X));                        
            
            % Now test removals
            X1 = X(1,:);
            self.nw.fit(X1);
            
            for n = 2:N
                other.remove_point(X(n,:));
            end
            
            self.assertAlmostEqualTo(other);
            assertElementsAlmostEqual(self.nw.pred_like(X), other.pred_like(X));
            
        end

    end
end
