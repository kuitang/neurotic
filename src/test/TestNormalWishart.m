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
            assertElementsAlmostEqual(self.nw.data_mean, other.data_mean);
            assertElementsAlmostEqual(self.nw.data_cov,  other.data_cov);
            
            assertElementsAlmostEqual(self.nw.pred_mean, other.pred_mean);
            assertElementsAlmostEqual(self.nw.pred_cov,  other.pred_cov);
            assertElementsAlmostEqual(self.nw.pred_dof,  other.pred_dof);
        end
        
        function testEdgeCases(self)
            self.nw.fit([]);            
            
            assertElementsAlmostEqual(self.nw.data_mean, zeros(1, self.D));
            assertElementsAlmostEqual(self.nw.data_cov,  zeros(self.D));
            
            self.nw.fit(self.x);
            
            assertElementsAlmostEqual(self.nw.data_mean, self.x);
            assertElementsAlmostEqual(self.nw.data_cov,  zeros(self.D));            
        end
        
        function testOnlineEdgeCases(self)
            self.nw.add_point(self.x);
            
            assertElementsAlmostEqual(self.nw.data_mean, self.x);
            assertElementsAlmostEqual(self.nw.data_cov,  zeros(self.D));
            
            self.nw.remove_point(self.x);
            
            assertElementsAlmostEqual(self.nw.data_mean, zeros(1, self.D));
            assertElementsAlmostEqual(self.nw.data_cov,  zeros(self.D));
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
