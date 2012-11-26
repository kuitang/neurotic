classdef TestNormalWishart < TestCase

    % Initialized by setup and used by methods (instead of passing around parameters)
    % Subclassing is useful when you have a lot of state you need to test.
    properties
        nw;
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
        end

        function tearDown(self)
            
        end

        function testOnlineUpdates(self)
            % Make sure the online updates and offline updates return the
            % same results
            N = 100;
            X = mvnrnd([2 2 2 2], 2*eye(4), N);
            
            other = copy(self.nw);
            
            self.nw.fit(X);
            
            for n = 1:N                
                other.add_point(X(n,:));
            end
                        
            assertElementsAlmostEqual(self.nw.data_mean, other.data_mean);
            assertElementsAlmostEqual(self.nw.data_cov,  other.data_cov);
            
            assertElementsAlmostEqual(self.nw.pred_mean, other.pred_mean);
            assertElementsAlmostEqual(self.nw.pred_cov,  other.pred_cov);
            assertElementsAlmostEqual(self.nw.pred_dof,  other.pred_dof);
            
            % Make sure both predictors return the same results
            assertElementsAlmostEqual(self.nw.pred_like(X), other.pred_like(X));            
        end

    end
end
