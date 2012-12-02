classdef TestProductDistribution < TestCase

    % Initialized by setup and used by methods (instead of passing around parameters)
    % Subclassing is useful when you have a lot of state you need to test.
    properties
        nw x D;        
    end

    % Tests are methods
    methods

        % Constructor. Takes name of test method and passes to the base constructor
        function self = TestProductDistribution(name)
            % This is the super() syntax in MATLAB.
            self = self@TestCase(name);
        end

        % Now the tests
        function setUp(self)

        end

        function tearDown(self)
            
        end

        function assertAlmostEqualTo(self, other)

        end

        function testIgnorance(self)
            % Configuring only a subset of indices should 
        end
        
        function testFitForwarding(self)
            % Test that fitting (on and offline) are correctly forwarded to
            % the underlying distributions. Make sure to test the
            % underlying distributions separately.            
        end
        
        function testProduct(self)
            % Test correctness of predictive likelihoods.
        end                    

    end
end
