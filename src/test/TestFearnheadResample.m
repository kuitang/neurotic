classdef TestFearnheadResample < TestCase

    % Generate some fixture data here...
    properties
        N = 100;
        N_test_iters = 10;
        expectation_tol = 0.05;
        weights;
    end

    % Tests are methods
    methods

        % Constructor. Takes name of test method and passes to the base constructor
        function self = TestFearnheadResample(name)
            % This is the super() syntax in MATLAB.
            self = self@TestCase(name);
        end

        % Now the tests
        function setUp(self)
            self.weights = rand(self.N, 1);
            self.weights = self.weights / sum(self.weights);
        end

        function tearDown(self)
            
        end
        
        function testStratifiedResampleEdgeCases(self)
            % Error
            zero_idxs = logical(zeros(size(self.weights)));
            sampling_nothing = @() stratified_resample(self.weights, zeros(size(self.weights)), 1);
            assertExceptionThrown(sampling_nothing, 'srs:subsetTooSmall');
            
            % Error: wrong indexing size
            wrong_dims = @() stratified_resample(self.weights, 1, 1);
            assertExceptionThrown(wrong_dims, 'srs:wrongDims');
            
            % Should return no indices
            resample_nothing_nothing = stratified_resample(self.weights, zero_idxs, 0);
            assertTrue(sum(resample_nothing_nothing) == 0);
            
            one_idxs = logical(zeros(size(self.weights)));
            one_idxs(1) = 1;
            
            resample_nothing_something = stratified_resample(self.weights, one_idxs, 0);
            assertTrue(sum(resample_nothing_something) == 0);
            
            % Should return [1 0 ... 0]
            resample_one_one = stratified_resample(self.weights, one_idxs, 1);
            assertTrue(all(resample_one_one == one_idxs))            
        end
        
        function testStratifiedResample(self)
            % Postconditions:
            %
            % N - L particles are resampled
            % E[number of times particle i resampled] \propto weight of i            
            
            n_samples = 1000;  
            
            for i = 1:self.N_test_iters
                subsetidxs     = logical(zeros(size(self.weights)));
                nnz_subsetidxs = 2 + randi(self.N - 2);
                n_particles    = floor(nnz_subsetidxs / 2); % so we can compute expectations in reasonable time
                
                inds = randsample(self.N, nnz_subsetidxs);
                subsetidxs(inds) = true;
                
                % Sanity check
                assert(sum(subsetidxs) == nnz_subsetidxs);
                assert(n_particles <= nnz_subsetidxs);
                                                                
                % Each sample is a column
                sampleidxs = logical(zeros(self.N, n_samples));
                for n = 1:n_samples
                    sampleidxs(:,n) = stratified_resample(self.weights, subsetidxs, n_particles);
                    
                    % n_particles particles are actually resampled
                    n_actual_particles = sum(sampleidxs(:,n));
                    assertTrue(abs(n_actual_particles - n_particles) <= 1);
                end
                
                % sampledixs is just an indicator matrix. Select the actual
                % samples indicated by multiplication.
                rep_weights = repmat(self.weights, 1, n_samples);                                
                sample_weights = rep_weights .* sampleidxs;                                
                
                % The statistical property is:
                %   E[ Q_i ] = q_i
                % 
                % Where Q is the random variable obtained by resampling.
                
                mean_weights = mean(sample_weights, 2);
                
                % You need some smarter way to diagnose this, because most
                % of the entries are sparse.
                assertElementsAlmostEqual(self.weights, mean_weights, 'absolute', self.expectation_tol);
            end
            
        end
                
        function testFearnheadEdgeCases(self)
            % Sampling something from nothing should error
            fail = @() fearnhead_resample([], 1);
            assertExceptionThrown(fail, 'frs:tooManyParticles');
            
            % Sampling nothing from nothing should error (because the
            % optimization equation  is not well-defined)
            zero_weights = zeros(size(self.weights));
            fail2 = @() fearnhead_resample(zero_weights, 0);
            assertExceptionThrown(fail2, 'frs:zeroParticles');            

            one_weights = zeros(size(self.weights));
            one_idxs(1) = 1;
            
            % Sampling nothing from something should return empty
            fail3 = @() fearnhead_resample(one_weights, 0);
            assertExceptionThrown(fail2, 'frs:zeroParticles');
            
            % Sampling one from one should return one
            assertTrue(all(isequal(fearnhead_resample(one_weights, 1), one_weights)));            
        end
        
        function testFearnheadStatistics(self)
            % Tests the conditions of Theorem 1 in [1]
            
            n_samples = 100;
            %for i = 1:self.N_test_iters
            for i = 1:10
                n_particles = floor(self.N / 10);
                
                new_weights = zeros(self.N, n_samples);
                for n = 1:n_samples
                    new_weights(:, n) = fearnhead_resample(self.weights, n_particles);
                                        
                    n_actual_particles = sum(new_weights(:, n) > 0);                    
                    assertTrue(abs(n_actual_particles - n_particles) <= 1);
                end
            end
            
            mean_new_weights = mean(new_weights, 2);
            
            assertElementsAlmostEqual(self.weights, mean_new_weights, 'absolute', self.expectation_tol);            
        end
    end
end
