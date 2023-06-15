classdef PolarIntervalTest < matlab.unittest.TestCase
    % Tests for the ciat.PolarInterval class

    properties
        % Define several polar intervals of different sizes on which to perform tests
        
        % 1x1 intervals
        polarInt1 = ciat.PolarInterval(1, 2, 0, pi/2);
        polarInt2 = ciat.PolarInterval(3, 4, pi/4, pi/2);
        polarInt3 = ciat.PolarInterval(5, 6, pi/6, pi/2);
        polarInt4 = ciat.PolarInterval(7, 8, pi/2, pi/2);

        % 1x2 intervals

        % 2x1 intervals

        % 2x2 intervals
    end
    
    methods(Test)

        % function constructor(testCase)
        % end

        function abs(testCase)
            testCase.verifyEqual(abs(testCase.polarInt1), ciat.RealInterval(1, 2));
        end

        function angle(testCase)
            testCase.verifyEqual(angle(testCase.polarInt1), ciat.RealInterval(0, pi/2));
        end

        % function real(testCase)
        % end

        % function imag(testCase)
        % end

        % function area(testCase)
        % end

        function eq(testCase)
            testCase.verifyEqual(testCase.polarInt1 == testCase.polarInt1, true);
            testCase.verifyEqual(testCase.polarInt1 == testCase.polarInt2, false);
        end

        function ne(testCase)
            testCase.verifyEqual(testCase.polarInt1 ~= testCase.polarInt1, false);
            testCase.verifyEqual(testCase.polarInt1 ~= testCase.polarInt2, true);
        end

        % function cast(testCase)
        % end


        % function union(testCase)
        % end

        % function intersection(testCase)
        % end

        % function uminus(testCase)
        % end

        function times(testCase)
            testCase.verifyEqual(testCase.polarInt1 .* testCase.polarInt2, ciat.PolarInterval(3, 8, pi/4, pi));
            testCase.verifyEqual(testCase.polarInt1 .* testCase.polarInt3, ciat.PolarInterval(5, 12, pi/6, pi));
            testCase.verifyEqual(testCase.polarInt1 .* testCase.polarInt4, ciat.PolarInterval(7, 16, pi/2, pi));
            testCase.verifyEqual(testCase.polarInt3 .* testCase.polarInt4, ciat.PolarInterval(35, 48, pi/2 + pi/6, pi));

            testCase.verifyEqual(2 .* testCase.polarInt1, ciat.PolarInterval(2, 4, 0, pi/2));
            testCase.verifyEqual(testCase.polarInt1 .* (2*exp(1j*pi/2)), ciat.PolarInterval(2, 4, pi/2, pi));
        end

        % function mtimes(testCase)
        % end
    end    
end
