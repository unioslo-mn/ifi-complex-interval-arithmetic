classdef CircularIntervalTest < matlab.unittest.TestCase
    % Tests for the ciat.CircularInterval class

    properties
        % Define several circular intervals of different sizes on which to perform tests
        
        % 1x1 intervals
        circInt1 = ciat.CircularInterval(0, 1);
        circInt2 = ciat.CircularInterval(1, 2);
        circInt3 = ciat.CircularInterval(2, 3);
        
        % 1x2 intervals
        circInt4 = ciat.CircularInterval([0, 2], [1, 3]);
        circInt5 = ciat.CircularInterval([1, 3], [2, 4]);
        
        % 2x1 intervals
        circInt6 = ciat.CircularInterval([0; 2], [1; 3]);
        circInt7 = ciat.CircularInterval([1; 3], [2; 4]);
        
        % 2x2 intervals
        circInt8 = ciat.CircularInterval([0, 2; 1, 3], [1, 3; 2, 4]);
        circInt9 = ciat.CircularInterval([1, 3; 2, 4], [2, 4; 3, 5]);
    end
    
    methods(Test)
        
        function real(testCase)
            testCase.verifyEqual(real(testCase.circInt1), ciat.RealInterval(-1, 1));
            testCase.verifyEqual(real(testCase.circInt2), ciat.RealInterval(-1, 3));
            testCase.verifyEqual(real(testCase.circInt3), ciat.RealInterval(-1, 5));
            testCase.verifyEqual(real(testCase.circInt4), ciat.RealInterval([-1, -1], [1, 5]));
            testCase.verifyEqual(real(testCase.circInt5), ciat.RealInterval([-1, -1], [3, 7]));
            testCase.verifyEqual(real(testCase.circInt6), ciat.RealInterval([-1; -1], [1; 5]));
            testCase.verifyEqual(real(testCase.circInt7), ciat.RealInterval([-1; -1], [3; 7]));
            testCase.verifyEqual(real(testCase.circInt8), ciat.RealInterval([-1, -1; -1, -1], [1, 5; 3, 7]));
            testCase.verifyEqual(real(testCase.circInt9), ciat.RealInterval([-1, -1; -1, -1], [3, 7; 5, 9]));
        end
        
        % function imag(testCase)
        %     testCase.verifyEqual(imag(testCase.circInt1), ciat.RealInterval(-1, 1));
        %     testCase.verifyEqual(imag(testCase.circInt4), ciat.RealInterval([-1, 1], [-1, 1]));
        %     testCase.verifyEqual(imag(testCase.circInt6), ciat.RealInterval([-1; 1], [-1; 1]));
        %     testCase.verifyEqual(imag(testCase.circInt8), ciat.RealInterval([-1, 1; -1, 1], [-1, 1; -1, 1]));
        % end
        
        function eq(testCase)
            testCase.verifyEqual(testCase.circInt1 == testCase.circInt1, true);
            testCase.verifyEqual(testCase.circInt1 == testCase.circInt2, false);
            testCase.verifyEqual(testCase.circInt2 == testCase.circInt3, false);

            testCase.verifyEqual(testCase.circInt4 == testCase.circInt4, true);
            testCase.verifyEqual(testCase.circInt4 == testCase.circInt5, false);

            testCase.verifyEqual(testCase.circInt6 == testCase.circInt6, true);
            testCase.verifyEqual(testCase.circInt6 == testCase.circInt7, false);
            
            testCase.verifyEqual(testCase.circInt8 == testCase.circInt8, true);
            testCase.verifyEqual(testCase.circInt8 == testCase.circInt9, false);
        end
        
        function ne(testCase)
            testCase.verifyEqual(testCase.circInt1 ~= testCase.circInt1, false);
            testCase.verifyEqual(testCase.circInt1 ~= testCase.circInt2, true);
            testCase.verifyEqual(testCase.circInt2 ~= testCase.circInt3, true);

            testCase.verifyEqual(testCase.circInt4 ~= testCase.circInt4, false);
            testCase.verifyEqual(testCase.circInt4 ~= testCase.circInt5, true);

            testCase.verifyEqual(testCase.circInt6 ~= testCase.circInt6, false);
            testCase.verifyEqual(testCase.circInt6 ~= testCase.circInt7, true);
            
            testCase.verifyEqual(testCase.circInt8 ~= testCase.circInt8, false);
            testCase.verifyEqual(testCase.circInt8 ~= testCase.circInt9, true);
        end
        
        % function cast(testCase)
        %     testCase.verifyEqual(ciat.CircularInterval(testCase.circInt1), ciat.CircularInterval(0, 1, 0, 0));
        %     testCase.verifyEqual(ciat.CircularInterval(testCase.circInt4), ciat.CircularInterval([0, 2], [1, 3], [0, 0], [0, 0]));
        %     testCase.verifyEqual(ciat.CircularInterval(testCase.circInt6), ciat.CircularInterval([0; 2], [1; 3], [0; 0], [0; 0]));
        %     testCase.verifyEqual(ciat.CircularInterval(testCase.circInt8), ciat.CircularInterval([0, 2; 1, 3], [1, 3; 2, 4], [0, 0; 0, 0], [0, 0; 0, 0]));
        % end
        
        % function plus(testCase)
        %     testCase.verifyEqual(testCase.circInt1 + testCase.circInt2, ciat.CircularInterval(1, 3, 2, 4));
        %     testCase.verifyEqual(testCase.circInt1 + testCase.circInt3, ciat.CircularInterval(2, 4, 3, 5));
        %     testCase.verifyEqual(testCase.circInt2 + testCase.circInt3, ciat.CircularInterval(3, 5, 4, 6));
        %     testCase.verifyEqual(testCase.circInt3 + testCase.circInt4, ciat.CircularInterval(2, 4, 4, 6));
        %     testCase.verifyEqual(testCase.circInt4 + testCase.circInt5, ciat.CircularInterval([1, 5], [2, 6], [2, 6], [4, 8]));
        %     testCase.verifyEqual(testCase.circInt5 + testCase.circInt6, ciat.CircularInterval([1, 5], [2, 6], [2, 6], [4, 8]));
        %     testCase.verifyEqual(testCase.circInt6 + testCase.circInt7, ciat.CircularInterval([1; 5], [2; 6], [2; 6], [4; 8]));
        %     testCase.verifyEqual(testCase.circInt7 + testCase.circInt8, ciat.CircularInterval([1; 5], [2; 6], [2; 6], [4; 8]));
        %     testCase.verifyEqual(testCase.circInt8 + testCase.circInt9, ciat.CircularInterval([1, 5; 3, 7], [2, 6; 4, 8], [2, 6; 4, 8], [4, 8; 6, 10]));
        % end
        
        % function minus(testCase)
        %     testCase.verifyEqual(testCase.circInt1 - testCase.circInt2, ciat.CircularInterval(-1, 1, -2, 0));
        %     testCase.verifyEqual(testCase.circInt1 - testCase.circInt3, ciat.CircularInterval(0, 2, 2, 4));
        %     testCase.verifyEqual(testCase.circInt2 - testCase.circInt3, ciat.CircularInterval(1, 3, 2, 4));
        %     testCase.verifyEqual(testCase.circInt3 - testCase.circInt4, ciat.CircularInterval(-1, 1, 2, 4));
        %     testCase.verifyEqual(testCase.circInt4 - testCase.circInt5, ciat.CircularInterval([-3, 1], [-2, 2], [-2, 2], [0, 4]));
        %     testCase.verifyEqual(testCase.circInt5 - testCase.circInt6, ciat.CircularInterval([-3, 1], [-2, 2], [-2, 2], [0, 4]));
        %     testCase.verifyEqual(testCase.circInt6 - testCase.circInt7, ciat.CircularInterval([-3; 1], [-2; 2], [-2; 2], [0; 4]));
        %     testCase.verifyEqual(testCase.circInt7 - testCase.circInt8, ciat.CircularInterval([-3; 1], [-2; 2], [-2; 2], [0; 4]));
        %     testCase.verifyEqual(testCase.circInt8 - testCase.circInt9, ciat.CircularInterval([-3, 1; -1, 3], [-2, 2; 0, 4], [-2, 2; 0, 4], [0, 4; 2, 6]));
        % end
        
        % function times(testCase)
        %     testCase.verifyEqual(testCase.circInt1 .* testCase.circInt2, ciat.CircularInterval(0, 2, 0, 4));
        %     testCase.verifyEqual(testCase.circInt1 .* testCase.circInt3, ciat.CircularInterval(0, 2, 0, 4));
        %     testCase.verifyEqual(testCase.circInt2 .* testCase.circInt3, ciat.CircularInterval(0, 4, 0, 8));
        %     testCase.verifyEqual(testCase.circInt3 .* testCase.circInt4, ciat.CircularInterval(0, 2, 0, 4));
        %     testCase.verifyEqual(testCase.circInt4 .* testCase.circInt5, ciat.CircularInterval([-4, 4], [-6, 6], [-6, 6], [-12, 12]));
        %     testCase.verifyEqual(testCase.circInt5 .* testCase.circInt6, ciat.CircularInterval([-4, 4], [-6, 6], [-6, 6], [-12, 12]));
        %     testCase.verifyEqual(testCase.circInt6 .* testCase.circInt7, ciat.CircularInterval([-4; 4], [-6; 6], [-6; 6], [-12; 12]));
        %     testCase.verifyEqual(testCase.circInt7 .* testCase.circInt8, ciat.CircularInterval([-4; 4], [-6; 6], [-6; 6], [-12; 12]));
        %     testCase.verifyEqual(testCase.circInt8 .* testCase.circInt9, ciat.CircularInterval([-4, 4; -2, 2], [-6, 6; -4, 4], [-6, 6; -4, 4], [-12, 12; -8, 8]));
        % end
        
        % function mtimes(testCase)
        %     testCase.verifyEqual(testCase.circInt1 * testCase.circInt2, ciat.CircularInterval(0, 2, 0, 4));
        %     testCase.verifyEqual(testCase.circInt1 * testCase.circInt3, ciat.CircularInterval(0, 2, 0, 4));
        %     testCase.verifyEqual(testCase.circInt2 * testCase.circInt3, ciat.CircularInterval(0, 4, 0, 8));
        %     testCase.verifyEqual(testCase.circInt3 * testCase.circInt4, ciat.CircularInterval(0, 2, 0, 4));
        %     testCase.verifyEqual(testCase.circInt4 * testCase.circInt5, ciat.CircularInterval([-4, 4], [-6, 6], [-6, 6], [-12, 12]));
        %     testCase.verifyEqual(testCase.circInt5 * testCase.circInt6, ciat.CircularInterval([-4, 4], [-6, 6], [-6, 6], [-12, 12]));
        %     testCase.verifyEqual(testCase.circInt6 * testCase.circInt7, ciat.CircularInterval([-4; 4], [-6; 6], [-6; 6], [-12; 12]));
        %     testCase.verifyEqual(testCase.circInt7 * testCase.circInt8, ciat.CircularInterval([-4; 4], [-6; 6], [-6; 6], [-12; 12]));
        %     testCase.verifyEqual(testCase.circInt8 * testCase.circInt9, ciat.CircularInterval([-4, 4; -2, 2], [-6, 6; -4, 4], [-6, 6; -4, 4], [-12, 12; -8, 8]));
        % end
        
        % function rdivide(testCase)
        %     testCase.verifyEqual(testCase.circInt1 ./ testCase.circInt2, ciat.CircularInterval(0, 1, 0, 0.5));
        %     testCase.verifyEqual(testCase.circInt1 ./ testCase.circInt3, ciat.CircularInterval(0, 1, 0, 0.5));
        %     testCase.verifyEqual(testCase.circInt2 ./ testCase.circInt3, ciat.CircularInterval(0, 2, 0, 1));
        %     testCase.verifyEqual(testCase.circInt3 ./ testCase.circInt4, ciat.CircularInterval(0, 1.5, 0, 0.75));
        %     testCase.verifyEqual(testCase.circInt4 ./ testCase.circInt5, ciat.CircularInterval([-2, 2], [-2, 2], [-2, 2], [-1, 1]));
        %     testCase.verifyEqual(testCase.circInt5 ./ testCase.circInt6, ciat.CircularInterval([-2, 2], [-2, 2], [-2, 2], [-1, 1]));
        %     testCase.verifyEqual(testCase.circInt6 ./ testCase.circInt7, ciat.CircularInterval([-2; 2], [-2; 2], [-2; 2], [-1; 1]));
        %     testCase.verifyEqual(testCase.circInt7 ./ testCase.circInt8, ciat.CircularInterval([-2; 2], [-2; 2], [-2; 2], [-1; 1]));
        %     testCase.verifyEqual(testCase.circInt8 ./ testCase.circInt9, ciat.CircularInterval([-2, 2; -1.5, 1.5], [-2, 2; -1, 1], [-2, 2; -1, 1], [-1, 1; -0.5, 0.5]));
        % end
        
        % function mrdivide(testCase)
        %     testCase.verifyEqual(testCase.circInt1 / testCase.circInt2, ciat.CircularInterval(0, 1, 0, 0.5));
        %     testCase.verifyEqual(testCase.circInt1 / testCase.circInt3, ciat.CircularInterval(0, 1, 0, 0.5));
        %     testCase.verifyEqual(testCase.circInt2 / testCase.circInt3, ciat.CircularInterval(0, 2, 0, 1));
        %     testCase.verifyEqual(testCase.circInt3 / testCase.circInt4, ciat.CircularInterval(0, 1.5, 0, 0.75));
        %     testCase.verifyEqual(testCase.circInt4 / testCase.circInt5, ciat.CircularInterval([-2, 2], [-2, 2], [-2, 2], [-1, 1]));
        %     testCase.verifyEqual(testCase.circInt5 / testCase.circInt6, ciat.CircularInterval([-2, 2], [-2, 2], [-2, 2], [-1, 1]));
        %     testCase.verifyEqual(testCase.circInt6 / testCase.circInt7, ciat.CircularInterval([-2; 2], [-2; 2], [-2; 2], [-1; 1]));
        %     testCase.verifyEqual(testCase.circInt7 / testCase.circInt8, ciat.CircularInterval([-2; 2], [-2; 2], [-2; 2], [-1; 1]));
        %     testCase.verifyEqual(testCase.circInt8 / testCase.circInt9, ciat.CircularInterval([-2, 2; -1.5, 1.5], [-2, 2; -1, 1], [-2, 2; -1, 1], [-1, 1; -0.5, 0.5]));
        % end 
    end
end