classdef RectangularIntervalTest < matlab.unittest.TestCase
    % Tests for the ciat.RectangularInterval class

    properties
        % Define several rectangular intervals of different sizes on which to perform tests
        
        % 1x1 intervals
        rectInt1 = ciat.RectangularInterval(1, 2, 3, 4);
        rectInt2 = ciat.RectangularInterval(3, 4, 5, 6);
        rectInt3 = ciat.RectangularInterval(-1, 1, -1, 1);
        rectInt4 = ciat.RectangularInterval(-2, -1, -2, -1);

        % 1x2 intervals
        rectInt5 = ciat.RectangularInterval([1, 3], [2, 4], [3, 6], [4, 8]);
        rectInt6 = ciat.RectangularInterval([3, 4], [4, 8], [5, 10], [6, 12]);

        % 2x1 intervals
        rectInt7 = ciat.RectangularInterval([1; 2], [3; 4], [5; 6], [7; 8]);
        rectInt8 = ciat.RectangularInterval([2; 3], [4; 5], [6; 7], [8; 9]);

        % 2x2 intervals
        rectInt9 = ciat.RectangularInterval([1, 2; 3, 4], [5, 6; 7, 8], [9, 10; 11, 12], [13, 14; 15, 16]);
        rectInt10 = ciat.RectangularInterval([17, 18; 19, 20], [21, 22; 23, 24], [25, 26; 27, 28], [29, 30; 31, 32]);


        % Real intervals

        % 1x1 intervals
        realInt1 = ciat.RealInterval(1,2);
        realInt2 = ciat.RealInterval(3,4);
        realInt3 = ciat.RealInterval(-1, 1);
        realInt4 = ciat.RealInterval(-2, -1)

        % 1x2 intervals
        realInt5 = ciat.RealInterval([1,3], [2,4]);
        realInt6 = ciat.RealInterval([3,4], [5,6]);

        % 2x1 intervals
        realInt7 = ciat.RealInterval([1;2], [3;4]);
        realInt8 = ciat.RealInterval([2;3], [4;5]);

        % 2x2 intervals
        realInt9 = ciat.RealInterval([1,2; 3,4], [5,6; 7,8]);
        realInt10 = ciat.RealInterval([2,4; 6,8], [10,12; 14,16]);

    end
    
    methods(Test)

        % function constructor(testCase)
            % Test different ways to construct a rectangular interval
            % 4 inputs
            % testCase.verifyEqual(ciat.RectangularInterval(1, 2, 3, 4).Real, [1, 2]);
            % testCase.verifyEqual(ciat.RectangularInterval(1, 2, 3, 4).Imag, [3, 4]);
            % testCase.verifyEqual(ciat.RectangularInterval(1, 2, 3, 4).Imag, [3, 4]);

        % end

        function real(testCase)
            testCase.verifyEqual(real(testCase.rectInt1), ciat.RealInterval(1, 2));
            testCase.verifyEqual(real(testCase.rectInt5), ciat.RealInterval([1, 3], [2, 4]));
            testCase.verifyEqual(real(testCase.rectInt7), ciat.RealInterval([1; 2], [3; 4]));
            testCase.verifyEqual(real(testCase.rectInt9), ciat.RealInterval([1, 2; 3, 4], [5, 6; 7, 8]));
        end

        function imag(testCase)
            testCase.verifyEqual(imag(testCase.rectInt1), ciat.RealInterval(3, 4));
            testCase.verifyEqual(imag(testCase.rectInt5), ciat.RealInterval([3, 6], [4, 8]));
            testCase.verifyEqual(imag(testCase.rectInt7), ciat.RealInterval([5; 6], [7; 8]));
            testCase.verifyEqual(imag(testCase.rectInt9), ciat.RealInterval([9, 10; 11, 12], [13, 14; 15, 16]));
        end

        % function abs(testCase)
        % end

        % function angle(testCase)
        % end

        % function area(testCase)
        % end

        function eq(testCase)
            testCase.verifyEqual(testCase.rectInt1 == testCase.rectInt1, true);
            testCase.verifyEqual(testCase.rectInt1 == testCase.rectInt2, false);

            testCase.verifyEqual(testCase.rectInt5 == testCase.rectInt5, true);
            testCase.verifyEqual(testCase.rectInt5 == testCase.rectInt6, false);

            testCase.verifyEqual(testCase.rectInt7 == testCase.rectInt7, true);
            testCase.verifyEqual(testCase.rectInt7 == testCase.rectInt8, false);

            testCase.verifyEqual(testCase.rectInt9 == testCase.rectInt9, true);
            testCase.verifyEqual(testCase.rectInt9 == testCase.rectInt10, false);
        end

        function ne(testCase)
            testCase.verifyEqual(testCase.rectInt1 ~= testCase.rectInt1, false);
            testCase.verifyEqual(testCase.rectInt1 ~= testCase.rectInt2, true);

            testCase.verifyEqual(testCase.rectInt5 ~= testCase.rectInt5, false);
            testCase.verifyEqual(testCase.rectInt5 ~= testCase.rectInt6, true);

            testCase.verifyEqual(testCase.rectInt7 ~= testCase.rectInt7, false);
            testCase.verifyEqual(testCase.rectInt7 ~= testCase.rectInt8, true);

            testCase.verifyEqual(testCase.rectInt9 ~= testCase.rectInt9, false);
            testCase.verifyEqual(testCase.rectInt9 ~= testCase.rectInt10, true);
        end

        function cast(testCase)
            testCase.verifyEqual(ciat.RectangularInterval(testCase.realInt1), ciat.RectangularInterval(1, 2, 0, 0));
            testCase.verifyEqual(ciat.RectangularInterval(testCase.realInt5), ciat.RectangularInterval([1, 3], [2, 4], [0, 0], [0, 0]));
            testCase.verifyEqual(ciat.RectangularInterval(testCase.realInt7), ciat.RectangularInterval([1; 2], [3; 4], [0; 0], [0; 0]));
            testCase.verifyEqual(ciat.RectangularInterval(testCase.realInt9), ciat.RectangularInterval([1, 2; 3, 4], [5, 6; 7, 8], [0, 0; 0, 0], [0, 0; 0, 0]));
        end
    
        % function sum(testCase)
        % end

        % function uminus(testCase)
        % end

        function plus(testCase)
            testCase.verifyEqual(testCase.rectInt1 + testCase.rectInt2, ciat.RectangularInterval(4, 6, 8, 10));
            testCase.verifyEqual(testCase.rectInt1 + testCase.rectInt3, ciat.RectangularInterval(0, 3, 2, 5));
            testCase.verifyEqual(testCase.rectInt1 + testCase.rectInt4, ciat.RectangularInterval(-1, 1, 1, 3));
            testCase.verifyEqual(testCase.rectInt3 + testCase.rectInt4, ciat.RectangularInterval(-3, 0, -3, 0));

            testCase.verifyEqual(testCase.rectInt1 + 2, ciat.RectangularInterval(3, 4, 3, 4));
            testCase.verifyEqual(-3j + 1 + testCase.rectInt2, ciat.RectangularInterval(4, 5, 2, 3));

            testCase.verifyEqual(testCase.rectInt5 + testCase.rectInt6, ciat.RectangularInterval([4, 7], [6, 12], [8, 16], [10, 20]));
            testCase.verifyEqual(testCase.rectInt5 + -3 + 1j, ciat.RectangularInterval([-2, 0], [-1, 1], [4, 7], [5, 9]));

            testCase.verifyEqual(testCase.rectInt9 + testCase.rectInt10, ciat.RectangularInterval([18, 20; 22, 24], [26, 28; 30, 32], [34, 36; 38, 40], [42, 44; 46, 48]));

            testCase.verifyEqual(testCase.rectInt1 + testCase.rectInt9, ciat.RectangularInterval([2, 3; 4, 5], [7, 8; 9, 10], [12, 13; 14, 15], [17, 18; 19, 20]));
            testCase.verifyEqual(testCase.rectInt9 + testCase.rectInt1, ciat.RectangularInterval([2, 3; 4, 5], [7, 8; 9, 10], [12, 13; 14, 15], [17, 18; 19, 20]));
        end

        function minus(testCase)
            testCase.verifyEqual(testCase.rectInt1 - testCase.rectInt2, ciat.RectangularInterval(-3, -1, -3, -1));
            testCase.verifyEqual(testCase.rectInt1 - testCase.rectInt3, ciat.RectangularInterval(0, 3, 2, 5));
            testCase.verifyEqual(testCase.rectInt1 - testCase.rectInt4, ciat.RectangularInterval(2, 4, 4, 6));
            testCase.verifyEqual(testCase.rectInt3 - testCase.rectInt4, ciat.RectangularInterval(0, 3, 0, 3));
        end
    
        function times(testCase)
            testCase.verifyEqual(testCase.rectInt1 .* testCase.rectInt2, ciat.RectangularInterval(-21, -7, 14, 28));
            testCase.verifyEqual(testCase.rectInt1 .* testCase.rectInt3, ciat.RectangularInterval(-6, 6, -6, 6));
            testCase.verifyEqual(testCase.rectInt1 .* testCase.rectInt4, ciat.RectangularInterval(-1, 7, -12, -4));
            testCase.verifyEqual(testCase.rectInt3 .* testCase.rectInt4, ciat.RectangularInterval(-4, 4, -4, 4));
        end

        % function mtimes(testCase)
        % end
        
        % function union(testCase)
        % end

        % function intersection(testCase)
        % end

        % function union(testCase)
        % end
    end    
end
