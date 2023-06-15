classdef RealIntervalTest < matlab.unittest.TestCase
    % Tests for the ciat.RealInterval class

    properties
        % Define several real intervals of different sizes on which to perform tests
        
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

        function constructor(testCase)
            % Test different ways to construct a real interval

            % One input
            testCase.verifyEqual(ciat.RealInterval(1), ciat.RealInterval(1,1));
            % testCase.verifyEqual(ciat.RealInterval([3, 4]), ciat.RealInterval([3,3], [4,4])); % AMBIGUITY OF CONSTRUCTOR, NEEDS CHECKING

            % Two inputs
            testCase.verifyEqual(testCase.realInt1.Infimum, 1);
            testCase.verifyEqual(testCase.realInt1.Supremum, 2);

            testCase.verifyEqual([testCase.realInt5.Infimum], [1,3]);
            
        end
        
        function inf(testCase)
            testCase.verifyEqual(inf(testCase.realInt1), 1);
            testCase.verifyEqual(inf(testCase.realInt2), 3);
            testCase.verifyEqual(inf(testCase.realInt3), -1);
            testCase.verifyEqual(inf(testCase.realInt4), -2);
            
            testCase.verifyEqual(inf(testCase.realInt5), [1,3]);
            testCase.verifyEqual(inf(testCase.realInt6), [3,4]);
            
            testCase.verifyEqual(inf(testCase.realInt7), [1;2]);
            testCase.verifyEqual(inf(testCase.realInt8), [2;3]);

            testCase.verifyEqual(inf(testCase.realInt9), [1,2; 3,4]);
            testCase.verifyEqual(inf(testCase.realInt10), [2,4; 6,8]);
        end

        function sup(testCase)
            testCase.verifyEqual(sup(testCase.realInt1), 2);
            testCase.verifyEqual(sup(testCase.realInt2), 4);
            testCase.verifyEqual(sup(testCase.realInt3), 1);
            testCase.verifyEqual(sup(testCase.realInt4), -1);
            
            testCase.verifyEqual(sup(testCase.realInt5), [2,4]);
            testCase.verifyEqual(sup(testCase.realInt6), [5,6]);
            
            testCase.verifyEqual(sup(testCase.realInt7), [3;4]);
            testCase.verifyEqual(sup(testCase.realInt8), [4;5]);

            testCase.verifyEqual(sup(testCase.realInt9), [5,6; 7,8]);
            testCase.verifyEqual(sup(testCase.realInt10), [10,12; 14,16]);
        end

        function mid(testCase)
            testCase.verifyEqual(mid(testCase.realInt1), 1.5);
            testCase.verifyEqual(mid(testCase.realInt2), 3.5);
            testCase.verifyEqual(mid(testCase.realInt3), 0);
            testCase.verifyEqual(mid(testCase.realInt4), -1.5);
            
            testCase.verifyEqual(mid(testCase.realInt5), [1.5,3.5]);
            testCase.verifyEqual(mid(testCase.realInt6), [4,5]);
            
            testCase.verifyEqual(mid(testCase.realInt7), [2;3]);
            testCase.verifyEqual(mid(testCase.realInt8), [3;4]);

            testCase.verifyEqual(mid(testCase.realInt9), [3,4; 5,6]);
            testCase.verifyEqual(mid(testCase.realInt10), [6,8; 10,12]);
        end

        function width(testCase)
            testCase.verifyEqual(width(testCase.realInt1), 1);
            testCase.verifyEqual(width(testCase.realInt2), 1);
            testCase.verifyEqual(width(testCase.realInt3), 2);
            testCase.verifyEqual(width(testCase.realInt4), 1);
            
            testCase.verifyEqual(width(testCase.realInt5), [1,1]);
            testCase.verifyEqual(width(testCase.realInt6), [2,2]);
            
            testCase.verifyEqual(width(testCase.realInt7), [2;2]);
            testCase.verifyEqual(width(testCase.realInt8), [2;2]);

            testCase.verifyEqual(width(testCase.realInt9), [4,4; 4,4]);
            testCase.verifyEqual(width(testCase.realInt10), [8,8; 8,8]);
        end
        
        function eq(testCase)
            testCase.verifyEqual(testCase.realInt1 == testCase.realInt1, true);
            testCase.verifyEqual(testCase.realInt1 == testCase.realInt2, false);

            testCase.verifyEqual(testCase.realInt5 == testCase.realInt5, true);
            testCase.verifyEqual(testCase.realInt5 == testCase.realInt6, false);

            testCase.verifyEqual(testCase.realInt9 == testCase.realInt9, true);
            testCase.verifyEqual(testCase.realInt9 == testCase.realInt10, false);
        end

        function ne(testCase)
            testCase.verifyEqual(testCase.realInt1 ~= testCase.realInt1, false);
            testCase.verifyEqual(testCase.realInt1 ~= testCase.realInt2, true);

            testCase.verifyEqual(testCase.realInt5 ~= testCase.realInt5, false);
            testCase.verifyEqual(testCase.realInt5 ~= testCase.realInt6, true);

            testCase.verifyEqual(testCase.realInt9 ~= testCase.realInt9, false);
            testCase.verifyEqual(testCase.realInt9 ~= testCase.realInt10, true);
        end

        % function sum(testCase)
        % end

        function uminus(testCase)
            testCase.verifyEqual(-testCase.realInt1, ciat.RealInterval(-2,-1));
            testCase.verifyEqual(-testCase.realInt2, ciat.RealInterval(-4,-3));
            testCase.verifyEqual(-testCase.realInt3, ciat.RealInterval(-1,1));
            testCase.verifyEqual(-testCase.realInt4, ciat.RealInterval(1,2));
            
            testCase.verifyEqual(-testCase.realInt5, ciat.RealInterval([-2,-4], [-1,-3]));
            testCase.verifyEqual(-testCase.realInt6, ciat.RealInterval([-5,-6], [-3,-4]));

            testCase.verifyEqual(-testCase.realInt7, ciat.RealInterval([-3;-4], [-1;-2]));
            testCase.verifyEqual(-testCase.realInt8, ciat.RealInterval([-4;-5], [-2;-3]));

            testCase.verifyEqual(-testCase.realInt9, ciat.RealInterval([-5,-6; -7,-8], [-1,-2; -3,-4]));
            testCase.verifyEqual(-testCase.realInt10, ciat.RealInterval([-10,-12; -14,-16], [-2,-4; -6,-8]));
            
        end

        function plus(testCase)
            testCase.verifyEqual(testCase.realInt1 + testCase.realInt2, ciat.RealInterval(4,6));
            testCase.verifyEqual(testCase.realInt1 + testCase.realInt3, ciat.RealInterval(0,3));
            testCase.verifyEqual(testCase.realInt1 + testCase.realInt4, ciat.RealInterval(-1,1));
            testCase.verifyEqual(testCase.realInt3 + testCase.realInt4, ciat.RealInterval(-3,0));
            
            testCase.verifyEqual(testCase.realInt1 + 2, ciat.RealInterval(3,4));
            testCase.verifyEqual(-3 + testCase.realInt2, ciat.RealInterval(0,1));


            testCase.verifyEqual(testCase.realInt5 + testCase.realInt6, ciat.RealInterval([4,7], [7,10]));
            testCase.verifyEqual(testCase.realInt5 + 3, ciat.RealInterval([4,6], [5,7]));


            testCase.verifyEqual(testCase.realInt9 + testCase.realInt10, ciat.RealInterval([3,6; 9,12], [15,18; 21,24]));


            testCase.verifyEqual(testCase.realInt1 + testCase.realInt9, ciat.RealInterval([2,3; 4,5], [7,8; 9,10]));
            testCase.verifyEqual(testCase.realInt9 + testCase.realInt1, ciat.RealInterval([2,3; 4,5], [7,8; 9,10]));

        end

        function minus(testCase)
            testCase.verifyEqual(testCase.realInt1 - testCase.realInt2, ciat.RealInterval(-3,-1));
            testCase.verifyEqual(testCase.realInt1 - testCase.realInt3, ciat.RealInterval(0,3));
            testCase.verifyEqual(testCase.realInt1 - testCase.realInt4, ciat.RealInterval(2,4));
            testCase.verifyEqual(testCase.realInt3 - testCase.realInt4, ciat.RealInterval(0,3));
        end

        function times(testCase)
            testCase.verifyEqual(testCase.realInt1 .* testCase.realInt2, ciat.RealInterval(3,8));
            testCase.verifyEqual(testCase.realInt1 .* testCase.realInt3, ciat.RealInterval(-2,2));
            testCase.verifyEqual(testCase.realInt1 .* testCase.realInt4, ciat.RealInterval(-4,-1));
            testCase.verifyEqual(testCase.realInt3 .* testCase.realInt4, ciat.RealInterval(-2,2));
        end

        % function mtimes(testCase)
        % end

        % function rdivide(testCase)
        % end

        % function mrdivide(testCase)
        % end

        % function recip(testCase)
        % end

        function division(testCase)
            testCase.verifyEqual(testCase.realInt1 / testCase.realInt2, ciat.RealInterval(1/4, 2/3));
            testCase.verifyEqual(testCase.realInt3 / testCase.realInt4, ciat.RealInterval(-1, 1));
        end

        % function abs(testCase)
        % end

        % function exp(testCase)
        % end

        % function log(testCase)
        % end

        % function log10(testCase)
        % end

        % function sqrt(testCase)
        % end

        % function sin(testCase)
        % end

        % function cos(testCase)
        % end

        % function union(testCase)
        % end

        % function intersection(testCase)
        % end

    end    
end