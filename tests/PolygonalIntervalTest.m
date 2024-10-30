classdef PolygonalIntervalTest < matlab.unittest.TestCase

    properties
        % Define several random polygonal intervals of different sizes on which to perform tests

        % Set uniform random bounds
        Size    = [1,1];
        RndPoints
    end
    
    methods
        function points = get.RndPoints(testCase)
            M = testCase.Size(1);
            N = testCase.Size(2);
            points = cell(M,N);
            for m = 1:M
            for n = 1:N
                L = randi(10)+2;
                points{m,n} = complex(unifrnd(-10,10,L,1), ...
                                      unifrnd(-10,10,L,1));
            end
            end
        end
    end

    methods(TestClassSetup)
        % Shared setup for the entire test class
        
    end
    
    methods(TestMethodSetup)
        % Setup for each test
    end
    
    methods(Test)
        % Test methods

        function eq(testCase)

            % 1x1 case
            testCase.Size = [1,1];
            p1 = ciat.PolygonalInterval(testCase.RndPoints);
            p2 = ciat.PolygonalInterval(testCase.RndPoints);
            testCase.verifyEqual(p1 == p1, true);
            testCase.verifyEqual(p2 == p2, true);
            testCase.verifyEqual(p1 == p2, false);
            testCase.verifyEqual(p2 == p1, false);

            % 1x2 case
            testCase.Size = [1,2];
            p1 = ciat.PolygonalInterval(testCase.RndPoints);
            p2 = ciat.PolygonalInterval(testCase.RndPoints);
            testCase.verifyEqual(all(p1 == p1), true);
            testCase.verifyEqual(all(p2 == p2), true);
            testCase.verifyEqual(all(p1 == p2), false);
            testCase.verifyEqual(all(p2 == p1), false);

            % 1x2 case
            testCase.Size = [2,1];
            p1 = ciat.PolygonalInterval(testCase.RndPoints);
            p2 = ciat.PolygonalInterval(testCase.RndPoints);
            testCase.verifyEqual(all(p1 == p1), true);
            testCase.verifyEqual(all(p2 == p2), true);
            testCase.verifyEqual(all(p1 == p2), false);
            testCase.verifyEqual(all(p2 == p1), false);

            % 1x2 case
            testCase.Size = [2,2];
            p1 = ciat.PolygonalInterval(testCase.RndPoints);
            p2 = ciat.PolygonalInterval(testCase.RndPoints);
            testCase.verifyEqual(all(p1 == p1,'all'), true);
            testCase.verifyEqual(all(p2 == p2,'all'), true);
            testCase.verifyEqual(all(p1 == p2,'all'), false);
            testCase.verifyEqual(all(p2 == p1,'all'), false);
        end
        
        function ne(testCase)

            % 1x1 case
            testCase.Size = [1,1];
            p1 = ciat.PolygonalInterval(testCase.RndPoints);
            p2 = ciat.PolygonalInterval(testCase.RndPoints);
            testCase.verifyEqual(p1 ~= p1, false);
            testCase.verifyEqual(p2 ~= p2, false);
            testCase.verifyEqual(p1 ~= p2, true);
            testCase.verifyEqual(p2 ~= p1, true);

            % 1x2 case
            testCase.Size = [1,2];
            p1 = ciat.PolygonalInterval(testCase.RndPoints);
            p2 = ciat.PolygonalInterval(testCase.RndPoints);
            testCase.verifyEqual(all(p1 ~= p1), false);
            testCase.verifyEqual(all(p2 ~= p2), false);
            testCase.verifyEqual(all(p1 ~= p2), true);
            testCase.verifyEqual(all(p2 ~= p1), true);

            % 1x2 case
            testCase.Size = [2,1];
            p1 = ciat.PolygonalInterval(testCase.RndPoints);
            p2 = ciat.PolygonalInterval(testCase.RndPoints);
            testCase.verifyEqual(all(p1 ~= p1), false);
            testCase.verifyEqual(all(p2 ~= p2), false);
            testCase.verifyEqual(all(p1 ~= p2), true);
            testCase.verifyEqual(all(p2 ~= p1), true);

            % 1x2 case
            testCase.Size = [2,2];
            p1 = ciat.PolygonalInterval(testCase.RndPoints);
            p2 = ciat.PolygonalInterval(testCase.RndPoints);
            testCase.verifyEqual(all(p1 ~= p1,'all'), false);
            testCase.verifyEqual(all(p2 ~= p2,'all'), false);
            testCase.verifyEqual(all(p1 ~= p2,'all'), true);
            testCase.verifyEqual(all(p2 ~= p1,'all'), true);
        end

        function cast(testCase)
            testCase.Size = [10,10];
            gI = ciat.PolygonalInterval(testCase.RndPoints);
            rI = ciat.RectangularInterval(gI);
            pI = ciat.PolarInterval(gI);
            cI = ciat.CircularInterval(gI);

            rIin = false(10);
            pIin = false(10);
            cIin = false(10);
            for m=1:10
            for n=1:10
                rIin(m,n) = all(rI(m,n).isin(gI(m,n).Points));
                pIin(m,n) = all(pI(m,n).isin(gI(m,n).Points));
                cIin(m,n) = all(cI(m,n).isin(gI(m,n).Points));
            end
            end
            testCase.verifyEqual(all(rIin,'all'),true);
            testCase.verifyEqual(all(pIin,'all'),true);
            testCase.verifyEqual(all(cIin,'all'),true);
        end

        function real(testCase)
            testCase.Size = [10,10];
            gI = ciat.PolygonalInterval(testCase.RndPoints);
            rI = ciat.RectangularInterval(gI);
            testCase.verifyEqual(all(gI.Real == rI.Real,'all'), true);
        end

        function imag(testCase)
            testCase.Size = [10,10];
            gI = ciat.PolygonalInterval(testCase.RndPoints);
            rI = ciat.RectangularInterval(gI);
            testCase.verifyEqual(all(gI.Imag == rI.Imag,'all'), true);
        end

        function abs(testCase)
            testCase.Size = [10,10];
            gI = ciat.PolygonalInterval(testCase.RndPoints);
            pI = ciat.PolarInterval(gI);
            testCase.verifyEqual(all(gI.Abs == pI.Abs,'all'), true);
        end

        function angle(testCase)
            testCase.Size = [10,10];
            gI = ciat.PolygonalInterval(testCase.RndPoints);
            pI = ciat.PolarInterval(gI);
            testCase.verifyEqual(all(gI.Angle == pI.Angle,'all'), true);
        end
        
        function plus(testCase)

            % Add two 1x1 intervals
            testCase.Size = [1,1];
            gI1 = ciat.PolygonalInterval(testCase.RndPoints);
            gI2 = ciat.PolygonalInterval(testCase.RndPoints);
            rI1 = ciat.RectangularInterval(gI1);
            rI2 = ciat.RectangularInterval(gI2);
            gIs = gI1 + gI2;
            rIs = rI1 + rI2;
            testCase.verifyEqual(all(rIs.isin(gIs.Points)), true);

            % Add two 2x1 intervals
            testCase.Size = [2,1];
            gI1 = ciat.PolygonalInterval(testCase.RndPoints);
            gI2 = ciat.PolygonalInterval(testCase.RndPoints);
            rI1 = ciat.RectangularInterval(gI1);
            rI2 = ciat.RectangularInterval(gI2);
            gIs = gI1 + gI2;
            rIs = rI1 + rI2;
            rIin = false(2,1);
            for m=1:2
                rIin(m) = all(rIs(m).isin(gIs(m).Points));
            end
            testCase.verifyEqual(all(rIin), true);

            % Add 2x1 and a 1x2 interval
            testCase.Size = [2,1];
            gI1 = ciat.PolygonalInterval(testCase.RndPoints);
            testCase.Size = [1,2];
            gI2 = ciat.PolygonalInterval(testCase.RndPoints);
            rI1 = ciat.RectangularInterval(gI1);
            rI2 = ciat.RectangularInterval(gI2);
            gIs = gI1 + gI2;
            rIs = rI1 + rI2;
            rIin = false(2,2);
            for m=1:2
            for n=1:2
                rIin(m,n) = all(rIs(m,n).isin(gIs(m,n).Points));
            end
            end
            testCase.verifyEqual(all(rIin,'all'), true);

            % Add two 2x2 intervals
            testCase.Size = [2,2];
            gI1 = ciat.PolygonalInterval(testCase.RndPoints);
            gI2 = ciat.PolygonalInterval(testCase.RndPoints);
            rI1 = ciat.RectangularInterval(gI1);
            rI2 = ciat.RectangularInterval(gI2);
            gIs = gI1 + gI2;
            rIs = rI1 + rI2;
            rIin = false(2,2);
            for m=1:2
            for n=1:2
                rIin(m,n) = all(rIs(m,n).isin(gIs(m,n).Points));
            end
            end
            testCase.verifyEqual(all(rIin,'all'), true);
        end

    end
    
end