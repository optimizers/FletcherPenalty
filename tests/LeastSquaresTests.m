classdef LeastSquaresTests < matlab.unittest.TestCase
    
   properties
      A
      b
      c
      Aprod
      Atprod
      delta = 0;
   end
    
   methods(TestMethodSetup)
      function setupProblem(testCase)
          m = 50;
          n = 100;
          testCase.A = sprand(n,m,0.2);
          testCase.b = rand(n,1);
          testCase.c = rand(m,1);
          testCase.Aprod = @(v) testCase.A*v;
          testCase.Atprod = @(v) testCase.A'*v;
      end
   end
    
   methods(Test)
       
      function testSNE(testCase)
         % Test without regularization
         options.regularized = false;
         lsq = least_squares.lssne(testCase.A, options);
         LeastSquaresTests.solveLSQ(testCase, lsq, testCase.delta);
         
         % Test with regularization
         testCase.delta = 1e-3;
         options.regularized = true;
         lsq = least_squares.lssne(testCase.A, options);
         lsq.delta = testCase.delta;
         LeastSquaresTests.solveLSQ(testCase, lsq, testCase.delta);
      end
       
      function testLDL(testCase)
         % Test without regularization
         options.regularized = false;
         
         lsq = least_squares.lsldl(testCase.A, options);
         LeastSquaresTests.solveLSQ(testCase, lsq, 0);
         
         % Test with regularization
         testCase.delta = 1e-3;
         options.regularized = true;
         lsq = least_squares.lsldl(testCase.A, options);
         lsq.delta = testCase.delta;
         LeastSquaresTests.solveLSQ(testCase, lsq, testCase.delta);
      end
      
      function testMINRES(testCase)
         % Test without regularization
         options.regularized = false;
         options.tol = 1e-14;
         
         [n,m] = size(testCase.A);
         lsq = least_squares.lsminres(m, n, options);
         LeastSquaresTests.solveLSQ(testCase, lsq, 0);
         
         % Test with regularization
         testCase.delta = 1e-1;
         options.regularized = true;
         lsq = least_squares.lsminres(m, n, options);
         lsq.delta = testCase.delta;
         LeastSquaresTests.solveLSQ(testCase, lsq, testCase.delta);
      end
      
      function testLNLQ(testCase)
         if ~exist('lnlq','file')
            fprintf('LNLQ not found on path. Skipping test.\n');
            return;
         end
          
         % Test without regularization
         options.regularized = false;
         options.tol = 1e-14;
         
         [n,m] = size(testCase.A);
         lsq = least_squares.lsminres(m, n, options);
         LeastSquaresTests.solveLSQ(testCase, lsq, 0);
         
         % Test with regularization
         testCase.delta = 1e-1;
         options.regularized = true;
         lsq = least_squares.lsminres(m, n, options);
         lsq.delta = testCase.delta;
         LeastSquaresTests.solveLSQ(testCase, lsq, testCase.delta);
      end
      
   end
   
   methods(Static)
       function solveLSQ(testCase, lsq, delta)
          lsq = preprocess(lsq, testCase.A, ...
             testCase.Aprod, testCase.Atprod, [], struct());
           
          [y,x] = lsq.lsq3(testCase.b, testCase.c);
         
          [n,m] = size(testCase.A);
          K = [speye(n) testCase.A; testCase.A' -delta^2*speye(m,m)];
          r = K*[x;y] - [testCase.b;testCase.c];
         
          testCase.verifyEqual(norm(r,'inf'), 0, 'AbsTol', 1e-12);
       end
   end
end