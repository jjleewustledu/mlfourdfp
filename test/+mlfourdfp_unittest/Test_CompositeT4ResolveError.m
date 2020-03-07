classdef Test_CompositeT4ResolveError < matlab.unittest.TestCase
	%% TEST_COMPOSITET4RESOLVEERROR 

	%  Usage:  >> results = run(mlfourdfp_unittest.Test_CompositeT4ResolveError)
 	%          >> result  = run(mlfourdfp_unittest.Test_CompositeT4ResolveError, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 29-May-2018 18:13:44 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlfourdfp/test/+mlfourdfp_unittest.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.
 	
	properties
        sessd
        sessf = 'HYGLY28'
 		testObj
        viewer
 	end

	methods (Test)
		function test_afun(this)
 			import mlfourdfp.*;
 			this.assumeEqual(1,1);
 			this.verifyEqual(1,1);
 			this.assertEqual(1,1);
 		end
        function test_meanAndStd(this)
            [m,s] = this.testObj.meanAndStd( ...
                this.testObj.summarizeFrames);
            this.verifyEqual(m, [nan nan nan]);
            this.verifyEqual(s, [nan nan nan]);
        end
        function test_summarizeComposite(this)
            s = this.testObj.summarizeComposite;
            disp(s)
            % verifying only subset (1,1:2)
            this.verifyEqual(s{1}(1,1:2), [], 'RelTol', 1e-10); 
            this.verifyEqual(size(s{1}), [2 2]);
            for m = 1:2
                this.verifyTrue(isnan(s{1}(m,m)));
            end
            pcolor(s{1});
            
            sd = this.sessd;
            sd.epoch = [];
            this.viewer.view(sd.tracerResolvedFinal)
        end
	end

 	methods (TestClassSetup)
		function setupCompositeT4ResolveError(this)
 			import mlfourdfp.*;
 			this.testObj_ = CompositeT4ResolveError;
 		end
	end

 	methods (TestMethodSetup)
		function setupCompositeT4ResolveErrorTest(this)            
 			import mlraichle.*;
            this.sessd = SessionData( ...
                'studyData', StudyData, 'sessionFolder', this.sessf, 'tracer', 'FDG', 'ac', true, 'rnumber', 2); % referenceTracer
 			this.testObj_ = mlfourdfp.CompositeT4ResolveError( ...
                'theImages', {}, ...
                'sessionData', this.sessd, ...
                'logPath', this.sessd.tracerLocation);
            this.viewer = mlfourdfp.Viewer;
 		end
	end

	properties (Access = private)
 		testObj_
 	end

	methods (Access = private)
		function cleanTestMethod(this)
 		end
	end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

