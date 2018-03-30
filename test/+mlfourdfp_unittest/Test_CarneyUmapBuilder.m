classdef Test_CarneyUmapBuilder < matlab.unittest.TestCase
	%% TEST_CARNEYUMAPBUILDER 

	%  Usage:  >> results = run(mlfourdfp_unittest.Test_CarneyUmapBuilder)
 	%          >> result  = run(mlfourdfp_unittest.Test_CarneyUmapBuilder, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 14-Dec-2016 17:51:20
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/test/+mlfourdfp_unittest.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.
 	

	properties
        groundTruth
        hyglyNN = 'HYGLY00'
        pwd0
 		registry
        sessd
        studyd
 		testObj
 	end

	methods (Test)
		function test_ctor(this)
            this.verifyClass(this.testObj, 'mlfourdfp.CarneyUmapBuilder');
 		end
		function test_buildUmap(this)
            [~,umap] = this.testObj_.buildUmap;
            umap = mlfourd.NIfTId.load([umap '.4dfp.ifh']);
            umap.fslview;
            this.verifyTrue(this.testObj_.isfinished);
        end
        function test_testGroundTruth(this)
            pwd1 = pushd(fileparts(this.groundTruth));
            umap = this.testObj.testGroundTruth(this.groundTruth);
            umap = mlfourdfp.NIfTId.load([umap '.4dfp.ifh']);
            umap.fslview;
            popd(pwd1);
        end
	end

 	methods (TestClassSetup)
		function setupCarneyUmapBuilder(this)
            [~,h] = mlbash('hostname');
            assert(lstrfind(h, 'william'));
            this.groundTruth = fullfile(getenv('PPG'), 'test', '1011024_T80_v15', 'ct_1011024_T80_v15');
            this.studyd = mlraichle.SynthStudyData;
            this.sessd = mlraichle.SynthSessionData( ...
                'studyData', this.studyd, ...
                'sessionPath', fullfile(this.studyd.subjectsDir, this.hyglyNN, ''), ...
                'tracer', 'FDG');
 			this.testObj_ = mlfourdfp.CarneyUmapBuilder( ...
                'sessionData', this.sessd); 
            %, ...
            %    'theImages', this.groundTruth, ...
            %    'indicesLogical', true);
 		end
	end

 	methods (TestMethodSetup)
		function setupCarneyUmapBuilderTest(this)
 			this.testObj = this.testObj_;
            this.pwd0 = pushd(this.sessd.sessionLocation);
 			this.addTeardown(@this.cleanFiles);
 		end
	end

	properties (Access = private)
 		testObj_
 	end

	methods (Access = private)
		function cleanFiles(this)
            if (lexist(this.testObj_.finished.finishedMarkerFilename, 'file'))
                delete(this.testObj_.finished.finishedMarkerFilename); end
            popd(this.pwd0);
 		end
	end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

