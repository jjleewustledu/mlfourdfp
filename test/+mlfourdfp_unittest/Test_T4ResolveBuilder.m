classdef Test_T4ResolveBuilder < matlab.unittest.TestCase
	%% TEST_T4RESOLVEBUILDER 

	%  Usage:  >> results = run(mlfourdfp_unittest.Test_T4ResolveBuilder)
 	%          >> result  = run(mlfourdfp_unittest.Test_T4ResolveBuilder, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 10-Mar-2016 21:30:00
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/test/+mlfourdfp_unittest.
 	%% It was developed on Matlab 9.0.0.307022 (R2016a) Prerelease for MACI64.
 	

	properties
 		studyd
        sessd
 		testObj
        t4ri
 	end

	methods (Test)
		function test_t4ResolveSubject(this)
            fprintf('Test_T4ResolveBuilder.test_t4ResolveSubject:  running t4ResolveSubject which may requires hours of processing time..........');
            this.testObj  = this.testObj.t4ResolveSubject;
            this.verifyTrue(~isempty(this.testObj.product));
            
            if (this.view)
                product.fdg.view;
            end
        end
        function test_t4ResolvePET(this)
%             this.touch_4dfp('ho2v1');
%             this.touch_4dfp('ho2v1_sumt');
%             this.touch_4dfp('ho2v1_sumt_g11');
%             mlbash('touch ho2v1_to_HYGLY09_mpr_t4');
%             this.touch_4dfp('ho2v1_on_HYGLY09_mpr');
%             mlbash('touch ho2v1_to_TRIO_Y_NDC_t4');
%             this.touch_4dfp('ho2v1_sumt_b55');
%             for f = 1:48
%                 frameFps{f} = sprintf('ho2v1_frame%i', f);
%                 this.touch_4dfp( frameFps{f});
%                 this.touch_4dfp([frameFps{f} '_b55']);
%             end
%             this.touch_4dfp('HYGLY09HO2_v1_NAC_frame4');
%             this.touch_4dfp('ho2v1_sumt_mskt');
%             for f = 4:48
%                 mlbash(sprintf('touch ho2v1_frame%i_to_resolved_t4', f));
%                 this.touch_4dfp(sprintf('ho2v1_frame%i', f));
%             end
%             this.touch_4dfp('ho2v1_frames4to31_resolved');
%             mlbash('touch resolved.mat0');
%             mlbash('touch resolved.sub');
%             this.touch_4dfp('ho2v1_sumt_mskt');
%             this.touch_4dfp('ho2v1_sumt_mskt_b55');
%             this.touch_4dfp('ho2v1_b55_mskt');
%             this.touch_4dfp('ho2v1_b55');
            
            this.testObj.t4ResolvePET;
        end
        function test_readFrameEnd(this)
            this.verifyEqual(this.testObj.readFrameEnd('HYGLY09FDG_v1_NAC'), 31);
        end
        function test_t4ResolveIteration(this)
            this.testObj.t4ResolveIteration( ...
                'HYGLY09FDG_v1_NAC', 'fdgv1', 'HYGLY09_mpr');
        end
        function test_msktgenInitial(this)
            this.touch_4dfp('fdgv1');
            this.touch_4dfp('fdgv1_sumt');
            this.touch_4dfp('fdgv1_sumt_g11');
            mlbash('touch fdgv1_to_HYGLY09_mpr_t4');
            this.touch_4dfp('fdgv1_on_HYGLY09_mpr');
            mlbash('touch fdgv1_to_TRIO_Y_NDC_t4');
            this.touch_4dfp('fdgv1_sumt_b55');
            this.testObj.msktgenInitial(this.t4ri);
        end
        function test_extractFrames(this)
            fps = this.testObj.extractFrames(this.t4ri);
            for f = 1:length(fps)
                fprintf('%s\n', fps{f});
            end
        end
        function test_frameReg(this)
            for f = 1:31
                frameFps{f} = sprintf('fdgv1_frame%i', f);
                this.touch_4dfp( frameFps{f});
                this.touch_4dfp([frameFps{f} '_b55']);
            end
            this.touch_4dfp('HYGLY09FDG_v1_NAC_frame4');
            this.touch_4dfp('fdgv1_sumt_mskt');
            this.testObj.frameReg(this.t4ri, frameFps);
        end
        function test_t4ResolveAndPaste(this)
            for f = 4:31
                mlbash(sprintf('touch fdgv1_frame%i_to_resolved_t4', f));
                this.touch_4dfp(sprintf('fdgv1_frame%i', f));
            end
            this.touch_4dfp('fdgv1_frames4to31_resolved');
            mlbash('touch resolved.mat0');
            mlbash('touch resolved.sub');
            this.testObj.t4ResolveAndPaste(this.t4ri);
        end
        function test_t4ResolveTeardown(this)
            this.testObj.t4ResolveTeardown(this.t4ri);
        end
	end

 	methods (TestClassSetup)
		function setupT4ResolveBuilder(this)
 			import mlfourdfp.*;
            
            this.studyd = mlpipeline.StudyDataSingletons.instance('test_raichle');
            this.sessd = mlraichle.SessionData( ...
                'studyData', this.studyd, 'sessionPath', fullfile(this.studyd.subjectsDir, 'HYGLY09', ''));            
            this.testObj_ = T4ResolveBuilder('sessionData', this.sessd);
            this.t4ri = struct( ...
                'fdfp0', 'HYGLY09FDG_v1', ...
                'fdfp1', 'fdgv1', ...
                'mprage', 'HYGLY09_mpr', ...
                'frame0', 4, ...
                'frameF', 31, ...
                'crop', 0.5, ...
                'atlas', 'TRIO_Y_NDC', ...
                'blur', 5.5);
            
            %setenv('REFDIR',  '/Volumes/InnominateHD3/Local/test/raichle/atlas');
            %setenv('RELEASE', '/Volumes/InnominateHD3/Local/test/raichle/lin64-tools');
            setenv('DEBUG', '0');
 		end
	end

 	methods (TestMethodSetup)
		function setupT4ResolveBuilderTest(this)
 			this.testObj = this.testObj_;            
            cd(fullfile(this.sessd.sessionPath, 'FDG_v1', ''));
 			this.addTeardown(@this.cleanFiles);
 		end
	end

	properties (Access = private)
 		testObj_
 	end

	methods (Access = private)
        function touch_4dfp(this, fp)
            mlbash(sprintf('touch %s.4dfp.hdr', fp));
            mlbash(sprintf('touch %s.4dfp.ifh', fp));
            mlbash(sprintf('touch %s.4dfp.img', fp));
            mlbash(sprintf('touch %s.4dfp.img.rec', fp));
        end
		function cleanFiles(this)
 		end
	end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

