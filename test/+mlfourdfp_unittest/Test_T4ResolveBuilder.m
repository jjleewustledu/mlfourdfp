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
 	%% It was developed on Matlab 9.0.0.307022 (R2016a) Prerelease for MACI64.  Copyright 2017 John Joowon Lee.
 	

	properties
        fourdfpv
        sessf = 'HYGLY09'
        ipr
        pwd0
 		studyd
        sessd
 		testObj
        tracer = 'FDG'
        
        view = false
        viewer
        quick = false
    end
    
    methods (Static)
        function constructTestData
        end
    end
    
	methods (Test)
        function test_setup(this)
            this.verifyClass(this.testObj, 'mlfourdfp.T4ResolveBuilder');
            this.verifyClass(this.testObj.studyData, 'mlraichle.SynthStudyData');
            this.verifyClass(this.testObj.sessionData, 'mlraichle.SynthSessionData');
        end
        function test_finished(this)
            touchfile = fullfile(this.testObj.getLogPath, '.test_mlfourdfp.T4ResolveBuilder_isfinished.touch');
            if (lexist(touchfile)); delete(touchfile); end 
            this.verifyEqual(this.testObj.finished.markerFilename, touchfile);           
            this.testObj.finished.markAsFinished;
            this.verifyTrue(lexist(this.testObj.finished.markerFilename, 'file'));
            this.verifyTrue(this.testObj.finished.isfinished);
            delete(touchfile); 
        end
        function test_reconstituteFramesAC2(this)
            mlbash('freeview E1/fdgv1e1r2_op_fdgv1e1to4r1_frame4.4dfp.img E2/fdgv1e2r2_op_fdgv1e1to4r1_frame4.4dfp.img E3/fdgv1e3r2_op_fdgv1e1to4r1_frame4.4dfp.img E4/fdgv1e4r1.4dfp.img E1to4/fdgv1e1to4r2_op_fdgv1e1to4r1_frame4_sumt.4dfp.img');
        end
        function test_resolve(this)
            tic
            this.testObj.NRevisions = 2;
            this.testObj = this.testObj.resolve('dest', 'testv1r1', 'source', 'testNativeCropped');
            [e,c] = this.fourdfpv.etaAndCurves('testv1r2_frame3_Test_T4ResolveBuilder', 'testv1r2_frame5_Test_T4ResolveBuilder');
            this.verifyEqual(e, 0.1631, 'RelTol', 0.05);
            %this.verifyEqual(c, [130 133 21 89 161 66], 'RelTol', 0.1);
            this.verifyEqual(this.testObj.product.niftid.entropy, 3.4182, 'RelTol', 0.05);
            this.verifyTrue(this.testObj.isfinished(this.testObj.sessionData));            
            fprintf('test_resolve:  ');
            toc
            if (this.view)
                this.testObj.product.niftid.fslview;
            end
        end
        function test_revise(this)
            if (this.quick)
                return
            end
            this.ipr.resolved = 'testNativeCropped';
            [~,this.testObj] = this.testObj.revise(this.ipr);
            this.verifyEqual(this.testObj.sessionData.rnumber, 2);
        end
        function test_imageReg(this)
            if (this.quick)
                return
            end
            
            this.testObj.imageReg(this.ipr);
            
            in3 = 'testv1r1_frame3';
            in5 = 'testv1r1_frame5';
            out = sprintf('%s_on_%s', in5, in3);
            this.fourdfpv.t4_mul([in3 '_to_' in5 '_t4'], [ in5 '_to_' in3 '_t4']);
            movefile([in3 '_to_' in3 '_t4'], [in3 '_to_' in3 '_t4.txt']);
            t4 = readtable([in3 '_to_' in3 '_t4'], ...
                'ReadVariableNames', 0, 'Delimiter', ' ', 'MultipleDelimsAsOne', 1, 'HeaderLines', 4, 'Format', '%f %f %f %f');
            this.verifyEqual(t4{1,1}, 1, 'AbsTol', 4e-3)
            this.verifyEqual(t4{1,2}, 0, 'AbsTol', 7e-2)
            this.verifyEqual(t4{1,3}, 0, 'AbsTol', 9e-3)
            this.verifyEqual(t4{1,4}, 0, 'AbsTol', 0.6)
            
            this.fourdfpv.t4img_4dfp(this.fourdfpv.filenameT4(in5, in3), in5, 'out', out, 'options', ['-O' in3]);
            if (this.view)
                mlbash(sprintf('fslview %s.4dfp.img %s.4dfp.img', in3, out));
            end
        end
        function test_t4ResolveAndPaste(this)
            if (this.quick)
                return
            end
            
            this.testObj.resolveLog = this.testObj.loggerFilename('', 'func', 'Test_T4ResolveBuilder_test_t4ResolveAndPaste');
            [ipr_,imgFns_] = this.testObj.resolveAndPaste(this.ipr);
            this.verifyEqual(ipr_.resolved, 'testv1r1_Test_T4ResolveBuilder');
            this.verifyEqual(imgFns_, 'testv1r1_frame3 testv1r1_frame5');
            if (this.view)
                mlbash('fslview testv1r1_frame3_Test_T4ResolveBuilder.4dfp.img testv1r1_frame5_Test_T4ResolveBuilder.4dfp.img');
            end
        end
        function test_teardownT4s(this)
            %this.testObj.teardownT4;
        end
        function test_teardownLogs(this)
            %this.testObj.teardownLogs;
        end
        function test_teardownRevision(this)
            %this.testObj.teardownRevision;
        end
        function test_teardownResolve(this)
            %this.testObj.teardownResolve;
        end
        function test_resolveTag(this)
        end
        function test_productOnFilesystem(this)
            pwd1 = pushd(fullfile(this.sessd.tracerLocation));
            
            fp = 'fdgv1r2_op_fdgv1e1to4r1_frame4';
            this.fourdfpv.imgblur_4dfp(fp, 11); % overwrite previous testing
            fdg = mlfourd.NumericalNIfTId.load([fp '.4dfp.hdr']);
            fdg = fdg.timeAveraged;
            fdg.filesuffix = '.4dfp.hdr';
            fdg.save; % overwrite previous testing
            this.viewer.view({[fp '_b110.4dfp.hdr'] fdg.filename});
            
            popd(pwd1);
        end
        
        %% UTILITY
        
        function test_fileprefixBlurred(this)
        end
        function test_fileprefixGaussed(this)
        end
        function test_fileprefixSumt(this)
        end
        function test_lazyExtractImage(this)
            this.ipr.currentIndex = 1;
            fqfp = this.testObj.imageComposite.lazyExtractImage(this.ipr);
            this.verifyEqual(mybasename(fqfp), 'testv1r1_frame1');
            if (this.view)
                mlbash('fslview testNativeCropped.4dfp.img -b 0,100 testv1r1_frame1.4dfp.img -b 0,100');
            end
        end
        function test_lazyStageImages(this)
            this.ipr.currentIndex = nan;
            
            fqfps = this.testObj.lazyStageImages(this.ipr);
            this.verifyEqual(length(fqfps), 2);
            this.verifyEqual(mybasename(fqfps{1}), 'testv1r1_frame3');
            this.verifyEqual(mybasename(fqfps{2}), 'testv1r1_frame5');
            if (this.view)
                mlbash('fslview testNativeCropped.4dfp.img -b 0,2000 testv1r1_frame3.4dfp.img -b 0,2000 testv1r1_frame5.4dfp.img -b 0,2000');
            end
        end
        function test_lazyBlurImage(this)
            if (this.quick)
                return
            end
            
            this.ipr.currentIndex = 1;
            fqfp = this.testObj.lazyBlurImage(this.ipr, 5.5);
            this.verifyEqual(mybasename(fqfp), 'testv1r1_frame1_b55');
            if (this.view)
                mlbash('fslview testNativeCropped.4dfp.img -b 0,100 testv1r1_frame1_b55.4dfp.img -b 0,100');
            end
        end
        function test_lazyBlurImages(this)
            if (this.quick)
                return
            end
            
            this.ipr.currentIndex = nan;
            
            fqfps = this.testObj.lazyBlurImages(this.ipr);
            this.verifyEqual(length(fqfps), 2);
            this.verifyEqual(mybasename(fqfps{1}), 'testv1r1_frame3_b55');
            this.verifyEqual(mybasename(fqfps{2}), 'testv1r1_frame5_b55');
            if (this.view)
                mlbash('fslview testNativeCropped.4dfp.img -b 0,2000 testv1r1_frame3_b55.4dfp.img -b 0,2000 testv1r1_frame5_b55.4dfp.img -b 0,2000');
            end
        end
        function test_lazyMaskForImages(this)
            if (this.quick)
                return
            end
            
            tic
            fqfp = this.testObj.lazyMaskForImages('maskForImages', 'testv1r1_frame3', 'testv1r1_frame5', 3, 5);
            this.verifyEqual(mybasename(fqfp), 'maskForImages_3_5');
            if (this.view)
                mlbash('fslview testv1r1_frame3.4dfp.img -b 0,2000 testv1r1_frame5.4dfp.img -b 0,2000 maskForImages_3_5.4dfp.img');
            end
            fprintf('test_lazyMaskForImages:  ');
            toc
        end        
        function test_frames(this)
            this.verifyEqual(sum(this.testObj.indicesLogical), 2);
        end
	end

 	methods (TestClassSetup)
		function setupT4ResolveBuilder(this)
            import mlraichle.*;
            this.studyd = StudyData;
            this.sessd = SynthSessionData( ...
                'studyData', this.studyd, ...
                'sessionFolder', this.sessf, ...
                'tracer', this.tracer);    
            
            this.pwd0 = pushd(this.sessd.tracerLocation);
            
            this.testObj_ = mlfourdfp.T4ResolveBuilder( ... 
                'sessionData', this.sessd, ...  
                'NRevisions', 1, ...
                'indicesLogical', this.indicesLogical, ... 
                'indexOfReference', 3);
            this.ipr = struct( ...
                'source', 'testNativeCropped', ...
                'dest', 'testv1r1', ...
                'maskForImages', 'maskForImages', ...
                'sourceBlur', [], ...
                'destBlur', [], ...
                'keepForensics', true);
            this.ipr.sourceBlur = this.blurArray;
            this.ipr.destBlur = this.blurArray;
            this.ensureTestImages;
            this.fourdfpv = mlfourdfp.FourdfpVisitor;
            this.viewer = mlfouredfp.Viewer;
            setenv('DEBUG', ''); % cf. dbbash
 			this.addTeardown(@this.teardownClass);
 		end
	end

 	methods (TestMethodSetup)
		function setupT4ResolveBuilderTest(this)
 			this.testObj = this.testObj_;
 			this.addTeardown(@this.teardownMethod);
 		end
    end
    
    %% PRIVATE

	properties (Access = private)
 		testObj_
 	end

	methods (Access = private)
        function b = blurArray(~)
            b = zeros(1, 5);
            b(3) = 5.5;
            b(5) = 5.5;
        end
        function ensureTestImages(this)
            if (~this.fourdfpv.lexist_4dfp('testNativeCropped'))
                this.fourdfpv.copy_4dfp('TEST_V1-LM-00-OP', 'testNativeCropped');
            end
            if (~this.fourdfpv.lexist_4dfp('testv1r1'))
                this.fourdfpv.copy_4dfp('testNativeCropped', 'testv1r1');
            end
        end
        function f = indicesLogical(~)
            f = false(1, 5);
            f(3) = true;           
            f(5) = true;
        end
		function teardownMethod(~)
 		end
		function teardownClass(this)
            if (~this.quick)
                this.testObj_.keepForensics = false;
                this.testObj_.teardownResolve(this.ipr);
                %delete('*_frame*');
                %delete('*Test_T4ResolveBuilder.4dfp.*');
                delete('*.nii');
                delete('*.nii.gz');
                if (isdir('Log')); rmdir('Log', 's'); end
                if (isdir('T4'));  rmdir('T4',  's'); end
            end
            delete('*.mat0');
            delete('*.sub');
            
            popd(this.pwd0);
 		end
	end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

