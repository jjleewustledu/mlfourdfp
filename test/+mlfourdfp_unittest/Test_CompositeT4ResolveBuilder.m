classdef Test_CompositeT4ResolveBuilder < matlab.unittest.TestCase
	%% TEST_COMPOSITET4RESOLVEBUILDER 

	%  Usage:  >> results = run(mlfourdfp_unittest.Test_CompositeT4ResolveBuilder)
 	%          >> result  = run(mlfourdfp_unittest.Test_CompositeT4ResolveBuilder, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 21-Jan-2017 13:47:31
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/test/+mlfourdfp_unittest.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.  Copyright 2017 John Joowon Lee.
 	

	properties
        fvisitor
        hyglyNN = 'HYGLY09'
        ipr
        pwd0
 		studyd
        resolveTag = 'Test_CompositeT4ResolveBuilder'
        sessd
 		testObj
        theImages = {'testNativeCropped1' 'testNativeCropped2' 'testNativeCropped3' 'testNativeCropped4' 'testNativeCropped5'}
        
        view = false
        quick = false
    end

	methods (Test)
        function test_setup(this)
            this.verifyClass(this.testObj, 'mlfourdfp.CompositeT4ResolveBuilder');
            this.verifyClass(this.testObj.studyData, 'mlraichle.SynthStudyData');
            this.verifyClass(this.testObj.sessionData, 'mlraichle.SynthSessionData');
        end
        function test_finished(this)
            touchfile = fullfile(this.testObj.logPath, '.test_mlfourdfp.CompositeT4ResolveBuilder_isfinished.touch');
            if (lexist(touchfile)); delete(touchfile); end 
            this.verifyEqual(this.testObj.finished.finishedMarkerFilename, touchfile);           
            this.testObj.finished.touchFinishedMarker;
            this.verifyTrue(lexist(this.testObj.finished.finishedMarkerFilename, 'file'));
            this.verifyTrue(this.testObj.finished.isfinished);
            delete(touchfile); 
        end
        function test_resolve(this)
            tic
            this.testObj.NRevisions = 2;
            this.testObj = this.testObj.resolve('source', {'testNativeCropped1' 'testNativeCropped2' 'testNativeCropped3' 'testNativeCropped4' 'testNativeCropped5'});
            [e,c] = this.fvisitor.etaAndCurves('testNativeCropped3r2_Test_CompositeT4ResolveBuilder', 'testNativeCropped5r2_Test_CompositeT4ResolveBuilder');
            this.verifyEqual(e, 0.1631, 'RelTol', 0.05);
            this.verifyEqual(c, [130 133 21 89 161 66], 'RelTol', 1);
            prod = this.testObj.product;
            entr = [3.6580 3.2309];
            for p = 1:length(prod)
                this.verifyEqual(prod{p}.niftid.entropy, entr(p), 'RelTol', 0.05);
            end            
            this.verifyTrue(this.testObj.isfinished(this.testObj.sessionData));            
            fprintf('test_resolve:  ');
            toc
            if (this.view)
                prod{p}.niftid.fslview;
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
            
            in3 = 'testNativeCropped3';
            in5 = 'testNativeCropped5';
            out = sprintf('%s_on_%s', in5, in3);
            this.fvisitor.t4_mul([in3 '_to_' in5 '_t4'], [ in5 '_to_' in3 '_t4']);
            movefile([in3 '_to_' in3 '_t4'], [in3 '_to_' in3 '_t4.txt']);
            t4 = readtable([in3 '_to_' in3 '_t4'], ...
                'ReadVariableNames', 0, 'Delimiter', ' ', 'MultipleDelimsAsOne', 1, 'HeaderLines', 4, 'Format', '%f %f %f %f');
            this.verifyEqual(t4{1,1}, 1,    'AbsTol', 0.01);
            this.verifyEqual(t4{1,2}, 0.1,  'AbsTol', 0.1);
            this.verifyEqual(t4{1,3}, 0.02, 'AbsTol', 0.02);
            this.verifyEqual(t4{1,4}, 0.2,  'AbsTol', 0.2);
            
            this.fvisitor.t4img_4dfp(this.fvisitor.filenameT4(in5, in3), in5, 'out', out, 'options', ['-O' in3]);
            if (this.view)
                mlbash(sprintf('fslview %s.4dfp.img %s.4dfp.img', in3, out));
            end
        end
        function test_t4ResolveAndPaste(this)
            if (this.quick)
                return
            end
            
            this.testObj.resolveLog = loggerFilename('', 'func', 'Test_CompositeT4ResolveBuilder_test_t4ResolveAndPaste');
            [ipr_,imgFns_] = this.testObj.resolveAndPaste(this.ipr);
            this.verifyEqual(ipr_.resolved, ...
                {'testNativeCropped1_Test_CompositeT4ResolveBuilder' ...
                 'testNativeCropped2_Test_CompositeT4ResolveBuilder' ...
                 'testNativeCropped3_Test_CompositeT4ResolveBuilder' ...
                 'testNativeCropped4_Test_CompositeT4ResolveBuilder' ...
                 'testNativeCropped5_Test_CompositeT4ResolveBuilder'});
            this.verifyEqual(imgFns_, 'testNativeCropped3 testNativeCropped5');
            if (this.view)
                mlbash('fslview testNativeCropped3_Test_CompositeT4ResolveBuilder.4dfp.img testNativeCropped5_Test_CompositeT4ResolveBuilder.4dfp.img');
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
        function test_product(this)
        end
        
        %% UTILITY
        
        function test_fileprefixBlurred(this)
        end
        function test_fileprefixGaussed(this)
        end
        function test_fileprefixSumt(this)
        end
        function test_lazyStageImages(this)
            this.ipr.currentIndex = nan;
            
            fqfps = this.testObj.lazyStageImages(this.ipr);
            this.verifyEqual(length(fqfps), 2);
            this.verifyEqual(mybasename(fqfps{1}), 'testNativeCropped3');
            this.verifyEqual(mybasename(fqfps{2}), 'testNativeCropped5');
            if (this.view)
                mlbash('fslview testNativeCropped.4dfp.img -b 0,2000 testNativeCropped3.4dfp.img -b 0,2000 testNativeCropped5.4dfp.img -b 0,2000');
            end
        end
        function test_lazyBlurImage(this)
            if (this.quick)
                return
            end
            
            this.ipr.currentIndex = 1;
            fqfp = this.testObj.lazyBlurImage(this.ipr, 5.5);
            this.verifyEqual(mybasename(fqfp), 'testNativeCropped1_b55');
            if (this.view)
                mlbash('fslview testNativeCropped.4dfp.img -b 0,100 testNativeCropped1_b55.4dfp.img -b 0,100');
            end
        end
        function test_lazyBlurImages(this)
            if (this.quick)
                return
            end
            
            this.ipr.currentIndex = nan;            
            fqfps = this.testObj.lazyBlurImages(this.ipr);
            this.verifyEqual(length(fqfps), 2);
            this.verifyEqual(mybasename(fqfps{1}), 'testNativeCropped3_b55');
            this.verifyEqual(mybasename(fqfps{2}), 'testNativeCropped5_b55');
            if (this.view)
                mlbash('fslview testNativeCropped.4dfp.img -b 0,2000 testNativeCropped3_b55.4dfp.img -b 0,2000 testNativeCropped5_b55.4dfp.img -b 0,2000');
            end
        end
        function test_lazyMaskForImages(this)
            if (this.quick)
                return
            end
            
            tic
            fqfp = this.testObj.lazyMaskForImages('maskForImages', 'testNativeCropped3', 'testNativeCropped5');
            this.verifyEqual(mybasename(fqfp), 'maskForImages_testNativeCropped3_testNativeCropped5');
            if (this.view)
                mlbash('fslview testNativeCropped3.4dfp.img -b 0,2000 testNativeCropped5.4dfp.img -b 0,2000 maskForImages_3_5.4dfp.img');
            end
            fprintf('test_lazyMaskForImages:  ');
            toc
        end
        function test_frames(this)
            this.verifyEqual(sum(this.testObj.indicesLogical), 2);
        end
	end

 	methods (TestClassSetup)
		function setupCompositeT4ResolveBuilder(this)
            import mlraichle.*;
            this.studyd = SynthStudyData;
            this.sessd = SynthSessionData( ...
                'studyData', this.studyd, ...
                'sessionPath', fullfile(this.studyd.subjectsDir, this.hyglyNN, ''), ...
                'tracer', 'TEST');    
            this.pwd0 = pushd(this.sessd.tracerLocation);        
            this.testObj_ = mlfourdfp.CompositeT4ResolveBuilder( ... 
                'sessionData', this.sessd, ...  
                'NRevisions', 1, ...
                'resolveTag', this.resolveTag, ...   
                'theImages', this.theImages, ...
                'indicesLogical', this.indicesLogical, ... 
                'indexOfReference', 3);
            this.ipr = struct( ...
                'dest', [], ...
                'source', [], ...
                'maskForImages', 'maskForImages', ...
                'sourceBlur', [], ...
                'destBlur', [], ...
                'keepForensics', true, ...
                'currentIndex', nan);
            this.ipr.dest = this.theImages;
            this.ipr.source = this.theImages;
            this.ipr.sourceBlur = this.blurArray;
            this.ipr.destBlur = this.blurArray;
            this.fvisitor = mlfourdfp.FourdfpVisitor;
            this.ensureImages;
            setenv('DEBUG', ''); % cf. dbbash
 			this.addTeardown(@this.teardownClass);
 		end
	end

 	methods (TestMethodSetup)
		function setupCompositeT4ResolveBuilderTest(this)
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
        function f = indicesLogical(~)
            f = false(1, 5);
            f(3) = true;           
            f(5) = true;
        end
        function ensureImages(this)
            if (~this.fvisitor.lexist_4dfp('testNativeCropped'))
                this.fvisitor.copy_4dfp('TEST_V1-LM-00-OP', 'testNativeCropped');
            end
            for idx = 1:5
                if (~this.fvisitor.lexist_4dfp(sprintf('testNativeCropped%i', idx)))
                    this.fvisitor.extract_frame_4dfp('testNativeCropped', idx);
                    this.fvisitor.move_4dfp(sprintf('testNativeCropped_frame%i', idx), ...
                                            sprintf('testNativeCropped%i', idx));
                end
            end
        end
		function teardownMethod(~)
 		end
		function teardownClass(this)
            if (~this.quick)
                this.testObj_.keepForensics = false;
                this.testObj_.teardownResolve(this.ipr);
                %delete('*Test_CompositeT4ResolveBuilder.4dfp.*');
                delete('*.nii');
                delete('*.nii.gz');
                if (isdir('Log')); rmdir('Log', 's'); end
                if (isdir('T4'));  rmdir('T4', 's');  end
            end
            delete('*.mat0');
            delete('*.sub');
            cd(this.pwd0);
 		end
	end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

