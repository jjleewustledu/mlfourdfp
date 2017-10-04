classdef Test_ImgRecLogger < matlab.unittest.TestCase
	%% TEST_IMGRECLOGGER 

	%  Usage:  >> results = run(mlfourdfp_unittest.Test_ImgRecLogger)
 	%          >> result  = run(mlfourdfp_unittest.Test_ImgRecLogger, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 19-Sep-2017 18:04:33 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/test/+mlfourdfp_unittest.
 	%% It was developed on Matlab 9.2.0.538062 (R2017a) for MACI64.  Copyright 2017 John Joowon Lee.
 	
	properties
        pwd0 = '~/Local/src/mlcvl/mlfourdfp/test/+mlfourdfp_unittest'
 		registry
 		testObj
 	end

	methods (Test)
        function test_cons(this)
 			import mlfourdfp.*;
            this.testObj.cons('preface 1');
            this.testObj.cons({'preface 2 line 1' 'preface 2 line 2'});
            this.testObj.cons('preface 3');
            this.testObj.save;
        end
        function test_consForTracers(this)
 			import mlfourdfp.*;
            hdrLine = 'Frame     	Length(msec)	Midpoint(sec)	Start(msec)	 Frame_Min	 Frame_Max	 Decay_Fac	Rescale';
            irp0 = ImgRecParser.loadx('test_ImgRecLogger0.4dfp.img.rec', '.4dfp.img.rec');
            irp  = ImgRecParser.loadx('test_ImgRecLogger.4dfp.img.rec', '.4dfp.img.rec');
            [~,hdr] = irp.findNextCell(hdrLine, 2);
            [~,ftr] = irp.findNextCell('endrec', hdr);
            irp2 = ImgRecParser.loadx('test_ImgRecLogger2.4dfp.img.rec', '.4dfp.img.rec');
            [~,hdr2] = irp.findNextCell(hdrLine, 2);
            [~,ftr2] = irp.findNextCell('endrec', hdr2);
            irp3 = ImgRecParser.loadx('test_ImgRecLogger3.4dfp.img.rec', '.4dfp.img.rec');
            [~,hdr3] = irp.findNextCell(hdrLine, 2);
            [~,ftr3] = irp.findNextCell('endrec', hdr3);
            
            this.testObj.cons(cell2str(irp0.cellContents));
            this.testObj.consNoHeadFoot(cell2str(irp3.cellContents(hdr3+1:ftr3-1)));
            this.testObj.consNoHeadFoot(cell2str(irp2.cellContents(hdr2+1:ftr2-1)));
            this.testObj.consNoHeadFoot(cell2str(irp.cellContents(hdr+1:ftr-1)));
            this.testObj.cons(hdrLine);
            this.testObj.save;
            
 		end
	end

 	methods (TestClassSetup)
		function setupImgRecLogger(this)
 			import mlfourdfp.*;
            cd(this.pwd0);
 			this.testObj_ = ImgRecLogger('test_ImgRecLogger4');
 			%this.addTeardown(@this.cleanFiles);
 		end
	end

 	methods (TestMethodSetup)
		function setupImgRecLoggerTest(this)
 			this.testObj = this.testObj_;
 		end
	end

	properties (Access = private)
 		testObj_
 	end

	methods (Access = private)
		function cleanFiles(this)
            deleteExisting(fullfile(this.pwd0, 'test_ImgRecLogger4.4dfp.img.rec'));
 		end
	end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

