classdef ImageFrames < mlfourdfp.AbstractImageComposite
	%% IMAGEFRAMES  

	%  $Revision$
 	%  was created 10-Jan-2017 22:15:44
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.  Copyright 2017 John Joowon Lee.
 	

    properties
        fractionalImageFrameThresh = 0 % eps is most permissive; 1 is most restrictive
    end
    
    properties (Dependent)
        referenceImage
        sessionData   
        sourceImage
        sourceImageTable
        theImages
    end
    
    methods %% GET, SET
        function g    = get.referenceImage(this)
            g = this.it4ResolveBuilder_.fileprefixIndexed( ...
                this.sessionData.tracerRevision('typ', 'fqfp'), this.indexOfReference);
        end
        function g    = get.sessionData(this)
            g = this.it4ResolveBuilder_.sessionData;
        end
        function g    = get.sourceImage(this)
            if (~isempty(this.theImages))
                g = this.theImages;
                return
            end
            g = this.sessionData.tracerSif('typ', 'fqfp');
        end
        function g    = get.sourceImageTable(this)
            assert(isa(this.sourceImageTable_, 'table'));
            g = this.sourceImageTable_;
        end        
        function g    = get.theImages(this)
            if (isempty(this.theImages_))
                g = this.lazyExtractImages(struct('dest', this.it4ResolveBuilder_.sessionData.tracerRevision('typ', 'fqfp')));
                return
            end
            g = this.theImages_;
        end
    end
    
	methods
 		function this  = ImageFrames(varargin)
 			%% IMAGEFRAMES
 			%  Usage:  this = ImageFrames()

            this = this@mlfourdfp.AbstractImageComposite(varargin{:});            
 			ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired( ip, 'it4rb', @(x) isa(x, 'mlfourdfp.IT4ResolveBuilder'));
            addParameter(ip, 'theImages', {}, @(x) iscell(x) || ischar(x));
            addParameter(ip, 'indicesLogical', true, @islogical);
            parse(ip, varargin{:}); 
            
            import mlfourdfp.*;
            this.theImages_ = FourdfpVisitor.ensureSafeFileprefix(ip.Results.theImages);
            if (isempty(this.theImages))
                this.length_ = this.readLength;
            else
                this.length_ = this.readLength(this.theImages);
            end
            if (this.fractionalImageFrameThresh < eps)
                this.nonEmptyImageIndices_ = true(1, this.length_);
            else
                this.nonEmptyImageIndices_ = this.nonEmptyImageIndices;
            end
            this.indicesLogical = ip.Results.indicesLogical;
            [this.indexMin_,this.indexMax_] = this.findIndexBounds;
            %this.sourceImageTable_ = this.readSourceImageTable__;
        end        
        function fqfp  = lazyExtractImage(this, ipr)
            %% LAZYEXTRACTIMAGE uses specifiers in ipr; will not replace any existing frame
            
            fqfp = this.it4ResolveBuilder_.fileprefixIndexed(ipr.dest, ipr.currentIndex);
            bv = mlfourdfp.FourdfpVisitor;
            if (~bv.lexist_4dfp(fqfp))
                bv.extract_frame_4dfp(ipr.dest, ipr.currentIndex, ['-o' fqfp]);
            end
        end
        function fqfps = lazyExtractImages(this, ipr)
            %% LAZYEXTRACTIMAGES uses specifiers in ipr; will not replace any existing indicesLogical
            
            fqfps = {};
            for f = 1:length(this.indicesLogical)
                if (this.indicesLogical(f))
                    ipr.currentIndex = f;
                    fqfps = [fqfps this.lazyExtractImage(ipr)]; %#ok<AGROW>
                end
            end
            assert(length(fqfps) == sum(this.indicesLogical));
        end
        function len   = length(this)
            len = this.length_;
        end
        function fr    = nonEmptyImageIndices(this, varargin)
            %% NONEMPTYIMAGEINDICES
            %  @param fqfn is the filename for the dynamic tracer image;
            %  default is this.it4ResolveBuilder.sessionData.tracerSif('typ', '4dfp.ifh').
            %  @param named fracThresh is < 1; 
            %  default is this.fractionalImageFrameThresh.
            %  @returns fr, a binary vector indicating nonempty indicesLogical
            %  @returns this
            
            ip = inputParser;
            addOptional(ip, 'fqfn', [this.sourceImage '.4dfp.ifh'], @(x) lexist(x, 'file'));
            addParameter(ip, 'fracThresh', this.fractionalImageFrameThresh, @isnumeric);
            parse(ip, varargin{:});
            
            if (~isempty(this.nonEmptyImageIndices_))
                fr = this.nonEmptyImageIndices_;
                return
            end
            
            cache = [myfileprefix(ip.Results.fqfn) '_nonEmptyImageIndices.nii.gz'];
            if (lexist(cache))
                tr = mlfourd.NIfTId.load(cache);
                fr = tr.img';
                return
            end
            
            tr = mlfourd.NumericalNIfTId.load(ip.Results.fqfn);
            tr = tr.volumeSummed;
            tr = tr > ip.Results.fracThresh*max(tr.img);
            tr.img = double(tr.img);
            tr.saveas(cache);
            fr = tr.img';
        end
        function len   = readLength(this, varargin)
            ip = inputParser;
            addOptional(ip, 'tracerSif', this.sourceImage, @lexist_4dfp);
            parse(ip, varargin{:});
            
            [~,len] = mlbash(sprintf('awk ''/matrix size \\[4\\]/{print $NF}'' %s.4dfp.ifh', ip.Results.tracerSif));
            len = str2double(len);
        end
        function mmpp  = readScalingFactors(this, varargin)
            ip = inputParser;
            addRequired(ip, 'fqfp', this.sourceImage, @lexist_4dfp);
            parse(ip, varargin{:});
            
            mmpp = zeros(1,3);
            for mu = 1:3
                [~,f] = mlbash(sprintf('awk ''/scaling factor \\(mm\\/pixel\\) \\[%i\\]/{print $NF}'' %s.4dfp.ifh', mu, ip.Results.fqfp));
                mmpp(mu) = str2double(f);
            end
        end
        function tab   = readSourceImageTable__(this, varargin)
            ip = inputParser;
            addOptional(ip, 'tracerSif', [this.sourceImage '.4dfp.img.rec'], @(x) lexist(x, 'file'));
            parse(ip, varargin{:});
            
            tab = readtable(...
                ip.Results.tracerSif, ...
                    detectImportOptions(ip.Results.tracerSif, 'FileType', 'text'), ...
                        'ReadVariableNames', true, 'ReadRowNames', true);
            tab = tab(1:end-1,:);
        end
    end 
    
    %% PROTECTED
    
    properties (Access = protected)
        length_
        nonEmptyImageIndices_
        sourceImageTable_
        theImages_
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

