classdef ImageFrames < mlfourdfp.AbstractImageComposite
	%% IMAGEFRAMES  

	%  $Revision$
 	%  was created 10-Jan-2017 22:15:44
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.  Copyright 2017 John Joowon Lee.
    
    properties (Dependent)
        fractionalImageFrameThresh
        referenceImage
        sessionData   
        sourceImage
        sourceImageTable
        theImages
    end
    
    methods 
        
        %% GET, SET
        
        function g    = get.fractionalImageFrameThresh(this)
            g = this.sessionData.fractionalImageFrameThresh;
        end
        function g    = get.referenceImage(this)
            ipr = struct( ...
                'dest', this.sessionData.tracerRevision('typ', 'fqfp'), ...
                'currentIndex', this.indexOfReference);
            g = this.it4ResolveBuilder_.fileprefixIndexed(ipr);
        end
        function g    = get.sessionData(this)
            g = this.it4ResolveBuilder_.sessionData;
        end
        function g    = get.sourceImage(this)
            if (~isempty(this.theImages))
                g = this.theImages;
                return
            end
            g = this.sessionData.tracerRevision('typ', 'fqfp');
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
        function this = set.theImages(this, s)
            assert(iscell(s) || ischar(s));
            this.theImages_ = mlfourdfp.FourdfpVisitor.ensureSafeFileprefix(s);
            this.length_ = this.readLength(this.theImages);
            this.indicesLogical = true;
            [this.indexMin_,this.indexMax_] = this.findIndexBounds;
        end
        
        %%
        
 		function this  = ImageFrames(varargin)
 			%% IMAGEFRAMES
 			%  Usage:  this = ImageFrames()

            this = this@mlfourdfp.AbstractImageComposite(varargin{:});            
 			ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired( ip, 'it4rb', @(x) isa(x, 'mlfourdfp.IT4ResolveBuilder'));
            addParameter(ip, 'theImages', {}, @(x) iscell(x) || ischar(x));
            addParameter(ip, 'indicesLogical', this.sessionData.indicesLogical, @islogical);
            parse(ip, varargin{:}); 
            
            this.theImages_ = this.it4ResolveBuilder_.ensureSafeFileprefix(ip.Results.theImages);
            if (isempty(this.theImages))
                this.length_ = this.readLength(this.sourceImage);
            else
                this.length_ = this.readLength(this.theImages);
            end
            this.indicesLogical = ip.Results.indicesLogical;
            [this.indexMin_,this.indexMax_] = this.findIndexBounds;
            %this.sourceImageTable_ = this.readSourceImageTable__;
        end       
        
        function fqfp  = lazyExtractImage(this, ipr)
            %% LAZYEXTRACTIMAGE uses specifiers in ipr; will not replace any existing frame
            
            fqfp = this.it4ResolveBuilder_.fileprefixIndexed(ipr);
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
        function len   = readLength(~, tracerSif)
            if iscell(tracerSif)
                tracerSif = tracerSif{1};
            end
            assert(istext(tracerSif));
            tracerSif = myfileprefix(tracerSif);
            
            [~,len] = mlbash(sprintf('awk ''/matrix size \\[4\\]/{print $NF}'' %s.4dfp.ifh', tracerSif));
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

