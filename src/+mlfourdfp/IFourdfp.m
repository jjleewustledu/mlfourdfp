classdef IFourdfp 
	%% IFOURDFP provide a minimal set of imaging properties, methods 

	%  $Revision$
 	%  was created 25-Apr-2016 23:01:59
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.0.0.341360 (R2016a) for MACI64.
 	
    
    properties (Constant) 
        FILETYPE     = '4DFP'
        FILETYPE_EXT = '.4dfp.img'
        SUPPORTED_EXT = {'.4dfp.ifh' '.4dfp.hdr' '.4dfp.img' '.4dfp'}
    end
    
	properties (Abstract)
%         img
%         
%         bitpix 
%         creationDate
%         datatype
%         descrip
%         entropy
%         hdxml
%         label
%         machine
%         mmppix
%         negentropy
%         orient
%         pixdim
%         seriesNumber        
    end 
    
	methods (Abstract) 
        
        %% for NIfTId and other concrete imaging classes
        
        %this = clone(this)
        %[tf,msg] = isequal(this, n)
        %[tf,msg] = isequaln(this, n)
        %this = makeSimilar(this)
        
        %% for AbstractNIfTId and other abstract imaging classes
        
%         char(this)
%         append_descrip(this, s)
%         prepend_descrip(this, s)
%         double(this)
%         duration(this)
%         append_fileprefix(this, s)
%         prepend_fileprefix(this, s)
%         fov(this)
%         matrixsize(this)
%         %ones(this)
%         rank(this)
%         scrubNanInf(this)
%         single(this)
%         size(this)
%         %zeros(this)
    end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

