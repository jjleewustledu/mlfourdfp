classdef SCANS_studies_txt
    %% 
    %
    %  Created 17-May-2022 16:21:25 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlfourdfp/src/+mlfourdfp.
    %  Developed on Matlab 9.10.0.1851785 (R2021a) Update 6 for MACI64.  Copyright 2022 John J. Lee.
    
    methods
        function this = SCANS_studies_txt(varargin)
            %% SCANS_STUDIES_TXT should be constructed with cwd containing SCANS.studies.txt
            %  Args:
            %      filename (file): of SCANS.studies.txt.
            
            ip = inputParser;
            addOptional(ip, "filename", 'SCANS.studies.txt', @istext)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            try
                this.the_table_ = readtable(ipr.filename);
            catch ME
                handwarning(ME);
                try
                    ipr.filename = 'scans.studies.txt';
                    this.the_table_ = readtable(ipr.filename);
                catch ME1
                    handwarning(ME1);
                    ipr.filename = fullfile(pwd, '.studies.txt');
                    this.the_table_ = readtable(ipr.filename);
                end
            end
            this.the_table_.Properties.VariableNames = {'scan_index', 'pulse_seq', 'scan_name', 'num_slices'};
        end
        
        function fqfn = dcm2niix(this, varargin)
            ip = inputParser;
            addParameter(ip, 'seq_pattern', '', @istext);
            addParameter(ip, 'name_pattern', '', @istext);
            addParameter(ip, 'match_exactly', false, @islogical);
            addParameter(ip, 'case_sensitive', true, @islogical);
            addParameter(ip, 'first', false, @islogical);
            addParameter(ip, 'last', false, @islogical);
            addParameter(ip, 'o', pwd, @isfolder);
            parse(ip, varargin{:});
            ipr = ip.Results;
            
            % find study* folder
            if ipr.case_sensitive
                t = this.table();
            else
                ipr.seq_pattern = lower(ipr.seq_pattern);
                ipr.name_pattern = lower(ipr.name_pattern);
                t = table(this.table.scan_index, ...
                          lower(this.table.pulse_seq), ...
                          lower(this.table.scan_name), ...
                          this.table.num_slices, ...
                          'VariableNames', {'scan_index', 'pulse_seq', 'scan_name', 'num_slices'});
            end

            if ipr.match_exactly
                found = strcmp(t.scan_name, ipr.name_pattern) & strcmp(t.pulse_seq, ipr.seq_pattern);                
            else
                found = contains(t.scan_name, ipr.name_pattern) & contains(t.pulse_seq, ipr.seq_pattern);
            end  
            si =  t.scan_index(found);
            if ipr.first
                si = si(1);
            end
            if ipr.last
                si = si(end);
            end
            if ~ipr.first && ~ipr.last
                [~,idx] = max(t.num_slices(found)); % select best z-resolution
                si = si(idx);
            end
            study = sprintf('study%i', si);
            assert(isfolder(study), 'expected study folder %s not found', study)
            
            % dcm2niix
            mlpipeline.Bids.dcm2niix(study, 'o', ipr.o, 'terse', true);
            
            % find fqfn
            g = glob(sprintf('%s%s*-%i.nii.gz', ipr.o, filesep, si));
            fqfn = fullfile(ipr.o, g{end});
            fprintf('%s wrote %s\n', clientname(false, 2), fqfn)
        end
        function [patt,fqfn] = dcm2niix_flair(this, varargin)
            try
                patt = 'FLAIR';
                fqfn = this.dcm2niix(varargin{:}, 'name_pattern', patt, 'case_sensitive', false);
            catch
                patt = 'T2';
                fqfn = this.dcm2niix(varargin{:}, 'seq_pattern', 'tse2d', 'name_pattern', patt, 'case_sensitive', false);
            end
        end
        function [patt,fqfn] = dcm2niix_mpr(this, varargin)
            patt = 'MPR';
            try
                fqfn = this.dcm2niix(varargin{:}, 'name_pattern', 'SAGT1MPR', 'case_sensitive', false);
            catch
                try
                    fqfn = this.dcm2niix(varargin{:}, 'name_pattern', 'T1MPR', 'case_sensitive', false);
                catch
                    fqfn = this.dcm2niix(varargin{:}, 'name_pattern', 'MPR', 'case_sensitive', false);
                end
            end
        end
        function [patt,fqfn] = dcm2niix_postGd(this, varargin)
            try
                patt = 'STEALTH';
                fqfn = this.dcm2niix(varargin{:}, 'name_pattern', 'T1GRESTEALTH');
            catch
                try
                    patt = 'STEALTH';
                    fqfn = this.dcm2niix(varargin{:}, 'name_pattern', 'TRA3DSTEALTH');
                catch
                    try
                        patt = 'POST';
                        fqfn = this.dcm2niix(varargin{:}, 'name_pattern', 'T1POST');
                    catch
                        try                        
                            patt = 'POST';
                            fqfn = this.dcm2niix(varargin{:}, 'name_pattern', 'T13DPOST');
                        catch 
                            patt = 'TRA_3D_GRE';
                            fqfn = this.dcm2niix(varargin{:}, 'name_pattern', 'TRA3DGRE', 'last', true);
                        end
                    end
                end
            end
        end
        function t = table(this)
            t = this.the_table_;
        end
    end
    
    properties (Access = private)        
        the_table_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
