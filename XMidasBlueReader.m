classdef XMidasBlueReader < handle
    %XMIDASBLUEREADER Read X-Midas BLUE file
    %   Reads X-Midas blue files
    
    properties
        hcb
        ext_header
    end
    
    properties(SetAccess = private)
        bluefile
        
        dataOffset
    end
    
    methods
        function obj = XMidasBlueReader(bluefile)
            %XMIDASBLUEREADER Construct a reader for the bluefile
            %   Creates a class that progressively reads samples from
            %   an X-Midas BLUE file
            obj.bluefile = bluefile;
            
            % Open the file, read headers
            [fid, msg] = fopen(obj.bluefile, 'r');
            if fid < 0
                disp(msg)
                return
            end
            obj.hcb = XMidasBlueReader.HCB(fid);
            obj.ext_header = XMidasBlueReader.ExtendedHeader(fid, obj.hcb);
            obj.resetRead();
            
            % Close file.
            fclose(fid);
        end
        
        function [samples] = read(obj, numSamples)
            %READ Read numSamples from the current position.
            %   Returns a vector or matrix of the data according to the
            %   hcb information.
            if nargin < 2
                numSmaples = 1;
            end
            
            % Open the file, skip to read pointer
            [fid, msg] = fopen(obj.bluefile, 'r');
            if fid < 0
                throw(MException(...
                    'XMidasBlueReader:Read', ...
                    sprintf('Unable to read file: %s', msg)));
            end
            
            
            % Convert format size to number of elements per sample
            elementsPerSample = XMidasBlueReader.FormatSize(obj.hcb.format(1));
            % Convert format type to the type string and size (bytes)
            [elementPrecisionType, elementPrecisionBytes] = ...
                XMidasBlueReader.FormatType(obj.hcb.format(2));
            
            % Update effectivePrecisionType for array inputs
            elementPrecisionEffectiveType = elementPrecisionType;
            if prod(elementsPerSample) > 1
                % Reformats 'char' to '8*char' if there are 8 elements per
                % sample, for example.
                elementPrecisionEffectiveType = sprintf(...
                    '%d*%s', prod(elementsPerSample), precision_type);
            end
            
            % Move read position.
            fseek(fid, obj.dataOffset, 'bof');
            
            % From header and sample calculations, calculate number of
            % bytes per sample, check remaining data length, then read if
            % possible.
            bytesPerSample = elementPrecisionBytes * prod(elementsPerSample);
            bytesRead = obj.dataOffset - obj.hcb.data_start;
            bytesRemaining = obj.hcb.data_size - bytesRead;
            if bytesPerSample * numSamples > bytesRemaining
                throw(MException(...
                    'XMidasBlueReader:Read', ...
                    sprintf(...
                        'Requested read is longer than remaining: (%d vs. %d)', ...
                        bytesPerSample * numSamples, bytesRemaining) ...
                        ));
            end

            % Read and move read pointer
            if prod(elementsPerSample) == 1
                samples = fread(fid, numSamples, elementPrecisionEffectiveType);
            else
                samples = fread(fid, numSamples, elementPrecisionEffectiveType, elementPrecisionBytes);
            end
            obj.dataOffset = obj.dataOffset + (bytesPerSample * numSamples);
            
            % Close file
            fclose(fid);
        end
        
        function resetRead(obj)
            %RESETREAD Resets the read pointer back to the start
            %    Moves data read pointer to the start of data.
            obj.dataOffset = obj.hcb.data_start;
        end
    end
    
    methods(Static)
        function [hcb] = HCB(fid)
            %HCB Reads a hcb control block structure
            %    Uses the file descriptor to read the hcb Control Block of an
            %    X-Midas BLUE file.  If no file handle is provided, it
            %    returns a generic HCB structure.
            hcb = struct( ...
                'version', '', ...
                'head_rep', '', ...
                'data_rep', '', ...
                'detached', 0, ...
                'protected', 0, ...
                'pipe', 0, ...
                'ext_start', 0, ...
                'ext_size', 0, ...
                'data_start', 0.0, ...
                'data_size', 0.0, ...
                'type', 0, ...
                'format', '', ...
                'flagmask', 0, ...
                'timecode', 0.0, ...
                'inlet', 0, ...
                'outlets', 0, ...
                'outmask', 0, ...
                'pipeloc', 0, ...
                'pipesize', 0, ...
                'in_byte', 0.0, ...
                'out_byte', 0.0, ...
                'outbytes', [], ...
                'keylength', 0, ...
                'keywords', int8([]), ...
                'adjunct', int8([]) ...
                );
            if fid < 0
                return
            end
            
            % Read the descriptor
            fseek(fid, 0, 'bof');
            hcb.version = fread(fid, 4, '*char')';
            hcb.head_rep = fread(fid, 4, '*char')';
            hcb.data_rep = fread(fid, 4, '*char')';
            hcb.detached = fread(fid, 1, 'int32');
            hcb.protected = fread(fid, 1, 'int32');
            hcb.pipe = fread(fid, 1, 'int32');
            hcb.ext_start = fread(fid, 1, 'int32');
            hcb.ext_size = fread(fid, 1, 'int32');
            hcb.data_start = fread(fid, 1, 'double');
            hcb.data_size = fread(fid, 1, 'double');
            hcb.type = fread(fid, 1, 'int32');

            % Check the type for ones we support.
            if ~ismember(hcb.type, [1000, 2000])
                e = MException(...
                    'BLUE:HCBType', ...
                    sprintf('Unsupported HCB.type: %d', hcb.type));
                throw(e);
            end

            hcb.format = fread(fid, 2, '*char')';
            hcb.flagmask = fread(fid, 1, 'int16');
            hcb.timecode = fread(fid, 1, 'double');
            hcb.inlet = fread(fid, 1, 'int16');
            hcb.outlets = fread(fid, 1, 'int16');
            hcb.outmask = fread(fid, 1, 'int32');
            hcb.pipeloc = fread(fid, 1, 'int32');
            hcb.pipesize = fread(fid, 1, 'int32');
            hcb.in_byte = fread(fid, 1, 'double');
            hcb.out_byte = fread(fid, 1, 'double');
            hcb.outbytes = fread(fid, 8, '8*double', 8); % 8x8 byte values
            hcb.keylength = fread(fid, 1, 'int32');
            hcb.keywords = fread(fid, 92, '*char')'; % 92 1-byte values

            % Interpret the adjunct
            switch hcb.type
                case 1000
                    hcb.adjunct = struct(...
                        'xstart', fread(fid, 1, 'double'), ...
                        'xdelta', fread(fid, 1, 'double'), ...
                        'xunits', fread(fid, 1, 'int32') ...
                        );
                case 2000
                    hcb.adjunct = struct(...
                        'xstart',  fread(fid, 1, 'double'), ...
                        'xdelta',  fread(fid, 1, 'double'), ...
                        'xunits',  fread(fid, 1, 'int32'), ...
                        'subsize', fread(fid, 1, 'int32'), ...
                        'ystart',  fread(fid, 1, 'double'), ...
                        'ydelta',  fread(fid, 1, 'double'), ...
                        'yunits',  fread(fid, 1, 'int32') ...
                        );
            end
        end % end HCB
    
        function [ext_header] = ExtendedHeader(fid, hcb)
            if fid < 0
                return;
            end
            
            
            % Move read pointer to ext_start*512
            fseek(fid, hcb.ext_start * 512, 'bof');
            
            ext_header = [];
            bytesRemaining = hcb.ext_size;            
            while bytesRemaining > 0
                key = struct(...
                    'lkey', fread(fid, 1, 'int32'),...  % 4 bytes
                    'lext', fread(fid, 1, 'int16'),...  % 2 bytes
                    'ltag', fread(fid, 1, 'int8'), ...  % 1 byte
                    'type', fread(fid, 1, '*char'), ... % 1 Byte
                    'value', 0, ... % Variable, lkey-lext * sizeof(type)
                    'tag', '');     % Variable, ltag bytes
                
                [valType, valTypeBytes] = XMidasBlueReader.FormatType(key.type);
                valCount = (key.lkey - key.lext) / valTypeBytes;
                key.value = fread(fid, valCount, valType)';
                key.tag = fread(fid, key.ltag, '*char')';
                
                % Calculate remainder for 8-byte boundary
                total = 4 + 2 + 1 + 1 + (key.lkey - key.lext) + key.ltag;
                remainder = 8 - rem(total, 8);
                if remainder > 0
                    fseek(fid, remainder, 'cof');
                end
                
                % Subtract from bytes remaining and add to keywords
                bytesRemaining = bytesRemaining - key.lkey;
                ext_header = [ext_header, key];
            end
        end
        
        function [typeStr, bytes] = FormatType(formatType)
            %FORMATTYPE Converts the format type character
            %     Convert format type to the fread string name and size in
            %     bytes.
            switch formatType
                case 'B'
                    % 8-bit int
                    typeStr = 'int8';
                    bytes = 1;
                case 'I'
                    % 16-bit int
                    typeStr = 'int16';
                    bytes = 2;
                case 'L'
                    % 32-bit int
                    typeStr = 'int32';
                    bytes = 4;
                case 'X'
                    % 64-bit int
                    typeStr = 'int64';
                    bytes = 8;
                case 'F'
                    % 32-bit float
                    typeStr = 'single';
                    bytes = 4;
                case 'D'
                    % 64-bit float
                    typeStr = 'double';
                    bytes = 8;
                case 'A'
                    % ASCII 8 characters or variable length
                    typeStr = '*char';
                    bytes = 1;
                otherwise
                    e = MException(...
                        'XMidasBlueReader:HCB:Format', ...
                        sprintf('Unsupported HCB format type: %s', formatType));
                    throw(e);
            end
        end % FormatType
        
        function [eps] = FormatSize(formatSize)
            %FORMATSIZE Converts HCB format size
            %    Convert format size to number of elements per sample
            switch formatSize
                case 'S'
                    % Scalar
                    eps = 1;
                case 'C'
                    % Complex (real, imaginary pairs)
                    eps = 2;
                case 'V'
                    % Vector (x,y,z)
                    eps = 3;
                case 'Q'
                    % Quad (x,y,z, time)
                    eps = 4;
                case 'M'
                    % Matrix 3x3 (9 elements)
                    eps = [3,3];
                case 'T'
                    % Transform Matrix 4x4 (16 elements)
                    eps = [4,4];
                case '1'
                    % 1 element per point
                    eps = 1;
                case '2'
                    % 2 elements per point
                    eps = 2;
                case '3'
                    % 3 elements per point
                    eps = 3;
                case '4'
                    % 4 elements per point
                    eps = 4;
                case '5'
                    % 5 elements per point
                    eps = 5;
                case '6'
                    % 6 elements per point
                    eps = 6;
                case '7'
                    % 7 elements per point
                    eps = 7;
                case '8'
                    % 8 elements per point
                    eps = 8;
                case '9'
                    % 9 elements per point
                    eps = 9;
                case 'X'
                    % 10 elements per point
                    eps = 10;
                case 'A'
                    % 32 elements per point
                    eps = 32;
                otherwise
                    e = MException(...
                        'XMidasBlueReader:HCB:Format', ...
                        sprintf('Unsupported HCB format size: %s', formatSize));
                    throw(e);
            end
        end
    end
end

