function fix_binary_3VCPI(infilename, outfilename)
%% This file fixes an issue with the new 5441X netburner code which writes a single extra word at the start of data communications.
% This file removes the first word after the 8 word MS timestamp, which are at the start of every 2048 word block(8 word timestamp added by Chris W's code).
% Then it moves the checksum that occurs for every 2048 word block from after the MS timestamp to before the MS timestamp
% THe resultsing file  is written to the same file name with "_reorder.2DSCPI_wordbegone" appended.

% Written by Chris R.  7/12/2019
% Rev 1.
% Modified by Joe Finlon 4/15/2020 for UIOOPS support

% Inputs:
%   infilename: full path to raw *.2DSCPI file
%   outfilename: full path to revized file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Open original *.2DSCPI file

fid = fopen(infilename, 'r');

%% Open a new file for data without timestamps

fid3 = fopen(outfilename, 'w'); % open the xls file


if fid > 0 
    %Read 2057 words = 4096 bytes + 2 byte checksum + 16 bytes of timestamp
    [data, cntfile] = fread(fid, 9, 'uint16');%,readloop);
    fwrite(fid3, data(1:8), 'uint16');        
    [data,cntfile] = fread(fid, 2057, 'uint16');%,readloop);
    while (cntfile == 2057)
        fwrite(fid3, data(1:2048), 'uint16'); 
        fwrite(fid3, data(2057), 'uint16'); 
        fwrite(fid3, data(2049:2056), 'uint16'); 

        [data,cntfile] = fread(fid,2057,'uint16');       
    end %end of "while (cntfile == 2057)"

    fclose(fid);
    fclose(fid3);
end % end if output to excel file

end


