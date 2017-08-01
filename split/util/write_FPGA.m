function suc = write_FPGA(sd,fid)

% scriptname = sprintf('user_probdata_FPGA.h');
% fid = fopen(scriptname,'wt');
% 
% fprintf(fid,'#ifndef __user_probData_FPGA_h__ \n');
% fprintf(fid,'#define __user_probData_FPGA_h__ \n');
% fprintf(fid, '#include <stdio.h>\n');
% fprintf(fid, '#include <stdlib.h>\n');
% fprintf(fid, '#include <string.h>\n');
for i=1:length(sd.data)
    
    if length(sd.data(i).x) > 1
        
        if (issparse(sd.data(i).x) == 1)
            sd.data(i).x = full(sd.data(i).x);
        end
        fprintf(fid,write_data_h('type_arg',sd.data(i).type,'1d', sd.data(i).x, sd.data(i).name ) );
        fprintf(fid,'\n');
    elseif length(sd.data(i).x) == 1
        fprintf(fid,write_data_h('type_arg',sd.data(i).type,'1d', sd.data(i).x, sd.data(i).name ) );
        fprintf(fid,'\n');
    end
end

% fprintf(fid,'#endif');
% fclose(fid);
suc = 1;
end