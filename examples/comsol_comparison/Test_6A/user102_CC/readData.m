function [ dat ] = readData(fname,nskip)
    fid = fopen(fname);
    
    dat = cell2mat(textscan(fid, '%f%f%f%f%f%f%f', 'headerlines', nskip));   %%assumes 7 columns
%     
%     % read the first 36 lines and do nothing to them
%     for k=1:nskip
%         tline = fgets(fid);
%     end
% 
%     dat = [];
%     cnt = 1;
%     while ischar(tline)
%         dat(cnt,:) = fgets(fid);
% %         tline = fgets(fid);
%         cnt = cnt + 1;
%     end
     fclose(fid);
end

