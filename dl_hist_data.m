function [Daily_Data] = dl_hist_data(stock_symbol,start_year)
% Get current date
[this_year, this_month, this_day] = datevec(now);

% Build URL string
url_string = 'http://ichart.finance.yahoo.com/table.csv?';
url_string = strcat(url_string, '&s=', upper(stock_symbol)   );
url_string = strcat(url_string, '&d=', num2str(this_month-1) );
url_string = strcat(url_string, '&e=', num2str(this_day)     );
url_string = strcat(url_string, '&f=', num2str(this_year)    );
url_string = strcat(url_string, '&g=d&a=0&b=1&c=', start_year);
url_string = strcat(url_string, '&ignore.csv');

% Open a connection to the URL and retrieve data into a buffer
buffer      = java.io.BufferedReader(...
    java.io.InputStreamReader(...
    openStream(...
    java.net.URL(url_string))));
dummy = readLine(buffer); % Read the first line (a header) and discard

% Read all remaining lines in buffer
ptr = 1;
while 1
    % Read line
    buff_line = char(readLine(buffer));
    
    % Break if this is the end
    if length(buff_line)<3, break; end
    
    % Find comma delimiter locations
    commas    = find(buff_line== ',');
    
    % Extract high, low, open, close, etc. from string
    DATEvar   = buff_line(1:commas(1)-1);
    OPENvar   = str2double( buff_line(commas(1)+1:commas(2)-1) );
    HIGHvar   = str2double( buff_line(commas(2)+1:commas(3)-1) );
    LOWvar    = str2double( buff_line(commas(3)+1:commas(4)-1) );
    CLOSEvar  = str2double( buff_line(commas(4)+1:commas(5)-1) );
    VOLvar    = str2double( buff_line(commas(5)+1:commas(6)-1) );
    adj_close = str2double( buff_line(commas(6)+1:end) );
    
    %Adjust for dividends, splits, etc.
    DATEtemp{ptr,1} = DATEvar;
    OPENtemp(ptr,1) = OPENvar  * adj_close / CLOSEvar;
    HIGHtemp(ptr,1) = HIGHvar  * adj_close / CLOSEvar;
    LOWtemp (ptr,1) = LOWvar   * adj_close / CLOSEvar;
    CLOSEtemp(ptr,1)= CLOSEvar * adj_close / CLOSEvar;
    VOLtemp(ptr,1)  = VOLvar;
    
    ptr = ptr + 1;
end
Daily_Data = [datenum(DATEtemp(:,1)), OPENtemp, HIGHtemp, LOWtemp, CLOSEtemp, VOLtemp];
Daily_Data = sortrows(Daily_Data,1);