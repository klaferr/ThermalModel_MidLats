function s = readStruct(file)
% READSTRUCT(FILE) Reads structure from text file, FILE.
% The first entries of each row are field names.
% The second entry in each row is the value.
% Anything after that is ignored as comments for that row.

% Open the specified file.
fid = fopen(file);

% Read in the first line and set the line number.
line = fgetl(fid);

% The function, fgetl, returns a -1 when an end-of-file is encountered.
while ~isequal(line, -1);
    
    % Initially, no items from this line
    items = {};
    
    % Calculate number of tokens in the line
    line2 = line;
    numberOfItems = 0;
    while line2
        numberOfItems = numberOfItems + 1;
        [test, line2] = strtok(line2);
    end
    
    % Only lines that don't start with % (those are assumed to be comments)
    if numberOfItems ~= 0 && line(1) ~= '%'
        
        % Only process non-empty lines
        if numberOfItems > 1
            
            % Get first two components
            [items{1}, line] = strtok(line);
            [items{2}, line] = strtok(line);
            
            
            % 2nd item on line is value of field. If first element of last item is a letter,
            % then assume that last item is a char array. Otherwise, assume it is a number.
            if isletter(items{2}(1))
                value = ['''', items{end}, '''']; % Wrap quotes around string.
            else
                value = items{2};
            end
            
            % Remove value from list of items (leaving only field names, the first element of items).
            items = items(1);
            
            fieldname = 's'; % Initialize construction of field name.
            for item = items
                fieldname = [fieldname, '.', item{:}]; % Append subfield names to field name.
            end
            
            % Can uncomment the 'line' here if getting stuck on a bad input to
            % figure out which one
            % line
            
            % Construct command to assign value to field name.
            command = [fieldname, ' = ', value, ';'];
            
            eval(command);
            
        elseif numberOfItems == 1;
            
            error ('Every non-empty line must contain must contain at least one field name and one value.')
            
        end
    end
    
    line = fgetl(fid);
    
end

fclose(fid);