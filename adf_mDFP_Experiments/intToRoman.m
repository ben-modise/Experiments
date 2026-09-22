function r = intToRoman(x)
% Validates and converts an integer to a Roman numeral
if x < 1 || x > 3999
    error('Roman numerals traditionally only support numbers from 1 to 3999.');
end

% Define mapping arrays
vals = [1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1];
syms = {'M', 'CM', 'D', 'CD', 'C', 'XC', 'L', 'XL', 'X', 'IX', 'V', 'IV', 'I'};

r = ''; 
for i = 1:numel(vals)
    while x >= vals(i)
        r = [r, syms{i}]; %#ok<AGROW>
        x = x - vals(i);
    end
end
end