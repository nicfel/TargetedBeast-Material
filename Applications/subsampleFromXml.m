clear
f = fopen('xmls_dta/h3n2_recent_Targeted_dta_10000_rep0.xml');
g = fopen('xmls_dta/h3n2_recent_Targeted_dta_1000_rep0.xml', 'w');
c = 1;
rng(1234)

% sample 1000 numbers between 1 and 10000
use_samples = randsample(10000, 1000, false);
remove_id = cell(0,0);
while ~feof(f)
    line = fgets(f);
    if contains(line, '<sequence') 
        tmp = strsplit(line, '"');
        Data(c).Header = tmp{6};
        Data(c).Sequence = tmp{10};
        if ismember(c, use_samples)
            fprintf(g, line);
        else
            remove_id{end+1} = tmp{6};
        end
        c=c+1;
    elseif contains(line, 'id="loc.traitSet.trait')
        c = 1;
        fprintf(g, line);
        line = fgets(f);

        while true
            tmp = strsplit(line, '=');
            tmp = regexprep(tmp{1}, '\t+', '');
            if ismember(tmp, remove_id)
                
            else
                fprintf(g, line);
                
            end
            line = fgets(f);
            if contains(line, '</traitSet>')
                fprintf(g, line);
                break;
            end               
        end

    else
        fprintf(g, line);
    end
end
fclose('all');


