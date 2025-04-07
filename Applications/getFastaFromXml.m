clear
f = fopen('xmls_long/h3n2_recent_Default_5000_rep0.xml');
c = 1;
while ~feof(f)
    line = fgets(f);
    if contains(line, '<sequence')
        tmp = strsplit(line, '"');
        Data(c).Header = tmp{6};
        Data(c).Sequence = tmp{10};
        c=c+1;
    end
end

fastawrite('h3n2_recent_5000.fasta', Data)
