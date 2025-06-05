clear


fastafiles = dir('fasta/*.fasta');
for i = 4:length(fastafiles)
    txtfile = strrep(fastafiles(i).name, '.fasta', '.txt');
    % if ~exist(['fasta/' txtfile])
        fasta = fastaread(['fasta/' fastafiles(i).name]);
        f = fopen(['fasta/' txtfile], 'w');
        for j = 1 :length(fasta)
            tmp = strsplit(fasta(j).Header, '|');
            fprintf(f,'PP_%s ', tmp{1});
        end
        fclose('all')
    % end
end

i=2
txtfile = strrep(fastafiles(i).name, '.fasta', '.txt');
fasta = fastaread(['fasta/' fastafiles(i).name])
        f = fopen(['fasta/' txtfile], 'w');

for j = 1 :length(fasta)
    tmp = strsplit(fasta(j).Header, '|')
    fprintf(f,'%s ', tmp{2});
end
        fclose('all')

