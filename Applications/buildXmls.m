clear
% get fasta files in fasta
fastafiles = dir("fasta/*.fasta");

percentages = [0.1 0.25 0.5 1];
template = {'Default', 'Targeted', 'Intervals', 'Topology'};

cl = [10 20 50 100]*10^5;
ll = [10 20 50 100]*10^2;
for i = 1 : length(fastafiles)
    fasta = fastaread(['fasta/' fastafiles(i).name]);
    fname = strsplit(fastafiles(i).name, '.');
    fname = fname{1};
    for j = 1 :length(percentages)
        indices = sort(randsample(length(fasta),round(length(fasta)*percentages(j))));

        for k = 1 :length(template)
            f = fopen(['Template/' template{k} '.xml']);
            g = fopen(['xmls/' fname '_' template{k} '_' num2str(percentages(j)) '.xml'], 'w');
            while ~feof(f)
                line = fgets(f);
                if contains(line, 'insert_sequences')
                    for l =1 : length(indices)
                        fprintf(g, '\t<sequence id="seq_%s" spec="Sequence" taxon="%s" totalcount="4" value="%s"/>\n',...
                            fasta(indices(l)).Header,fasta(indices(l)).Header,fasta(indices(l)).Sequence);
                    end
                elseif contains(line, 'insert_chain_length')
                    fprintf(g, strrep(line, 'insert_chain_length', num2str(cl(j))));
                elseif contains(line, 'insert_logEvery')
                    fprintf(g, strrep(line, 'insert_logEvery', num2str(ll(j))));
                else
                    fprintf(g, line);
                end
            end
            fclose('all');
               
        end
    end
end