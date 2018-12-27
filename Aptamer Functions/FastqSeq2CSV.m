function list = FastqSeq2CSV(input_fastq_name, output_csv_name);
%1st input: fastq file name 'xxx.fastq'
%2nd input : csv file name you want to get in 'xxxx.csv'
%output: a written csv and seqs (very long), use';'

[headers, seqs, quals] = fastqread(input_fastq_name);
list = string(seqs);

%%% write cell to csv line by line(string cells)
fid = fopen(output_csv_name, 'w');
fprintf(fid,'%s \n', list(1:end));
fclose(fid);