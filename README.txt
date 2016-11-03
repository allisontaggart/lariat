Make sure to have bowtie and bedtools in path.


usage:
perl find_lariats.pl -f file.fastq -i bowtieindex -o outputdirectory

optional flags:
-h this help message
-m mininmum fragment length
-l read lengths


Results are in file "lariat_data_table.txt"
columns are
1-directory
2-inverted alignment type
3-read ID
4-raw read sequence
5-chromosome
6-5'ss
7-3'ss
8-BP
9-raw branch site sequence

indexes included:

human (hg19): -i ./hg19/hg19
mouse (mm9): -i ./mm9/mm9
s. pombe (EF2): -i ./EF2/genome
