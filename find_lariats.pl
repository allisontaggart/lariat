#find_lariats.pl

#use warnings;
use Getopt::Std;

#use lib ('/home/ataggart/lariat_scripts');
use lib ('.');
use aligner; use splicemap; use filter; use analyzer;


my %options=();
getopts("hf:i:m:l:o:",\%options);

if (defined $options{h}){
    &displayHelp();
}

elsif ((defined $options{f})&&(defined $options{i})&&(defined $options{o})){

    $file = $options{f};
    $index = $options{i};
    $outdir = $options{o};
    
    if (defined $options{m}){$minfraglen = $options{m};}
    else {$minfraglen = 8;}
    
    if (defined $options{l}){$readlen = $options{l};}
    else {$readlen = 76;}

    ##Align fragments to index
    
    $aligner = aligner->new($file,$index,$readlen,$minfraglen,$outdir);
    $aligner->align();
    
    ##Filter reads
    
    $filter = filter->new($file,$outdir);
    $filter->outOfOrder();
    
    ##Map alignments to splicemap
    
    $splicemap = splicemap->new($outdir,$index);
    
    #1-outoforder, both map
    $splicemap->mapSS("outoforder.txt");
    $splicemap->sameTranscript("outoforder_ss.txt");
    $splicemap->findLariats("outoforder_ss_filter.txt");
    $splicemap->resolveGaps("outoforder_ss_filter_lariats_gap.txt");
   
    
    
    #2 - halfmap
    
    ##Remap half-mapped alignments.  Build bowtie indices of specific transcripts where one half is near ish a splice site
    
    
    ##Format data into bed files of bp/read alignments, and fasta files of bpseq
    
    $analyzer = analyzer->new($outdir,$index,"outoforder_ss_filter_lariats.txt","outoforder_ss_filter_lariats_overlap.txt","outoforder_ss_filter_lariats_gap_truelariats.txt");
    $analyzer->exact();
    $analyzer->overlap();
    $analyzer->gap();
    
}

else {
    &displayHelp();
}

sub displayHelp {
    print "\n\nusage:\nperl find_lariats.pl -f file.fastq -i bowtieindex -o outputdirectory\n\n";
    print "optional flags:\n-h this help message\n-m mininmum fragment length\n-l read lengths\n\n";
}

