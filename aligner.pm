#aligner.pm

package aligner;

use warnings;

sub new
{
    my $class = shift;
    my $self = {
        _rawfile => shift,
        _index => shift,
        _readlen => shift,
        _minfraglen => shift,
        _outdir => shift,
    };
    
    print "\n\n****************\nparsing fastq files\n\n";
    
    $self->{_rawfile} =~ /(.+).fa?s?t?q/;
    $self->{_base} = $1;
    $self->{_file}=$self->{_outdir}."/seq.fastq";
    
    $mkdir = "mkdir ".$self->{_outdir};
    $cp = "cp ".$self->{_rawfile}." ".$self->{_file};
    
    system($mkdir);
    system($cp);
    
    bless $self, $class;
    return $self;
    
}

sub forwardAlign {
    
    my $self = shift;
    my $file = $self->{_file};
    my $index = $self->{_index};
    
    print "\n\n****************\nconducting forward alignment: \n";
    
    $command = "bowtie -v3 -p8 -k1 ".$index." ".$file." --un ".$self->{_outdir}."/unaligned.fastq > ".$self->{_outdir}."/aligned.txt";
    print $command,"\n\n";
    system($command);

    $cp1 = "cp ".$self->{_outdir}."/unaligned.fastq ".$self->{_outdir}."/unaligned_left.fastq";
    $cp2 = "cp ".$self->{_outdir}."/unaligned.fastq ".$self->{_outdir}."/unaligned_right.fastq";
    
    system($cp1);
    system($cp2);
    
}

sub fragmentAlignLeft {
    
    my $self = shift;
    
    print "\n\n****************\nprocessing left side: \n";
    
    my $file = $self->{_outdir}."/unaligned_left.fastq";
    my $index = $self->{_index};
    
    for (my $j = $self->{_minfraglen}; $j<= ($self->{_readlen}-$self->{_minfraglen}); $j++){
        $file =~ /(.*left)/;
        $remainderfile = $1."_".$j.".fastq";
        $maxfile = $1."_multalignments_".$j.".fastq";
        $outfile = $1."_".$j.".txt";
        $alignedfq = $self->{_outdir}."/left_aligned_".$j.".fastq";
        $command = "bowtie -v0 -p8 -a -m1 --trim3 ".$j." ".$index." ".$file." --un ".$remainderfile." --max ".$maxfile." --al ".$alignedfq." > ".$outfile;
        print "\n\n*****\nexecuting: ",$command,"\n";
        system($command);
        
        $prevrem = $file;
        
        $rmmax = "rm ".$maxfile;
        $rmprevrem = "rm ".$prevrem;
        
        system($rmmax);
        system($rmprevrem);
        
        $file = $remainderfile;
    }
    
}

sub fragmentAlignRight {
    
    my $self = shift;
    
    print "\n\n****************\nprocessing right side: \n";
    
    my $file = $self->{_outdir}."/unaligned_right.fastq";
    my $index = $self->{_index};
    
    for (my $j = $self->{_minfraglen}; $j<= ($self->{_readlen}-$self->{_minfraglen}); $j++){
        $file =~ /(.*right)/;
        $remainderfile = $1."_".$j.".fastq";
        $maxfile = $1."_multalignments_".$j.".fastq";
        $outfile = $1."_".$j.".txt";
        $alignedfq = $self->{_outdir}."/right_aligned_".$j.".fastq";
        $command = "bowtie -v0 -p8 -a -m1 --trim5 ".$j." ".$index." ".$file." --un ".$remainderfile." --max ".$maxfile." --al ".$alignedfq." > ".$outfile;
        print "\n\n*****\nexecuting: ",$command,"\n";
        system($command);
        
        $prevrem = $file;
        
        $rmmax = "rm ".$maxfile;
        $rmprevrem = "rm ".$prevrem;
        
        system($rmmax);
        system($rmprevrem);
        
        $file = $remainderfile;
        
    }
    
}


sub align {
    
    my $self = shift;
    
    &forwardAlign($self);
    &fragmentAlignLeft($self);
    &fragmentAlignRight($self);
    &mergeSides($self);
    &rmTempFiles($self);
    
}


sub mergeSides {
    
    my $self = shift;
    
    print "\n\n****************\ncompiling left and right sides: \n";
    
    $alldatafile = $self->{_outdir}."/new_alignments.txt";
    open (FHout, ">$alldatafile") or die $!;
    
    my %lefthits;
    
    for (my $j = $self->{_minfraglen}; $j<= ($self->{_readlen}-$self->{_minfraglen}); $j++){
        $leftsam = $self->{_outdir}."/unaligned_left_".$j.".txt";
        open (FH, "$leftsam");
        while (<FH>){
            $data = $_;
            chomp $data;
            $id = (split /\t/,$data)[0];
            if ($id =~ /\@/){}
            else {
                $id = (split /\t/,$data)[0];
                $lefthits{$id}[0] = 0;
                $lefthits{$id}[1] = $data;
                $lefthits{$id}[2] = "0\t0\t0\t0\t0\t0\t0\t0";
            }
        }
        close FH;
    }

    for (my $j = $self->{_minfraglen}; $j<= ($self->{_readlen}-$self->{_minfraglen}); $j++){
        $leftfq = $self->{_outdir}."/left_aligned_".$j.".fastq";
        if (-e $leftfq){
        open (FH,"$leftfq");
        $line = 1;
        while (<FH>){
            chomp $_;
            if ($line ==1){
                $id = $_;
                $id =~ s/^@//g;
            }
            elsif (($line==2)&&($_ =~ /([ACGTNacgtn]+)/)){
                $seq = $1;
                $lefthits{$id}[0] = $seq;
            }
            $line++;
            if ($line==5){
                $line = 1;
            }
        }}
    }
    
    for (my $j = $self->{_minfraglen}; $j<= ($self->{_readlen}-$self->{_minfraglen}); $j++){
        $rightsam = $self->{_outdir}."/unaligned_right_".$j.".txt";
        open (FH, "$rightsam");
        while (<FH>){
            $data = $_;
            chomp $data;
            $id = (split /\t/,$data)[0];
            if ($id =~ /\@/){}
            else {
                $id = (split /\t/,$data)[0];
                if (exists $lefthits{$id}[0]){
                    $lefthits{$id}[2] = $data;
                }
                else {
                    $lefthits{$id}[0] = 0;
                    $lefthits{$id}[2] = $data;
                    $lefthits{$id}[1] = "0\t0\t0\t0\t0\t0\t0\t0";
                }
            }
        }
        close FH;
    }

    for (my $j = $self->{_minfraglen}; $j<= ($self->{_readlen}-$self->{_minfraglen}); $j++){
        $rightfq = $self->{_outdir}."/right_aligned_".$j.".fastq";
        if (-e $rightfq){
        open (FH,"$rightfq");
        $line = 1;
        while (<FH>){
            chomp $_;
            if ($line ==1){
                $id = $_;
                $id =~ s/^@//g;
            }
            elsif (($line==2)&&($_ =~ /([ACGTNacgtn]+)/)){
                $seq = $1;
                $lefthits{$id}[0] = $seq;
            }
            $line++;
            if ($line==5){
                $line = 1;
            }
        }}
    }
    
    foreach my $key (keys %lefthits){
        print FHout $key,"\t";
        print FHout $lefthits{$key}[0],"\t";
        print FHout $lefthits{$key}[1],"\t";
        print FHout $lefthits{$key}[2],"\n";
    }
    close FHout;
}

sub rmTempFiles {
    
    my $self = shift;
    
    $rm = "rm ".$self->{_outdir}."/*.fastq ".$self->{_outdir}."/unaligned*.txt";
    print "\n",$rm,"\n\n";
    system($rm);
    
}

1;