#analyzer.pm

package analyzer;
#use warnings;


sub new
{
    my $class = shift;
    my $self = {
        _outdir => shift,
        _index => shift,
        _exact => shift,
        _overlap => shift,
        _gap => shift,
    };
    
    $rm = "rm ".$self->{_outdir}."/lariat_data_table.txt";
    system($rm);
    
    bless $self, $class;
    return $self;  
}

sub exact
{
    my $self = shift;
    $file = $self->{_outdir}."/".$self->{_exact};
    $fileout = $self->{_outdir}."/lariat_data_table.txt";
    open (FH, $file) or die $!;
    open (FHout, ">>$fileout") or die $!;
    while (<FH>){
        $line = $_;
        $id = (split /\t/,$line)[0];
        $seq = (split /\t/,$line)[1];
        $readstrand = (split /\t/,$line)[3];
        $chrom = (split /\t/,$line)[4];
        $headseq = (split /\t/,$line)[6];
        $tailseq = (split /\t/,$line)[14];
        $headcoord1 = (split /\t/,$line)[5];
        $headcoord2 = length($headseq)+$headcoord1;
        $tailcoord1 = (split /\t/,$line)[13];
        $tailcoord2 = length($tailseq)+$tailcoord1;
        my @ss;
        for (my $j = 18; $j<=25; $j++){
                push(@ss,((split /\t/,$line)[$j]));
            }
        for (my $j = 0; $j<=8; $j++){
            if ($ss[$j] =~ /\w+:(\d+)_([+-])/){
                $coord = $1;
                $strand = $2;
                $ss_coord[$j][0] = $coord;
                $ss_coord[$j][1] = $strand;
            }
            else {
                $ss_coord[$j][0] = undef;
                $ss_coord[$j][1] = undef;
            }
        }
        #sense, positive strand
        if (($ss_coord[2][0]==$tailcoord1)&&($ss_coord[2][1] eq $readstrand)&&($readstrand eq "+")&&($ss_coord[2][1] ne undef)){
            $threeprss = $ss_coord[5][0]; 
            print FHout $self->{_outdir},"\texact\t",$id,"\t",$seq,"\t",$chrom,"\t+\t",$tailcoord1,"\t",$threeprss,"\t",$headcoord2,"\t";
            $bpseq1 = $headcoord2-5;
            $bpseq2 = $headcoord2+5;
            $fhbed = $self->{_outdir}."/tempcoord.bed";
            $fhfasta = $self->{_outdir}."/tempseq.fasta";
            open (FHbed, ">$fhbed") or die $!;
            print FHbed $chrom,"\t",$bpseq1,"\t",$bpseq2,"\t",$id,"\t0\t+\n";
            close FHbed;
            $command = "bedtools getfasta -s -fi ".$self->{_index}.".fa -bed ".$fhbed." -fo ".$fhfasta;
            system($command);
            open (FHfasta, $fhfasta) or die $!;
            while (<FHfasta>){
                $_ =~s/\s//g;
                if ($_ =~ /chr/){}
                elsif ($_ =~ /^([ACGTUNacgtun]+)$/){
                    $seq = $1;
                    $seq =~ tr/a-z/A-Z/;
                    print FHout $seq;
                }
            }
             print FHout "\n";
            close FHfasta;
            $command1 = "rm ".$self->{_outdir}."/tempcoord.bed";
            $command2 = "rm ".$self->{_outdir}."/tempseq.fasta";
            system($command1);
            system($command2);     
        }
        #sense, negative strand
        if (($ss_coord[3][0]==$tailcoord2)&&($ss_coord[3][1] eq $readstrand)&&($readstrand eq "-")&&($ss_coord[3][1] ne undef)){
            $threeprss = $ss_coord[4][0]; 
            print FHout $self->{_outdir},"\texact\t",$id,"\t",$seq,"\t",$chrom,"\t-\t",$tailcoord2,"\t",$threeprss,"\t",$headcoord1,"\t";
            $bpseq1 = $headcoord1-5;
            $bpseq2 = $headcoord1+5;
            $fhbed = $self->{_outdir}."/tempcoord.bed";
            $fhfasta = $self->{_outdir}."/tempseq.fasta";
            open (FHbed, ">$fhbed") or die $!;
            print FHbed $chrom,"\t",$bpseq1,"\t",$bpseq2,"\t",$id,"\t0\t-\n";
            close FHbed;
            $command = "bedtools getfasta -s -fi ".$self->{_index}.".fa -bed ".$fhbed." -fo ".$fhfasta;
            system($command);
            open (FHfasta, $fhfasta) or die $!;
            while (<FHfasta>){
                $_ =~s/\s//g;
                if ($_ =~ /chr/){}
                elsif ($_ =~ /^([ACGTUNacgtun]+)$/){
                    $seq = $1;
                    $seq =~ tr/a-z/A-Z/;
                    print FHout $seq;
                }
            }
             print FHout "\n";
            close FHfasta;
            $command1 = "rm ".$self->{_outdir}."/tempcoord.bed";
            $command2 = "rm ".$self->{_outdir}."/tempseq.fasta";
            system($command1);
            system($command2);
        }
        #antisense, positive strand
        if (($ss_coord[0][0]==$headcoord1)&&($ss_coord[0][1] ne $readstrand)&&($readstrand eq "-")&&($ss_coord[0][1] ne undef)){
            $threeprss = $ss_coord[7][0]; 
            print FHout $self->{_outdir},"\texact\t",$id,"\t",$seq,"\t",$chrom,"\t+\t",$headcoord1,"\t",$threeprss,"\t",$tailcoord2,"\t";
            $bpseq1 = $tailcoord2-5;
            $bpseq2 = $tailcoord2+5;
            $fhbed = $self->{_outdir}."/tempcoord.bed";
            $fhfasta = $self->{_outdir}."/tempseq.fasta";
            open (FHbed, ">$fhbed") or die $!;
            print FHbed $chrom,"\t",$bpseq1,"\t",$bpseq2,"\t",$id,"\t0\t+\n";
            close FHbed;
            $command = "bedtools getfasta -s -fi ".$self->{_index}.".fa -bed ".$fhbed." -fo ".$fhfasta;
            system($command);
            open (FHfasta, $fhfasta) or die $!;
            while (<FHfasta>){
                $_ =~s/\s//g;
                if ($_ =~ /chr/){}
                elsif ($_ =~ /^([ACGTUNacgtun]+)$/){
                    $seq = $1;
                    $seq =~ tr/a-z/A-Z/;
                    print FHout $seq;
                }
            }
             print FHout "\n";
            close FHfasta;
            $command1 = "rm ".$self->{_outdir}."/tempcoord.bed";
            $command2 = "rm ".$self->{_outdir}."/tempseq.fasta";
            system($command1);
            system($command2);   
        }
        #antisense, negative strand
        if (($ss_coord[1][0]==$headcoord2)&&($ss_coord[1][1] ne $readstrand)&&($readstrand eq "+")&&($ss_coord[1][1] ne undef)){
            $threeprss = $ss_coord[6][0]; 
            print FHout $self->{_outdir},"\texact\t",$id,"\t",$seq,"\t",$chrom,"\t-\t",$headcoord2,"\t",$threeprss,"\t",$tailcoord1,"\t";
            $bpseq1 = $tailcoord1-5;
            $bpseq2 = $tailcoord1+5;
            $fhbed = $self->{_outdir}."/tempcoord.bed";
            $fhfasta = $self->{_outdir}."/tempseq.fasta";
            open (FHbed, ">$fhbed") or die $!;
            print FHbed $chrom,"\t",$bpseq1,"\t",$bpseq2,"\t",$id,"\t0\t-\n";
            close FHbed;
            $command = "bedtools getfasta -s -fi ".$self->{_index}.".fa -bed ".$fhbed." -fo ".$fhfasta;
            system($command);
            open (FHfasta, $fhfasta) or die $!;
            while (<FHfasta>){
                $_ =~s/\s//g;
                if ($_ =~ /chr/){}
                elsif ($_ =~ /^([ACGTUNacgtun]+)$/){
                    $seq = $1;
                    $seq =~ tr/a-z/A-Z/;
                    print FHout $seq;
                }
            }
             print FHout "\n";
            close FHfasta;
            $command1 = "rm ".$self->{_outdir}."/tempcoord.bed";
            $command2 = "rm ".$self->{_outdir}."/tempseq.fasta";
            system($command1);
            system($command2);
        }
    }
    close FHout; close FH;
} 

sub overlap {
    $self = shift;
    $file = $self->{_outdir}."/".$self->{_overlap};
    $fileout = $self->{_outdir}."/lariat_data_table.txt";
    open (FH, $file) or die $!;
    open (FHout, ">>$fileout") or die $!;
    while (<FH>){
        $line = $_;
        $id = (split /\t/,$line)[0];
        $seq = (split /\t/,$line)[1];
        $readstrand = (split /\t/,$line)[3];
        $chrom = (split /\t/,$line)[4];
        $headseq = (split /\t/,$line)[6];
        $tailseq = (split /\t/,$line)[14];
        $headcoord1 = (split /\t/,$line)[5];
        $headcoord2 = length($headseq)+$headcoord1;
        $tailcoord1 = (split /\t/,$line)[13];
        $tailcoord2 = length($tailseq)+$tailcoord1;
        $overlap = (length($headseq)+length($tailseq))-length($seq);
        

        
        
        my @ss;
        for (my $j = 18; $j<=25; $j++){
                push(@ss,((split /\t/,$line)[$j]));
            }
        for (my $j = 0; $j<=8; $j++){
            if ($ss[$j] =~ /\w+:(\d+)_([+-])/){
                $coord = $1;
                $strand = $2;
                $ss_coord[$j][0] = $coord;
                $ss_coord[$j][1] = $strand;
            }
            else {
                $ss_coord[$j][0] = undef;
                $ss_coord[$j][1] = undef;
            }
        }
        #sense, positive strand
        #if ((abs($ss_coord[2][0]-$tailcoord1)<=$overlap)&&($ss_coord[2][1] eq $readstrand)&&($readstrand eq "+")&&($ss_coord[2][1] ne undef)){
        if ((abs($ss_coord[2][0]-$tailcoord1)<=$overlap)&&($ss_coord[2][0]>=$tailcoord1)&&($ss_coord[2][1] eq $readstrand)&&($readstrand eq "+")&&($ss_coord[2][1] ne undef)){ 
            $threeprss = $ss_coord[5][0];
            $excess = $ss_coord[2][0] - $tailcoord1;
            
            $tailcoord1 = $ss_coord[2][0];
            $headcoord2 = $headcoord2-($overlap-$excess); 
            print FHout $self->{_outdir},"\toverlap ",$overlap,"\t",$id,"\t",$seq,"\t",$chrom,"\t+\t",$tailcoord1,"\t",$threeprss,"\t",$headcoord2,"\t";
            $bpseq1 = $headcoord2-5;
            $bpseq2 = $headcoord2+5;
            $fhbed = $self->{_outdir}."/tempcoord.bed";
            $fhfasta = $self->{_outdir}."/tempseq.fasta";
            open (FHbed, ">$fhbed") or die $!;
            print FHbed $chrom,"\t",$bpseq1,"\t",$bpseq2,"\t",$id,"\t0\t+\n";
            close FHbed;
            $command = "bedtools getfasta -s -fi ".$self->{_index}.".fa -bed ".$fhbed." -fo ".$fhfasta;
            system($command);
            open (FHfasta, $fhfasta) or die $!;
            while (<FHfasta>){
                $_ =~s/\s//g;
                if ($_ =~ /chr/){}
                elsif ($_ =~ /^([ACGTUNacgtun]+)$/){
                    $seq = $1;
                    $seq =~ tr/a-z/A-Z/;
                    print FHout $seq;
                }
            }
             print FHout "\n";
            close FHfasta;
            $command1 = "rm ".$self->{_outdir}."/tempcoord.bed";
            $command2 = "rm ".$self->{_outdir}."/tempseq.fasta";
            system($command1);
            system($command2);     
        }
        #sense, negative strand
        if (((-1)*($ss_coord[3][0]-$tailcoord2)<=$overlap)&&($tailcoord2>=$ss_coord[3][0])&&($ss_coord[3][1] eq $readstrand)&&($readstrand eq "-")&&($ss_coord[3][1] ne undef)){
        
            $threeprss = $ss_coord[4][0]; 
            $excess = $tailcoord2-$ss_coord[3][0];
            $tailcoord2 = $ss_coord[3][0];
            $headcoord1 = $headcoord1+($overlap-$excess);
            print FHout $self->{_outdir},"\toverlap ",$overlap,"\t",$id,"\t",$seq,"\t",$chrom,"\t-\t",$tailcoord2,"\t",$threeprss,"\t",$headcoord1,"\t";
            $bpseq1 = $headcoord1-5;
            $bpseq2 = $headcoord1+5;
            $fhbed = $self->{_outdir}."/tempcoord.bed";
            $fhfasta = $self->{_outdir}."/tempseq.fasta";
            open (FHbed, ">$fhbed") or die $!;
            print FHbed $chrom,"\t",$bpseq1,"\t",$bpseq2,"\t",$id,"\t0\t-\n";
            close FHbed;
            $command = "bedtools getfasta -s -fi ".$self->{_index}.".fa -bed ".$fhbed." -fo ".$fhfasta;
            system($command);
            open (FHfasta, $fhfasta) or die $!;
            while (<FHfasta>){
                $_ =~s/\s//g;
                if ($_ =~ /chr/){}
                elsif ($_ =~ /^([ACGTUNacgtun]+)$/){
                    $seq = $1;
                    $seq =~ tr/a-z/A-Z/;
                    print FHout $seq;
                }
            }
             print FHout "\n";
            close FHfasta;
            $command1 = "rm ".$self->{_outdir}."/tempcoord.bed";
            $command2 = "rm ".$self->{_outdir}."/tempseq.fasta";
            system($command1);
            system($command2);
        }
        #antisense, positive strand
        if ((($ss_coord[0][0]-$headcoord1)<=$overlap)&&($ss_coord[0][0] >=$headcoord1)&&($ss_coord[0][1] ne $readstrand)&&($readstrand eq "-")&&($ss_coord[0][1] ne undef)){
            $threeprss = $ss_coord[7][0];
            $excess = $ss_coord[0][0]-$headcoord1;
            $headcoord1 = $ss_coord[0][0];
            $tailcoord2 = $tailcoord2-($overlap-$excess); 
            print FHout $self->{_outdir},"\toverlap ",$overlap,"\t",$id,"\t",$seq,"\t",$chrom,"\t+\t",$headcoord1,"\t",$threeprss,"\t",$tailcoord2,"\t";
            $bpseq1 = $tailcoord2-5;
            $bpseq2 = $tailcoord2+5;
            $fhbed = $self->{_outdir}."/tempcoord.bed";
            $fhfasta = $self->{_outdir}."/tempseq.fasta";
            open (FHbed, ">$fhbed") or die $!;
            print FHbed $chrom,"\t",$bpseq1,"\t",$bpseq2,"\t",$id,"\t0\t+\n";
            close FHbed;
            $command = "bedtools getfasta -s -fi ".$self->{_index}.".fa -bed ".$fhbed." -fo ".$fhfasta;
            system($command);
            open (FHfasta, $fhfasta) or die $!;
            while (<FHfasta>){
                $_ =~s/\s//g;
                if ($_ =~ /chr/){}
                elsif ($_ =~ /^([ACGTUNacgtun]+)$/){
                    $seq = $1;
                    $seq =~ tr/a-z/A-Z/;
                    print FHout $seq;
                }
            }
             print FHout "\n";
            close FHfasta;
            $command1 = "rm ".$self->{_outdir}."/tempcoord.bed";
            $command2 = "rm ".$self->{_outdir}."/tempseq.fasta";
            system($command1);
            system($command2);   
        }
        #antisense, negative strand
        if (((-1)*($ss_coord[1][0]-$headcoord2)<=$overlap)&&($headcoord2>=$ss_coord[1][0])&&($ss_coord[1][1] ne $readstrand)&&($readstrand eq "+")&&($ss_coord[1][1] ne undef)){
            $threeprss = $ss_coord[6][0];
            $excess = $headcoord2-$ss_coord[1][0];
            $headcoord2 = $ss_coord[1][0];
            $tailcoord1 = $tailcoord1+($overlap-$excess); 
            print FHout $self->{_outdir},"\toverlap ",$overlap,"\t",$id,"\t",$seq,"\t",$chrom,"\t-\t",$headcoord2,"\t",$threeprss,"\t",$tailcoord1,"\t";
            $bpseq1 = $tailcoord1-5;
            $bpseq2 = $tailcoord1+5;
            $fhbed = $self->{_outdir}."/tempcoord.bed";
            $fhfasta = $self->{_outdir}."/tempseq.fasta";
            open (FHbed, ">$fhbed") or die $!;
            print FHbed $chrom,"\t",$bpseq1,"\t",$bpseq2,"\t",$id,"\t0\t-\n";
            close FHbed;
            $command = "bedtools getfasta -s -fi ".$self->{_index}.".fa -bed ".$fhbed." -fo ".$fhfasta;
            system($command);
            open (FHfasta, $fhfasta) or die $!;
            while (<FHfasta>){
                $_ =~s/\s//g;
                if ($_ =~ /chr/){}
                elsif ($_ =~ /^([ACGTUNacgtun]+)$/){
                    $seq = $1;
                    $seq =~ tr/a-z/A-Z/;
                    print FHout $seq;
                }
            }
             print FHout "\n";
            close FHfasta;
            $command1 = "rm ".$self->{_outdir}."/tempcoord.bed";
            $command2 = "rm ".$self->{_outdir}."/tempseq.fasta";
            system($command1);
            system($command2);
        }
    }
    close FHout; close FH;
    
    
    
}

sub gap {
    $self = shift;
    $file = $self->{_outdir}."/".$self->{_gap};
    $fileout = $self->{_outdir}."/lariat_data_table.txt";
    open (FH, $file) or die $!;
    open (FHout, ">>$fileout") or die $!;
    while (<FH>){
        $line = $_;
        $id = (split /\t/,$line)[0];
        $seq = (split /\t/,$line)[1];
        $readstrand = (split /\t/,$line)[3];
        $chrom = (split /\t/,$line)[4];
        $headseq = (split /\t/,$line)[6];
        $tailseq = (split /\t/,$line)[14];
        $headcoord1 = (split /\t/,$line)[5];
        $headcoord2 = length($headseq)+$headcoord1;
        $tailcoord1 = (split /\t/,$line)[13];
        $tailcoord2 = length($tailseq)+$tailcoord1;
        $gap = length($seq)-(length($headseq)+length($tailseq));
        my @ss;
        for (my $j = 18; $j<=25; $j++){
                push(@ss,((split /\t/,$line)[$j]));
            }
        for (my $j = 0; $j<=8; $j++){
            if ($ss[$j] =~ /\w+:(\d+)_([+-])/){
                $coord = $1;
                $strand = $2;
                $ss_coord[$j][0] = $coord;
                $ss_coord[$j][1] = $strand;
            }
            else {
                $ss_coord[$j][0] = undef;
                $ss_coord[$j][1] = undef;
            }
        }
        #sense, positive strand
        if ((abs($ss_coord[2][0]-$tailcoord1)<=$gap)&&($ss_coord[2][1] eq $readstrand)&&($readstrand eq "+")&&($ss_coord[2][1] ne undef)){
            $threeprss = $ss_coord[5][0];
            $remainder = $tailcoord1-$ss_coord[2][0];
            $tailcoord1 = $ss_coord[2][0];
            $headcoord2 = $headcoord2+($gap-$remainder); 
            print FHout $self->{_outdir},"\tgap ",$gap,"\t",$id,"\t",$seq,"\t",$chrom,"\t+\t",$tailcoord1,"\t",$threeprss,"\t",$headcoord2,"\t";
            $bpseq1 = $headcoord2-5;
            $bpseq2 = $headcoord2+5;
            $fhbed = $self->{_outdir}."/tempcoord.bed";
            $fhfasta = $self->{_outdir}."/tempseq.fasta";
            open (FHbed, ">$fhbed") or die $!;
            print FHbed $chrom,"\t",$bpseq1,"\t",$bpseq2,"\t",$id,"\t0\t+\n";
            close FHbed;
            $command = "bedtools getfasta -s -fi ".$self->{_index}.".fa -bed ".$fhbed." -fo ".$fhfasta;
            system($command);
            open (FHfasta, $fhfasta) or die $!;
            while (<FHfasta>){
                $_ =~s/\s//g;
                if ($_ =~ /chr/){}
                elsif ($_ =~ /^([ACGTUNacgtun]+)$/){
                    $seq = $1;
                    $seq =~ tr/a-z/A-Z/;
                    print FHout $seq;
                }
            }
            print FHout "\n";
            close FHfasta;
            $command1 = "rm ".$self->{_outdir}."/tempcoord.bed";
            $command2 = "rm ".$self->{_outdir}."/tempseq.fasta";
            system($command1);
            system($command2);     
        }
        #sense, negative strand
        if ((abs($ss_coord[3][0]-$tailcoord2)<=$gap)&&($ss_coord[3][1] eq $readstrand)&&($readstrand eq "-")&&($ss_coord[3][1] ne undef)){
            $threeprss = $ss_coord[4][0]; 
            $remainder = $ss_coord[3][0]-$tailcoord2;
            $tailcoord2 = $ss_coord[3][0];
            $headcoord1 = $headcoord1-($gap-$remainder);
            print FHout $self->{_outdir},"\tgap ",$gap,"\t",$id,"\t",$seq,"\t",$chrom,"\t-\t",$tailcoord2,"\t",$threeprss,"\t",$headcoord1,"\t";
            $bpseq1 = $headcoord1-5;
            $bpseq2 = $headcoord1+5;
            $fhbed = $self->{_outdir}."/tempcoord.bed";
            $fhfasta = $self->{_outdir}."/tempseq.fasta";
            open (FHbed, ">$fhbed") or die $!;
            print FHbed $chrom,"\t",$bpseq1,"\t",$bpseq2,"\t",$id,"\t0\t-\n";
            close FHbed;
            $command = "bedtools getfasta -s -fi ".$self->{_index}.".fa -bed ".$fhbed." -fo ".$fhfasta;
            system($command);
            open (FHfasta, $fhfasta) or die $!;
            while (<FHfasta>){
                $_ =~s/\s//g;
                if ($_ =~ /chr/){}
                elsif ($_ =~ /^([ACGTUNacgtun]+)$/){
                    $seq = $1;
                    $seq =~ tr/a-z/A-Z/;
                    print FHout $seq;
                }
            }
            print FHout "\n";
            close FHfasta;
            $command1 = "rm ".$self->{_outdir}."/tempcoord.bed";
            $command2 = "rm ".$self->{_outdir}."/tempseq.fasta";
            system($command1);
            system($command2);
        }
        #antisense, positive strand
        if ((abs($ss_coord[0][0]-$headcoord1)<=$gap)&&($ss_coord[0][1] ne $readstrand)&&($readstrand eq "-")&&($ss_coord[0][1] ne undef)){
            $threeprss = $ss_coord[7][0];
            $remainder = $headcoord1-$ss_coord[0][0];
            $headcoord1 = $ss_coord[0][0];
            $tailcoord2 = $tailcoord2+($gap-$remainder); 
            print FHout $self->{_outdir},"\tgap ",$gap,"\t",$id,"\t",$seq,"\t",$chrom,"\t+\t",$headcoord1,"\t",$threeprss,"\t",$tailcoord2,"\t";
            $bpseq1 = $tailcoord2-5;
            $bpseq2 = $tailcoord2+5;
            $fhbed = $self->{_outdir}."/tempcoord.bed";
            $fhfasta = $self->{_outdir}."/tempseq.fasta";
            open (FHbed, ">$fhbed") or die $!;
            print FHbed $chrom,"\t",$bpseq1,"\t",$bpseq2,"\t",$id,"\t0\t+\n";
            close FHbed;
            $command = "bedtools getfasta -s -fi ".$self->{_index}.".fa -bed ".$fhbed." -fo ".$fhfasta;
            system($command);
            open (FHfasta, $fhfasta) or die $!;
            while (<FHfasta>){
                $_ =~s/\s//g;
                if ($_ =~ /chr/){}
                elsif ($_ =~ /^([ACGTUNacgtun]+)$/){
                    $seq = $1;
                    $seq =~ tr/a-z/A-Z/;
                    print FHout $seq;
                }
            }
             print FHout "\n";
            close FHfasta;
            $command1 = "rm ".$self->{_outdir}."/tempcoord.bed";
            $command2 = "rm ".$self->{_outdir}."/tempseq.fasta";
            system($command1);
            system($command2);   
        }
        #antisense, negative strand
        if ((abs($ss_coord[1][0]-$headcoord2)<=$gap)&&($ss_coord[1][1] ne $readstrand)&&($readstrand eq "+")&&($ss_coord[1][1] ne undef)){
            $threeprss = $ss_coord[6][0];
            $remainder = $ss_coord[1][0]-$headcoord2;
            $headcoord2 = $ss_coord[1][0];
            $tailcoord1 = $tailcoord1-($gap-$remainder); 
            print FHout $self->{_outdir},"\tgap ",$gap,"\t",$id,"\t",$seq,"\t",$chrom,"\t-\t",$headcoord2,"\t",$threeprss,"\t",$tailcoord1,"\t";
            $bpseq1 = $tailcoord1-5;
            $bpseq2 = $tailcoord1+5;
            $fhbed = $self->{_outdir}."/tempcoord.bed";
            $fhfasta = $self->{_outdir}."/tempseq.fasta";
            open (FHbed, ">$fhbed") or die $!;
            print FHbed $chrom,"\t",$bpseq1,"\t",$bpseq2,"\t",$id,"\t0\t-\n";
            close FHbed;
            $command = "bedtools getfasta -s -fi ".$self->{_index}.".fa -bed ".$fhbed." -fo ".$fhfasta;
            system($command);
            open (FHfasta, $fhfasta) or die $!;
            while (<FHfasta>){
		$_ =~s/\s//g;
                if ($_ =~ /chr/){}
                elsif ($_ =~ /^([ACGTUNacgtun]+)$/){
                    $seq = $1;
                    $seq =~ tr/a-z/A-Z/;
                    print FHout $seq;
                }
            }
             print FHout "\n";
            close FHfasta;
            $command1 = "rm ".$self->{_outdir}."/tempcoord.bed";
            $command2 = "rm ".$self->{_outdir}."/tempseq.fasta";
            system($command1);
            system($command2);
        }
    }
    close FHout; close FH;
}     
        
1;
