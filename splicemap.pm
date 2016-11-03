#splicemap.pm

package splicemap;



#use warnings;

sub new 
{
    my $class = shift;
    my $self = {
        _outdir => shift,
        _index => shift,
    };
    
    my $fh = $self->{_index}."_ss_table.txt";
    
    print "\n\n****************\nbuilding splice site hash map\n\n";
    
    my %splice5pr=();
    my %splice3pr=();
    

    open (FH,"$fh") or die $!;
    while (<FH>){

        $chr = (split /\t/,$_)[1];
        $id = (split /\t/,$_)[0];
        $strand = (split /\t/,$_)[2];
        $exoncount = (split /\t/,$_)[5];
        
        if ($exoncount =~ /\d+/){
            $exonstarts = (split /\t/,$_)[6];
            $exonstops = (split /\t/,$_)[7];
            @starts = (split /\,/,$exonstarts);
            @stops = (split /\,/,$exonstops);
            chomp $exonstops;
        
            for (my $j = 1; $j<$exoncount; $j++){
                if ($strand eq "+"){
                    $fiveprid = $chr.":".$stops[$j-1]."_".$strand;
                    $threeprid = $chr.":".$starts[$j]."_".$strand;
                }
                else {
                    $fiveprid = $chr.":".$starts[$j]."_".$strand;
                    $threeprid = $chr.":".$stops[$j-1]."_".$strand;
                }
        
                if (exists $splice5pr{$fiveprid}){
                    $value = $splice5pr{$fiveprid};
                    $newvalue = $value.",".$id;
                    $splice5pr{$fiveprid} = $newvalue;
                }
                else {
                    $splice5pr{$fiveprid} = $id;
                }
        
                if (exists $splice3pr{$threeprid}){
                    $value = $splice3pr{$threeprid};
                    $newvalue = $value.",".$id;
                    $splice3pr{$threeprid} = $newvalue;
                }
                else {
                    $splice3pr{$threeprid} = $id;
                }
            }
        }
    }

    $self->{_ss5pr} = \%splice5pr;
    $self->{_ss3pr} = \%splice3pr;
    
    bless $self, $class;
    return $self;
    
}

sub mapSS
{
    $self = shift;
    $group = shift;
    
    $group =~ /(.+).txt/;
    $base = $1;
    
    $fh = $self->{_outdir}."/".$group;
    $fhout = $self->{_outdir}."/".$base."_ss.txt";
    
    %ss5 = %{$self->{_ss5pr}};
    %ss3 = %{$self->{_ss3pr}};
    
    open (FH, $fh) or die $!;
    open (FHout, ">$fhout") or die $!;
    
    while (<FH>){
        #loop through and lookup splice hash as you go
        $data = $_;
        chomp $data;
        
        $seq = (split /\t/,$data)[1];
        $length = length($seq);
        
        $headchrom = (split /\t/,$data)[4];
        $headstrand = (split /\t/,$data)[3];    
        $headcoord1 = (split /\t/,$data)[5];
        $headseq = (split /\t/,$data)[6];
        $headlen = length($headseq);
        $headcoord2 = $headcoord1+$headlen;

        
        $tailstrand = (split /\t/,$data)[11];
        $tailchrom = (split /\t/,$data)[12];
        $tailcoord1 = (split /\t/,$data)[13];
        $tailseq = (split /\t/,$data)[14];
        $taillen = length($tailseq);
        $tailcoord2 = $tailcoord1+$taillen;


        my @lookup;
        
        print FHout $data,"\t";
        
        $lookupkey[0] = $headchrom.":".$headcoord1."_".$headstrand;
        $lookupkey[1] = $headchrom.":".$headcoord2."_".$headstrand;
        $lookupkey[2] = $tailchrom.":".$tailcoord1."_".$tailstrand;
        $lookupkey[3] = $tailchrom.":".$tailcoord2."_".$tailstrand;


        
        
        for (my $j = 0; $j<=3; $j++){
            $foundflag = 0;
            $dist = 0;
            if ($lookupkey[$j] =~ /(\w+):(\d+)_([+-])/){
            $chrom = $1; $coord = $2; $strand = $3;
            while (($dist<=1000)&&($foundflag ==0)){
                $newcoord = $coord+$dist;
                $lookup = $chrom.":".$newcoord."_+";
                if ((exists $ss5{$lookup})&&($foundflag==0)){
                    print FHout ($lookup.":".$ss5{$lookup}),"\t";
                    $foundflag = 1;
                }
                $newcoord = $coord-$dist;
                $lookup = $chrom.":".$newcoord."_+";
                if ((exists $ss5{$lookup})&&($foundflag==0)){
                    print FHout ($lookup.":".$ss5{$lookup}),"\t";
                    $foundflag = 1;
                }
                $newcoord = $coord+$dist;
                $lookup = $chrom.":".$newcoord."_-";
                if ((exists $ss5{$lookup})&&($foundflag==0)){
                    print FHout ($lookup.":".$ss5{$lookup}),"\t";
                    $foundflag = 1;
                }
                $newcoord = $coord-$dist;
                $lookup = $chrom.":".$newcoord."_-";
                if ((exists $ss5{$lookup})&&($foundflag==0)){
                    print FHout ($lookup.":".$ss5{$lookup}),"\t";
                    $foundflag = 1;
                }
                $dist++;
            }}
            else {print "error: $lookupkey[$j]\n";}
            if ($foundflag ==0){print FHout "null\t";}
        } 
        
        for (my $j = 0; $j<=3; $j++){
            $foundflag = 0;
            $dist = 0;
            if ($lookupkey[$j] =~ /(\w+):(\d+)_([+-])/){
            $chrom = $1; $coord = $2; $strand = $3;
            while (($dist<=1000)&&($foundflag ==0)){
                $newcoord = $coord+$dist;
                $lookup = $chrom.":".$newcoord."_+";
                if ((exists $ss3{$lookup})&&($foundflag==0)){
                    print FHout ($lookup.":".$ss3{$lookup}),"\t";
                    $foundflag = 1;
                }
                $newcoord = $coord-$dist;
                $lookup = $chrom.":".$newcoord."_+";
                if ((exists $ss3{$lookup})&&($foundflag==0)){
                    print FHout ($lookup.":".$ss3{$lookup}),"\t";
                    $foundflag = 1;
                }
                $newcoord = $coord+$dist;
                $lookup = $chrom.":".$newcoord."_-";
                if ((exists $ss3{$lookup})&&($foundflag==0)){
                    print FHout ($lookup.":".$ss3{$lookup}),"\t";
                    $foundflag = 1;
                }
                $newcoord = $coord-$dist;
                $lookup = $chrom.":".$newcoord."_-";
                if ((exists $ss3{$lookup})&&($foundflag==0)){
                    print FHout ($lookup.":".$ss3{$lookup}),"\t";
                    $foundflag = 1;
                }
                $dist++;
            }}
            else {print "error: $lookupkey[$j]\n";}
            if ($foundflag ==0){print FHout "null\t";}
        } 
        print FHout "\n";
            
    }
}

sub sameTranscript {
    $self = shift;
    $group = shift;
    
    $group =~ /(.+).txt/;
    $fhout = $self->{_outdir}."/".$1."_filter.txt";
    $fh = $self->{_outdir}."/".$group;
    open (FH, $fh) or die $!;
    open (FHout, ">$fhout") or die $!;
    
    
    while (<FH>){
        $line = $_;
        my @ss; 
        my @genes;
        if ($_ =~ /chr/){
            for (my $j = 18; $j<=25; $j++){
                push(@ss,((split /\t/,$line)[$j]));
            }
            for (my $j = 0; $j<8; $j++){
                if ($ss[$j] =~ /chr\w+:\d+.+:(.+)/){
                    $genelist = $1;
                    $genes[$j] = $genelist;
                }
                elsif ($ss[$j] =~ /null/) {
                    $genes[$j] = "";
                }
            }
            $headgenes = $genes[0].",".$genes[1].",".$genes[4].",".$genes[5];
            $tailgenes = $genes[2].",".$genes[3].",".$genes[6].",".$genes[7];
        
            @headgenes = (split /\,/,$headgenes);
            @tailgenes = (split /\,/,$tailgenes);
        
            $numheadgenes = @headgenes;
            $numtailgenes = @tailgenes;
        
            $flag = 0;
        
            for (my $j = 0; $j<$numheadgenes; $j++){
                for (my $k = 0; $k<$numtailgenes; $k++){
                    if ($headgenes[$j] eq $tailgenes[$k]){
                        $flag = 1;
                    }
                }
            }
        
            if ($flag==1){
                print FHout $line;
            }
        }
    }
    close FH;
    
    
}
        
        
        

sub findLariats {
    $self = shift;
    $group = shift;
    
    $group =~ /(.+).txt/;
    $fhout = $self->{_outdir}."/".$1."_lariats.txt";
    $fhoutoverlap = $self->{_outdir}."/".$1."_lariats_overlap.txt";
    $fhoutgap = $self->{_outdir}."/".$1."_lariats_gap.txt";
    $fh = $self->{_outdir}."/".$group;
    
    open (FH, $fh) or die $!;
    open (FHout,">$fhout") or die $!;
    open (FHoutoverlap,">$fhoutoverlap") or die $!;
    open (FHoutgap, ">$fhoutgap") or die $!;
    
    print "opening: $fh\n\n";
    
    while (<FH>){
        $line = $_;
        chomp $line;
        $seq = (split /\t/,$line)[1];
        $headseq = (split /\t/,$line)[6];
        $tailseq = (split /\t/,$line)[14];
        $readstrand = (split /\t/,$line)[3];
        #exact alignment
        if ((length($seq)) == ((length($headseq))+(length($tailseq)))){
            my @ss; my @ss_coord;
            $headcoord1 = (split /\t/,$_)[5];
            $headcoord2 = length($headseq)+$headcoord1;
            $tailcoord1 = (split /\t/,$line)[13];
            $tailcoord2 = length($tailseq)+$tailcoord1;
            for (my $j = 18; $j<=25; $j++){
                push(@ss,((split /\t/,$line)[$j]));
            }
            for (my $j = 0; $j<8; $j++){
                if ($ss[$j] =~ /chr\w+:(\d+)_([+-])/){
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
            $lariatflag = 0;
            #sense, positive strand
            if (($ss_coord[2][0]==$tailcoord1)&&($ss_coord[2][1] eq $readstrand)&&($readstrand eq "+")&&($ss_coord[2][1] ne undef)){
                $lariatflag = 1;
            }
            #sense, negative strand
            if (($ss_coord[3][0]==$tailcoord2)&&($ss_coord[3][1] eq $readstrand)&&($readstrand eq "-")&&($ss_coord[3][1] ne undef)){
                $lariatflag = 1;
            }
            #antisense, positive strand
            if (($ss_coord[0][0]==$headcoord1)&&($ss_coord[0][1] ne $readstrand)&&($readstrand eq "-")&&($ss_coord[0][1] ne undef)){
                $lariatflag = 1;
            }
            #antisense, negative strand
            if (($ss_coord[1][0]==$headcoord2)&&($ss_coord[1][1] ne $readstrand)&&($readstrand eq "+")&&($ss_coord[1][1] ne undef)){
                $lariatflag = 1;
            }
        
            if ($lariatflag==1){
                print FHout $line,"\n";
            }
        }
        #overlapped alignment
        elsif ((length($seq)) < ((length($headseq))+(length($tailseq)))){
            my @ss; my @ss_coord;
            $headcoord1 = (split /\t/,$_)[5];
            $headcoord2 = length($headseq)+$headcoord1;
            $tailcoord1 = (split /\t/,$line)[13];
            $tailcoord2 = length($tailseq)+$tailcoord1;
            $overlap = (length($headseq)+(length($tailseq)))-length($seq);
            for (my $j = 18; $j<=25; $j++){
                push(@ss,((split /\t/,$line)[$j]));
            }
            for (my $j = 0; $j<8; $j++){
                if ($ss[$j] =~ /chr\w+:(\d+)_([+-])/){
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
            $lariatflag = 0;
            #sense, positive strand
            #if ((abs($ss_coord[2][0]-$tailcoord1)<=$overlap)&&($ss_coord[2][1] eq $readstrand)&&($readstrand eq "+")&&($ss_coord[2][1] ne undef)){
            if ((abs($ss_coord[2][0]-$tailcoord1)<=$overlap)&&($ss_coord[2][0]>=$tailcoord1)&&($ss_coord[2][1] eq $readstrand)&&($readstrand eq "+")&&($ss_coord[2][1] ne undef)){
                $lariatflag = 1;
            }
            #sense, negative strand
            #if ((abs($ss_coord[3][0]-$tailcoord2)<=$overlap)&&($ss_coord[3][1] eq $readstrand)&&($readstrand eq "-")&&($ss_coord[3][1] ne undef)){
            if (((-1)*($ss_coord[3][0]-$tailcoord2)<=$overlap)&&($tailcoord2>=$ss_coord[3][0])&&($ss_coord[3][1] eq $readstrand)&&($readstrand eq "-")&&($ss_coord[3][1] ne undef)){
                   $lariatflag = 1;
            }
            #antisense, positive strand
            #if ((abs($ss_coord[0][0]-$headcoord1)<=$overlap)&&($ss_coord[0][1] ne $readstrand)&&($readstrand eq "-")&&($ss_coord[0][1] ne undef)){
            if ((($ss_coord[0][0]-$headcoord1)<=$overlap)&&($ss_coord[0][0] >=$headcoord1)&&($ss_coord[0][1] ne $readstrand)&&($readstrand eq "-")&&($ss_coord[0][1] ne undef)){
                 $lariatflag = 1;
            }
            #antisense, negative strand
            #if ((abs($ss_coord[1][0]-$headcoord2)<=$overlap)&&($ss_coord[1][1] ne $readstrand)&&($readstrand eq "+")&&($ss_coord[1][1] ne undef)){
            if (((-1)*($ss_coord[1][0]-$headcoord2)<=$overlap)&&($headcoord2>=$ss_coord[1][0])&&($ss_coord[1][1] ne $readstrand)&&($readstrand eq "+")&&($ss_coord[1][1] ne undef)){
                $lariatflag = 1;
            }
        
            if ($lariatflag==1){
                print FHoutoverlap $line,"\n";
            }
        }
        #gapped alignment
        elsif ((length($seq)) > ((length($headseq))+(length($tailseq)))){
            my @ss; my @ss_coord;
            $headcoord1 = (split /\t/,$_)[5];
            $headcoord2 = length($headseq)+$headcoord1;
            $tailcoord1 = (split /\t/,$line)[13];
            $tailcoord2 = length($tailseq)+$tailcoord1;
            $gap = length($seq)-(length($headseq)+(length($tailseq)));
            for (my $j = 18; $j<=25; $j++){
                push(@ss,((split /\t/,$line)[$j]));
            }
            for (my $j = 0; $j<8; $j++){
                if ($ss[$j] =~ /chr\w+:(\d+)_([+-])/){
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
            $lariatflag = 0;
            #sense, positive strand
            if ((abs($ss_coord[2][0]-$tailcoord1)<=$gap)&&($ss_coord[2][1] eq $readstrand)&&($readstrand eq "+")&&($ss_coord[2][1] ne undef)){
                $lariatflag = 1;
            }
            #sense, negative strand
            if ((abs($ss_coord[3][0]-$tailcoord2)<=$gap)&&($ss_coord[3][1] eq $readstrand)&&($readstrand eq "-")&&($ss_coord[3][1] ne undef)){
                $lariatflag = 1;
            }
            #antisense, positive strand
            if ((abs($ss_coord[0][0]-$headcoord1)<=$gap)&&($ss_coord[0][1] ne $readstrand)&&($readstrand eq "-")&&($ss_coord[0][1] ne undef)){
                $lariatflag = 1;
            }
            #antisense, negative strand
            if ((abs($ss_coord[1][0]-$headcoord2)<=$gap)&&($ss_coord[1][1] ne $readstrand)&&($readstrand eq "+")&&($ss_coord[1][1] ne undef)){
                $lariatflag = 1;
            }
        
            if ($lariatflag==1){
                print FHoutgap $line,"\n";
            }
        }
        
    }
    close FH; close FHoutoverlap; close FHout; close FHoutgap;
    
}

sub resolveGaps {
    
    $self = shift;
    $group = shift;
    
    $group =~ /(.+).txt/;
    $base = $1;
    $fhhead = $self->{_outdir}."/".$base."_heads.bed";
    $fhtail = $self->{_outdir}."/".$base."_tails.bed";
    $fh = $self->{_outdir}."/".$group;
    
    open (FH, $fh) or die $!;
    open (FHhead,">$fhhead") or die $!;
    open (FHtail,">$fhtail") or die $!;
    
     while (<FH>){
        $line = $_;
        chomp $line;
        $id = (split /\t/,$line)[0];
        $chrom = (split /\t/,$line)[4];
        $seq = (split /\t/,$line)[1];
        $headseq = (split /\t/,$line)[6];
        $tailseq = (split /\t/,$line)[14];
        $readstrand = (split /\t/,$line)[3];
    
        my @ss; my @ss_coord;
        $headcoord1 = (split /\t/,$_)[5];
        $headcoord2 = length($headseq)+$headcoord1;
        $tailcoord1 = (split /\t/,$line)[13];
        $tailcoord2 = length($tailseq)+$tailcoord1;
        $gap = length($seq)-(length($headseq)+(length($tailseq)));
    
        for (my $j = 18; $j<=25; $j++){
            push(@ss,((split /\t/,$line)[$j]));
        }
        for (my $j = 0; $j<8; $j++){
            if ($ss[$j] =~ /chr\w+:(\d+)_([+-])/){
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
    
        #EDITED DISTANCE5ss >=0 from >0 IN ALL CASES
        #sense, positive strand
        if ((abs($ss_coord[2][0]-$tailcoord1)<=$gap)&&($ss_coord[2][1] eq $readstrand)&&($readstrand eq "+")&&($ss_coord[2][1] ne undef)){
            $dist5ss = $tailcoord1-$ss_coord[2][0];
            if ($dist5ss >= 0){
                $tailcoord1 = $ss_coord[2][0];
                $headcoord2 = $headcoord2+($gap-$dist5ss);
                print FHhead $chrom,"\t",$headcoord1,"\t",$headcoord2,"\t",$id,"\t0\t",$readstrand,"\n";
                print FHtail $chrom,"\t",$tailcoord1,"\t",$tailcoord2,"\t",$id,"\t0\t",$readstrand,"\n";
            }

        }
        #sense, negative strand
        elsif ((abs($ss_coord[3][0]-$tailcoord2)<=$gap)&&($ss_coord[3][1] eq $readstrand)&&($readstrand eq "-")&&($ss_coord[3][1] ne undef)){
            $dist5ss = $ss_coord[3][0] - $tailcoord2;
            if ($dist5ss >= 0){
                $tailcoord2 = $ss_coord[3][0];
                $headcoord1 = $headcoord1 - ($gap-$dist5ss);
                print FHhead $chrom,"\t",$headcoord1,"\t",$headcoord2,"\t",$id,"\t0\t",$readstrand,"\n";
                print FHtail $chrom,"\t",$tailcoord1,"\t",$tailcoord2,"\t",$id,"\t0\t",$readstrand,"\n";
            }
        }
        #antisense, positive strand
        elsif ((abs($ss_coord[0][0]-$headcoord1)<=$gap)&&($ss_coord[0][1] ne $readstrand)&&($readstrand eq "-")&&($ss_coord[0][1] ne undef)){
            $dist5ss = $headcoord1-$ss_coord[0][0];
            if ($dist5ss >= 0){
                $headcoord1 = $ss_coord[0][0];
                $tailcoord2 = $tailcoord2+($gap-$dist5ss);
                print FHhead $chrom,"\t",$headcoord1,"\t",$headcoord2,"\t",$id,"\t0\t",$readstrand,"\n";
                print FHtail $chrom,"\t",$tailcoord1,"\t",$tailcoord2,"\t",$id,"\t0\t",$readstrand,"\n";
            }
        }
        #antisense, negative strand
        elsif ((abs($ss_coord[1][0]-$headcoord2)<=$gap)&&($ss_coord[1][1] ne $readstrand)&&($readstrand eq "+")&&($ss_coord[1][1] ne undef)){
            $dist5ss = $ss_coord[1][0] - $headcoord2;
            if ($dist5ss >= 0){
                $headcoord2 = $ss_coord[1][0];
                $tailcoord1 = $tailcoord1 - ($gap-$dist5ss);
                print FHhead $chrom,"\t",$headcoord1,"\t",$headcoord2,"\t",$id,"\t0\t",$readstrand,"\n";
                print FHtail $chrom,"\t",$tailcoord1,"\t",$tailcoord2,"\t",$id,"\t0\t",$readstrand,"\n";
            }
        }
    }
    close FH; close FHhead; close FHtail;
    
    $command1 = "bedtools getfasta -fi ".$self->{_index}.".fa -bed ".$fhhead." -fo ".$self->{_outdir}."/".$base."_heads.fasta -name -tab -s";
    $command2 = "bedtools getfasta -fi ".$self->{_index}.".fa -bed ".$fhtail." -fo ".$self->{_outdir}."/".$base."_tails.fasta -name -tab -s";
    
    print "\n\n$command1\n\n$command2\n\n";
    
    system($command1);
    system($command2);
    
    my %gappedreads;
    $fhhead =~ s/bed/fasta/;
    $fhtail =~ s/bed/fasta/;
    
    print "opening:\t$fhhead\n";
    print "opening:\t$fhtail\n";
    
    open (FHhead, $fhhead) or die $!;
    open (FHtail, $fhtail) or die $!;
    while (<FHhead>){
        ($id,$seq) = (split /\t/,$_);
        chomp $seq;
        $gappedreads{$id}[0] = $seq;
    }
    while (<FHtail>){
        ($id,$seq) = (split /\t/,$_);
        chomp $seq;
        $gappedreads{$id}[1] = $seq;
    }
    close FHhead; close FHtail;
    
    open (FH, $fh) or die $!;
    $fhout = $self->{_outdir}."/".$base."_truelariats.txt";
    open (FHout, ">$fhout") or die $!;
    
    while (<FH>){
        $line = $_;
        chomp $line;
        
        $id = (split /\t/,$line)[0];
        $seq = (split /\t/,$line)[1];
        $headseq = (split /\t/,$line)[6];
        $tailseq = (split /\t/,$line)[14];
        $headcoord1 = (split /\t/,$_)[5];
        $headcoord2 = length($headseq)+$headcoord1;
        $tailcoord1 = (split /\t/,$line)[13];
        $tailcoord2 = length($tailseq)+$tailcoord1;
        
        
        ####print output file
        if (exists $gappedreads{$id}[0]){
            $genomic = $gappedreads{$id}[0].$gappedreads{$id}[1];
            $stringdist = &levenshtein($seq,$genomic);
            if ($stringdist<=3){
                print FHout $line,"\t",$genomic,"\t",$stringdist,"\n";
            }
        }
    }
    
}





sub levenshtein
{
    my ($s1, $s2) = @_;
    my ($len1, $len2) = (length $s1, length $s2);

    return $len2 if ($len1 == 0);
    return $len1 if ($len2 == 0);

    my %mat;

    for (my $i = 0; $i <= $len1; ++$i)
    {
        for (my $j = 0; $j <= $len2; ++$j)
        {
            $mat{$i}{$j} = 0;
            $mat{0}{$j} = $j;
        }

        $mat{$i}{0} = $i;
    }

    my @ar1 = split(//, $s1);
    my @ar2 = split(//, $s2);

    for (my $i = 1; $i <= $len1; ++$i)
    {
        for (my $j = 1; $j <= $len2; ++$j)
        {
            my $cost = ($ar1[$i-1] eq $ar2[$j-1]) ? 0 : 1;

            $mat{$i}{$j} = min([$mat{$i-1}{$j} + 1,
                                $mat{$i}{$j-1} + 1,
                                $mat{$i-1}{$j-1} + $cost]);
        }
    }

    return $mat{$len1}{$len2};
}


sub min
{
    my @list = @{$_[0]};
    my $min = $list[0];

    foreach my $i (@list)
    {
        $min = $i if ($i < $min);
    }

    return $min;
}
    

1;
