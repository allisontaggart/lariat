#filter.pm

package filter;

use warnings;

sub new
{
    my $class = shift;
    my $self = {
        _rawfile => shift,
        _outdir => shift,
    };
    
    $self->{_rawfile} =~ /(.+).fa?s?t?q/;
    $self->{_file}=$self->{_outdir}."/new_alignments.txt";
    
    bless $self, $class;
    return $self;
    
}

sub outOfOrder
{
    my $self = shift;
    print "\n\n****************\nfiltering out of order reads\n\nopening file: ";
    print $self->{_file},"\n\n";
    
    open (FH, $self->{_file});

    $halffile = $self->{_outdir}."/halfmap.txt";
    $outoforderfile = $self->{_outdir}."/outoforder.txt";
    
    open (FHhalf,">$halffile") or die $!;
    open (FHoutoforder,">$outoforderfile") or die $!;
    
    while (<FH>){
    
        $line = $_;
        $readid = (split /\t/,$_)[0];
        $readseq = (split /\t/,$_)[1];
        $headstrand = (split /\t/,$_)[3];
        $headchrom = (split /\t/,$_)[4];
        $headcoord = (split /\t/,$_)[5];
        $tailstrand = (split /\t/,$_)[11];
        $tailchrom = (split /\t/,$_)[12];
        $tailcoord = (split /\t/,$_)[13];
        
        if (($headchrom ne 0)&&($tailchrom ne 0)){
            if (($headchrom eq $tailchrom)&&($headstrand eq $tailstrand)){
                if (($headstrand eq "+")&&($headcoord > $tailcoord)){
                    print FHoutoforder $line;
                }
                elsif (($headstrand eq "-")&&($headcoord < $tailcoord)){
                    print FHoutoforder $line;
                }
            }
        } 
        else {
            print FHhalf $line;
        }              
    }
    close FH;
    close FHhalf;
    closeFHoutoforder;
}

sub setFile
{
    my $self = shift;
    $self->{_file} = shift;

}

sub getFile
{
    my $self = shift;
    return $self->{_file};
}
    
1;
