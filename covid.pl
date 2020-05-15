#!/usr/bin/perl

use DateTime;
use DateTime::Duration;
use Time::Piece;



my $SQRT2PI = 2.506628274631;

sub pdf {
  my ( $x, $m, $s ) = ( 0, 0, 1 );
  $x = shift if @_;
  $m = shift if @_;
  $s = shift if @_;

  if( $s <= 0 ) {
    croak( "Can't evaluate Math::Gauss:pdf for \$s=$s not strictly positive" );
  }

  my $z = ($x-$m)/$s;

  return exp(-0.5*$z*$z)/($SQRT2PI*$s);
}

sub cdf {
  my ( $x, $m, $s ) = ( 0, 0, 1 );
  $x = shift if @_;
  $m = shift if @_;
  $s = shift if @_;

  # Abramowitz & Stegun, 26.2.17
  # absolute error less than 7.5e-8 for all x

  if( $s <= 0 ) {
    croak( "Can't evaluate Math::Gauss:cdf for \$s=$s not strictly positive" );
  }

  my $z = ($x-$m)/$s;

  my $t = 1.0/(1.0 + 0.2316419*abs($z));
  my $y = $t*(0.319381530
          + $t*(-0.356563782
            + $t*(1.781477937
              + $t*(-1.821255978
                + $t*1.330274429 ))));
  if( $z > 0 ) {
    return 1.0 - pdf( $z )*$y;
  } else {
    return pdf( $z )*$y;
  }
}

sub gaussian {
    # adjust for end effects
    my ($x, $mean, $var) = @_;
    if ($x>($mean+(2*$var))) { return 0; }
    my $e=cdf(1, $mean, $var) + (1- cdf($mean+(2*$var),$mean,$var));
    return pdf($x,$mean,$var) * ($e+1);
}


my @rcases;

sub run {
    my ($verbose, $R0_pop, $R0_ld, $importperday) = @_;

    my $lockdownday=75;
    my $deconfineday=136;
    my $R0_dld=1.2;

    my $day = Time::Piece->strptime("20191227", "%Y%m%d"); # day 0

    my @casereports;
    my @infected;
    $infected[0]=1;     # start with 1 cases

    my $infectday_m=4;  # syptoms appear ~ day 5, 1 day before is most infectious
    my $infectday_v=5;

    my $unreported=0.9; # number of unreported cases 90%
    # Spain has ~50M 5% have antibodies, 2.5M, reported number 250k


    my $casedelay_m=25;  #day cases are notified (e.g. delay to get tested, delay for result, delay to record)
    #Can be measured (lockdown to peak)
    my $casedelay_v=1;

    my $cases=0;
    my $score=0;



    $oldcases=0;
    $verbose and print '"day","date","New infections","recorded cases","delta"'."\n";
    foreach $d (1..200) {

        if ($rcases[$d]) {
            $score+=(abs($rcases[$d] - $cases));
        }
        
        if ($d == $lockdownday) {
            $R0_pop=$R0_ld;
        }
        if ($d == $deconfineday) {
            $R0_pop=$R0_dld;
        }
        
        $verbose and print $d.", ".$day->dmy('/').", $infected[0], $cases, ".($cases-$oldcases)."\n";
        $oldcases=$cases;
        $day+=60*60*24;
        # day +1
        foreach $ii (1..100) {
            $i=101-$ii;
            $infected[$i] = $infected[$i-1];
            $casereports[$i] = $casereports[$i-1];
        }

        $infected[0]=0;
        $casereports[0]=0;

        $infected[0]+=$importperday;

        foreach $i (1..100) {

            $nc=$R0_pop * $infected[$i] * gaussian($i, $infectday_m, $infectday_v);

            $infected[0] += $nc;
            $casereports[0] += $nc*(1-$unreported);

            if ($i<$casedelay_m) {
                $cr=$casereports[$i]*gaussian($i,$casedelay_m, $casedelay_v);
            } else {
                $cr=$casereports[$i]*gaussian($casedelay_m,$casedelay_m, $casedelay_v);
            }
            $cr*=1-(gaussian(($d+(1.5))%7,3,0.4) * 0.6); # 40% weekend effect.
            $cases+=$cr;
            $casereports[$i]-=$cr;
        }
    }

    $verbose and print "Score: $score\n";
    return $score;
}


open (my $data, '<', 'france.csv') or die "Cant find csv\n";
while (my $line = <$data>) {
  chomp $line;
 
  my @fields = split "," , $line;
  $rcases[$fields[0]]=$fields[1];
}



sub findlocalmin {
    my ($f_ref) = @_;

    my $R0=@{$f_ref}[0];
    my $R0ld=@{$f_ref}[1];
    my $imp=@{$f_ref}[2];
    $best=run(0, @{$f_ref}[0], @{$f_ref}[1], @{$f_ref}[2]);

    foreach $i (1..4) {
        foreach $f (0..2) {
            my $d=10;
            while (abs($d)<100) {

                $old=@{$f_ref}[$f];
                @{$f_ref}[$f] += @{$f_ref}[$f]*(1/$d);
                $s=run(0, @{$f_ref}[0], @{$f_ref}[1], @{$f_ref}[2]);
                if ($s < $best) {
                    $best=$s;
                    $d=abs($d*2);
                } else {
                    @{$f_ref}[$f]=$old;
                    if ($d>0) {
                        $d*=-1;
                    } else {
                        $d=abs($d*2);
                    }
                }
            }
        }
    }
    return $best;
}


if (scalar @ARGV == 3) {
    run(1, @ARGV[0], $ARGV[1], $ARGV[2]);
    exit(-1);
}
    
my $R0=2;#1.88;
my $R0ld=0.5;#0.7;
my $imp=5;#10;
my $bestscore=run(0,1,0.1,0);

my @facs=(2,0.5,5);
my @facs_min=(1.0, 0.1,  0);
my @facs_max=(5.0, 2.0, 50);
my @bestfacs;
foreach $t (1..80) {
    $temp=(1 - $t/100 );
    foreach $i (0..2) {
        my $max=$facs[$i]+(($facs_max[$i]-$facs[$i]) * $temp);
        my $min=$facs[$i]-(($facs[$i]-$facs_min[$i]) * $temp);
        #        $facs[$i] = rand($facs_max[$i]-$facs_min[$i])+$facs_min[$i];
        $facs[$i] = rand($max-$min)+$min;
    }
    my $olds=$s;
    my $s=findlocalmin(\@facs);
    print "$t) tried R0 $facs[0], R0 lockdown $facs[1], Import/Day $facs[2] - score : $s\n";
    if ($s<$bestscore) {
        print "Found one\n";
        $bestscore=$s;
        @bestfacs=@facs;
    }
}

print "R0 $bestfacs[0], R0 lockdown $bestfacs[1], Import/Day $bestfacs[2]\n";



# R0 1.83510069563299, R0 lockdown 0.733478820021408, Import/Day 13.511478613536 - score : 142742.515308867
