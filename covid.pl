#!/usr/bin/perl
use strict;
use DateTime;
use DateTime::Duration;
use Time::Piece;
use threads;
use threads::shared;



my @facs_names=("Base R0 for virus", "Mobility multiply", "Daily Imported cases", "Infectious day mean", "Infectious day variance", "Case reporting delay mean", "Case reporting delay variance", "Social Distancing effect", "Start day");
my @facs_min=(1.5, 0.5, 1,  3.5,  4, 15, 3, 0.5, 0  );
my @facs=    (1.8, 1,   10,  5,   5, 19, 4, 0.8, 1  );
my @facs_max=(2.5, 2.5, 20,  7,   7, 25, 6, 1.2, 30  );


if (scalar @ARGV ==0) {
    die "$0 country [factors]\n";
}
my $country=shift @ARGV;


my %rcases;
my %mobility;
my $lastmobility=0;
my $lastmobilityday;
my $weekendday=0;
my $weekendeffect=0;


# Read in the data we have from mobility and cases
my @casesperday;
open (my $data, '<', "$country.csv") or die "Cant find csv\n";
my $oldcases=0;
while (my $line = <$data>) {
  chomp $line;

  my @fields = split "," , $line;
  $rcases{$fields[0]}=int($fields[1]);

  my $day = Time::Piece->strptime($fields[0], "%Y-%m-%d"); # day 0
  $casesperday[$day->day_of_week]+=int($fields[1])-$oldcases;
  $oldcases=$fields[1];
}
close ($data);
my $total=0;
foreach my $d (0..6) {
    ($casesperday[$d]<$casesperday[$weekendday]) and $weekendday=$d;
    $total+=$casesperday[$d];
}
$weekendeffect=$casesperday[$weekendday]/($total/7);

open (my $data, '<', "mobility-$country.csv") or die "Cant find csv\n";
while (my $line = <$data>) {
  chomp $line;

  my @fields = split "," , $line;
  if ($fields[4]) {
      my $m=($fields[5]+$fields[6]+$fields[7]+$fields[8]+$fields[9]-$fields[10])/6.0;
      $mobility{$fields[4]}=$m;
      $lastmobility=(($lastmobility*2)+$m)/3;
      $lastmobilityday=Time::Piece->strptime($fields[4],"%Y-%m-%d");
  }
}
close ($data);



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
    # adjust for end effects, assume numbers go 1-upwards
    my ($x, $mean, $var) = @_;
    if ($x>($mean+(2*$var))) { return 0; }
    my $e=cdf(1, $mean, $var) + (1- cdf($mean+(2*$var),$mean,$var));
    return pdf($x,$mean,$var) * ($e+1);
}


sub run {
    my ($verbose, $f_ref) = @_;
    my @facs=@{$f_ref};

    my $R0_pop=$facs[0];
    my $mob_fac=$facs[1];
    my $importperday=$facs[2];

    my $SD_introday=75;

    #    my $day = Time::Piece->strptime("20191227", "%Y%m%d"); # day 0
    my $day = Time::Piece->strptime("20200101", "%Y%m%d"); # day 0
    $day+=60*60*24*$facs[8];

    my @casereports;
    my @infected;

    my $infectday_m=$facs[3];  # syptoms appear ~ day 5, 1 day before is most infectious
    my $infectday_v=$facs[4];

    my $unreported=0.9; # number of unreported cases 90%
    # Spain has ~50M 5% have antibodies, 2.5M, reported number 250k


    my $casedelay_m=$facs[5];

    my $casedelay_v=$facs[6];
    my $SD_effect=$facs[7];

    my $cases=0;
    my $score=0;


    my $oldcases=0;
    $verbose and print '"day","date","Estimated R0", "New infections","simulated cases","Actual cases", "Error", "Daily simulated delta"'."\n";
    my $mobilityday=0;;
    foreach my $d (1..200) {

        my $err=0;
        if ($rcases{$day->ymd('-')}) {
            $err=($rcases{$day->ymd('-')} - $cases);
            $score+=($err*$err)/(200-$d);
        }

            
        my $R0=$R0_pop;
        my $SDe=1;
        ($d>=$SD_introday) and $SDe=$SD_effect;
        if ($day > $lastmobilityday) {
            $R0=$R0_pop * (((100+$lastmobility) * $mob_fac * $SDe)/100);#$R0;#_current;
        } else {
            if ($mobility{$day->ymd('-')}) {
                $R0=$R0_pop * (((100+$mobility{$day->ymd('-')}) * $mob_fac * $SDe)/100);
            } else {
                $R0=$R0_pop * (((100) * $mob_fac * $SDe)/100);
            }
        }

        foreach my $ii (1..100) {
            my $i=101-$ii;
            $infected[$i] = $infected[$i-1];
            $casereports[$i] = $casereports[$i-1];
        }

        $infected[0]=0;
        $casereports[0]=0;

        $infected[0]+=$importperday;

        my $dow=$day->day_of_week(); #0 is Sunday

        foreach my $i (1..100) {

            my $nc=$R0 * $infected[$i] * gaussian($i, $infectday_m, $infectday_v);

            $infected[0] += $nc;
            $casereports[0] += $nc*(1-$unreported);

            my $cr;
            if ($i<$casedelay_m) {
                $cr=$casereports[$i]*gaussian($i,$casedelay_m, $casedelay_v);
            } else {
                $cr=$casereports[$i]*gaussian($casedelay_m,$casedelay_m, $casedelay_v);
            }
            if ($dow==$weekendday) { $cr*=$weekendeffect; }
            $cases+=$cr;
            $casereports[$i]-=$cr;
        }
        
        $verbose and print $d.", ".$day->dmy('/').", $R0,  $infected[0], $cases, $rcases{$day->ymd('-')}, $err, ".($cases-$oldcases)."\n";
        $oldcases=$cases;
        $day+=60*60*24;
        # day +1
    }

    $verbose and print "Score: $score\n";
    return $score;
}


sub findlocalmin {
    my ($f_ref, $best) = @_;
    my $dinit=10;
    limitloop: foreach my $i (1..10) {
        my $changes=0;
        foreach my $f (0..(scalar @{$f_ref} -1)) {
            if ($facs_min[$f]!=$facs_max[$f]) {
                my $d=$dinit;
                my $found=0;
                while (abs($d)<=10000 && $found<10) {
                    my $old=@{$f_ref}[$f];
                    @{$f_ref}[$f] += @{$f_ref}[$f]*(1/$d);
                    my $s=run(0, $f_ref);
                    if ($s < $best) {
                        $changes+=($best-$s);
                        $best=$s;
                        $found+=1;
                        if ($found>=5) {
                            $d/=2;
                            $found=0;
                            if (abs($d)<2) {
                                return $best;
                            }
                        }
                    } else {
                        @{$f_ref}[$f]=$old;
                        if (($d>0) && $found==0) {
                            $d*=-1;
                        } else {
                            $found=0;
                            $d=abs($d*10);
                        }
                    }
                }
            }
        }
        ($changes<10) and  last limitloop;
#        print "$best\n";
        $dinit*=2;
    }
    return $best;
}



my $lock :shared=0;
my @bestfacs :shared;
my $bestscore :shared;
my $testid :shared;

sub runandtests
{
    my $temp=1;
    
    while ($testid < 95000) {# && $temp>0.1) {
        my @test;
        {
            lock($lock);
            $testid+=1;
            printf "\r%6d ",$testid;
            STDOUT->flush();
            my $temp=(1 - $testid/100000 );
            @test=@facs;
            my @current=@bestfacs;
            foreach my $i (0..(scalar @facs-1)) {
                if ($facs_min[$i]==$facs_max[$i]) {
                    $test[$i]=$facs_max[$i];
                } else {
                    my $max=$current[$i]+(($facs_max[$i]-$current[$i]) * $temp);
                    my $min=$current[$i]-(($current[$i]-$facs_min[$i]) * $temp);
                    $test[$i] = rand($max-$min)+$min;
                }
            }
        }
        my $locals=run(0, \@test);
        {
            if ($locals < $bestscore*50) {
                my $s=findlocalmin(\@test, $locals);
                foreach my $i (0..(scalar @facs)-1) {
                    printf "%-*f\t",length($facs_names[$i]), $test[$i];
                }
                print " Score= $s (".($locals/$s).") ";
                {
                    lock($lock);
                    if ($s<$bestscore) {

                        print "(*)";
                        $bestscore=$s;
                        @bestfacs=@test;

#                        $temp*=0.95;
                    }
                }
                print "\n";
            }
        }
    }
}

if (scalar @ARGV > 2 && $ARGV[0] ne "-c") {
    my @test=@ARGV;
    my $s=run(1, \@test);
    exit(-1);
}

if (scalar @ARGV > 2 && $ARGV[0] eq "-c") {
    shift @ARGV;
    @facs=@ARGV;
    print "Continuing...\n";
}

@bestfacs  =@facs;
$bestscore =run(0, \@facs);

print "       ",join("\t",@facs_names),"\n";
my @threads;
foreach my $thrid (1..8) {
    push @threads, threads->create(\&runandtests);
}
foreach my $thr (@threads) {
    $thr->join();
}

print join(' ',@facs_names),"\n";
print join(' ',@bestfacs),"\n";
print " Score= $bestscore\n";

