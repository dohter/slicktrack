#!/bin/perl -w
# (1) quit unless we have the correct number of command-line args
$num_args = $#ARGV + 1;
if ($num_args != 2) {
    print "\nUsage: scan.pl latname np_s\n";
    exit;
}
# (2) we got two command line args, so assume they are the
# first name and last name
$name=$ARGV[0];
$np_s=$ARGV[1];

#$E00=9.25;
$E00=17.6;
$ECAV=10.0;
$NSTEP=250/$np_s;
$DE0=0.002;
print "Scanning $name on $np_s processes\n";

$inFile="./template.in";
print "$inFile\n";
for( $n = 1; $n <= $np_s; $n = $n+1){
  open(FILE,"$inFile") || die "cannot open file $inFile!" ;
  open(OUTFILE,">./$name.in.$n") || die "cannot open file $name.in.$n!" ;
  # Setup the outfile
  while( <FILE> )
    {
      $_ =~ s/NNN/$NSTEP/;
      $_ =~ s/PPPP/$n/;
      $_ =~ s/EEEE/$E00/;
      $_ =~ s/DDDDD/$DE0/;
      $_ =~ s/LLLL/$name/;
      $_ =~ s/CCCCC/$ECAV/;
      print OUTFILE $_;
    }
  close( OUTFILE );
  close( FILE );
  open(OUTFILE,">./slicktrack$n.sh") || die "cannot open file ./slicktrack$n.sh!" ;
  print OUTFILE "\#!/bin/bash\n";
  print OUTFILE "/home/obeznosov/git/slicktrack/slicktrack < /home/obeznosov/git/slicktrack/rundir/$name.in.$n\n";
  system("qsub -q all.q slicktrack$n.sh");
  print "Started slicktrack process number $n with starting energy $E00, $NSTEP steps\n";
  $E00=$n*$NSTEP*$DE0
  }
exit

#
#}
#
#exit
