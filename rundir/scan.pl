#!/usr/bin/perl -w
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
$E00=17.54;
$NSTEP=320/$np_s;
$DE0=0.002;
print "Scanning $name on $np_s processes\n";

$inFile="./$name.in.T";
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
      print OUTFILE $_;
    }
  close( OUTFILE );
  close( FILE );
  system("nohup ../slicktrack < ./$name.in.$n > nohup.$name.$n&");
  print "Started slicktrack process number $n with starting energy $E00, $NSTEP steps\n";
  $E00+=($NSTEP+1)*$DE0
  }
exit

#
#}
#
#exit
