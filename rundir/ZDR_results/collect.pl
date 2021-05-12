#!/usr/bin/perl -w
# (1) quit unless we have the correct number of command-line args
$num_args = $#ARGV + 1;
if ($num_args != 1) {
    print "\nUsage: collect.pl np_s\n";
    exit;
}
# (2) we got two command line args, so assume they are the
# first name and last name
$np_s=$ARGV[0];

print "Collecting output from $np_s processes\n";
$name="mctdep";
$se="_";
open(OUTFILE,">./$name.out") || die "cannot open file $name.in.$n!" ;
for( $n = 1; $n <= $np_s; $n = $n+1){
  open(FILE,"$name$se$n.out") || die "cannot open file $name$se$n.out!" ;
  # Setup the outfile
  while( <FILE> )
  {
    print OUTFILE $_;
  }
  close( FILE );
  }
close( OUTFILE );

$name="pol";
open(OUTFILE,">./$name.out") || die "cannot open file $name.in.$n!" ;
for( $n = 1; $n <= $np_s; $n = $n+1){
  open(FILE,"$name$se$n.out") || die "cannot open file $name$se$n.out!" ;
  # Setup the outfile
  while( <FILE> )
  {
    print OUTFILE $_;
  }
  close( FILE );
  }
close( OUTFILE );
exit

#
#}
#
#exit
