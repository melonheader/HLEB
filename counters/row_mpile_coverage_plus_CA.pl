#!/usr/bin/perl


while(<>)
  {
    chomp;
    
    @tmp = split(/\s+/,$_);
    @row = ();
    next unless ($tmp[2] eq 'C' || $tmp[2] eq 'c' || $tmp[2] eq 'G' || $tmp[2] eq 'g');

print $tmp[0],"\t",$tmp[1],"\t",$tmp[2],"\t",$tmp[3],"\t";

	my $c=0;
       if($tmp[2] eq 'C' || $tmp[2] eq 'c')
	 {
	   	$tmp[4]=~ s/([aA])/$1; $c++/eg;
	   
	 }
       else
	 {
	   	$tmp[4]=~ s/([Tt])/$1; $c++/eg;

	 }

	print $c;

    print "\n";
  }

