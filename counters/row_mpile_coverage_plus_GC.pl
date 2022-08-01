#!/usr/bin/perl


while(<>)
  {
    chomp;
    
    @tmp = split(/\s+/,$_);
    @row = ();
    next unless ($tmp[2] eq 'G' || $tmp[2] eq 'g' || $tmp[2] eq 'C' || $tmp[2] eq 'c');

print $tmp[0],"\t",$tmp[1],"\t",$tmp[2],"\t",$tmp[3],"\t";

	my $c=0;
       if($tmp[2] eq 'G' || $tmp[2] eq 'g')
	 {
	   	$tmp[4]=~ s/([cC])/$1; $c++/eg;
	   
	 }
       else
	 {
	   	$tmp[4]=~ s/([Gg])/$1; $c++/eg;

	 }

	print $c;

    print "\n";
  }

