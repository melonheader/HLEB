#!/usr/bin/perl


while(<>)
  {
    chomp;
    
    @tmp = split(/\s+/,$_);
    @row = ();
    next unless ($tmp[2] eq 'T' || $tmp[2] eq 't' || $tmp[2] eq 'A' || $tmp[2] eq 'a');

print $tmp[0],"\t",$tmp[1],"\t",$tmp[2],"\t",$tmp[3],"\t";

	my $c=0;
       if($tmp[2] eq 'T' || $tmp[2] eq 't')
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

