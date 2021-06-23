use strict;
use warnings;


my $divdir = "D:/dummy/vbox_share/bioo/database/pdb/";
my $fasfile = "results/clustered.dat.fas";

open(IN,$fasfile);
my %candidates;
while(my $ss = <IN>){
	if($ss =~ />[\s]*([^\s]+)/){
		my $name = $1;
		if($name =~ /([0-9][0-9a-zA-Z][0-9a-zA-Z][0-9a-zA-Z])[:_]([0-9a-zA-Z]+)/){
			my $code = lc $1;
			$code =~ /^.(..)./;
			my $dir = $1;
			my $path = $divdir."/".$dir."/pdb".$code.".ent.gz";
			if(-f $path){
				$candidates{$path} = 1000;
			}else{
				print $path." was not found.\n";
			}
			
			
		}else{
			print $name." was not parsed.\n";
		}
	}else{
		if($ss =~ />/){
			print $ss." was not parsed.\n";
		}
	}
}

close(IN);


foreach my $kk(keys %candidates){
	print $kk."\n";
}