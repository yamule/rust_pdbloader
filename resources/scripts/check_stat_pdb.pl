use strict;
use warnings;



sub getAllFiles{
    my $targetdir = $_[0];
    $targetdir =~ s/[\\\/]+$//g;
    $targetdir .= "/";
    opendir(TARGET,$targetdir);
    my @allfiles = grep(!/^[\.]+$/,readdir(TARGET));
    closedir(TARGET);
    
    my @ret;
    foreach my $ff(@allfiles){
        if(-f $targetdir.$ff){
            push(@ret,$targetdir.$ff);
        }
        if(-d  $targetdir.$ff){
            push(@ret,getAllFiles( $targetdir.$ff));
        }
    }
    return @ret;
       
}


my @allfiles = getAllFiles("D:/dummy/vbox_share/bioo/database/pdb/");

my $cpunum = 8;
for(my $cc = 0;$cc < $cpunum;$cc++){
	my $pid = fork();
	
	if($pid){
	}else{
		my $fh;
		open($fh,"> res$cc.dat");
		for(my $ii = 0;$ii <= $#allfiles;$ii++){
			if($ii%$cpunum == $cc){
				my $fname = $allfiles[$ii];
				if($fname =~ /\.gz$/){
					my $ress = `python pdb_prop_check.py $fname`;
					$ress =~ s/[\r\n]+$//g;
					print $fh $ress."\n";
				}
			}
			#last;
		}
		close($fh);
		last;
	}
}

