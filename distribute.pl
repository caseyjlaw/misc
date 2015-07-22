#!/usr/bin/perl
#
# Script to distribute processing jobs on cluster.
#
# v01: Adds capability to submit only one scan per node, with dynamic
#      subsequent submission as jobs end on nodes.
# v02: Distributes DM ranges across nodes in addition to scans.
# v03: Produces pickle and ms files before processing job submission.
#      Also runs batch_plot.py.
# v04: Now passes actual obs scan to batch_*.py, rather than
#      0-index scan.
# v05: Now ignores ready/finished files from other sessions (allowing parallel jobs) and takes nodelist as argument
#
#use strict;
#use warnings;

use Term::ReadKey;
use POSIX;
#use Net::SSH::Perl;


  # Read inputs and set up DM ranges
$asdm_file  = $ARGV[0];
$ms_base = $ARGV[1];
$node_file = $ARGV[2];       #"nodelist.txt";
#@dm_scans = (0,44,75,100,119);     # good for 1900x1900x90int images
#@dm_scans = (0,70,98,119);
@dm_scans = (-1,0);                 # new format script interprets starting -1 to do all DMs
$user = `whoami`;chomp($user);
$working_dir = `pwd`;chomp($working_dir);

  # Acquire user password
  # !!! Password entry may not be necessary
  # !!! here if ssh tunneling is enabled (I
  # !!! think it is)
#$pass = get_password();



# * * * * * * * * * * * * * * * * * * *
# (0)  INIT DIRS + BASIC ERROR CHECK
# * * * * * * * * * * * * * * * * * * *

  # Test right number of command-line inputs
if ($#ARGV != 2){
    throw_error(0);
}

  # Test existence of source file and scan info
if (!-e "$asdm_file"){
    throw_error(1);
}
if (!-e "$asdm_file/Scan.xml"){
    throw_error(6);
}

  # Test if ms files already exist with requested basename
  # !!! This error testing is done by Casey but should I do it too?
  # !!! Also check for pkl files?
#$ms_base_exists = `find ./ -d 1 -name "$ms_base*" | wc -l | awk '{print \$1}'`;#chomp($ms_base_exists);
#if ($ms_base_exists){
#    throw_error(2);
#}

  # Test existence of valid node list
  # !!! Here is the fill point for autoquerying available nodes
  # !!! *** SHOULD ALSO DO A TEST HERE TO CHECK NODE ACTUALLY CONTACTABLE 
if (!-e "$node_file" || $node_file eq ''){
    throw_error(3);
}

  # Test that scripts exist where they need to
if (!-e "batch\_leanpipedt.py"){
    throw_error(4);
}

if (!-e "batch\_filler.py"){
    throw_error(7);
}

if (!-e "batch\_plot.py"){
    throw_error(8);
}

  # Create tracking directory if it doesn't exist
if (!-e "tracking_dir"){
    print "\nTracking directory tracking_dir does not exist; I will create it.\n";
    mkdir "tracking_dir";
}



# * * * * * * * * * * * * * * * * * * *
# (1) DETERMINE VALID PROCESSING NODES
# * * * * * * * * * * * * * * * * * * *
#
# !!! Here is the other filler point for autoquerying available nodes.
# !!! Currently this is done from text file. In the future this could
# !!! be done automatically but I don't know where to query.
#

  # Read node list
@node_list = ();
open(INF, "$node_file") || die "\n\tUNANTICIPATED ERROR\n\tCould not open file $node_file\n\n";
while(<INF>){
    chomp;
    if ($_ eq ""){next;}
    push(@node_list,$_);
}
close(INF);

  # Test N_nodes>0
$N_nodes = @node_list;
$N_failed_nodes = 0;
if ($N_nodes == 0){
    throw_error(5);
}



# * * * * * * * * * * * * * * * * * * *
#(2) DETERMINE SCAN, TARGET, and DM LIST
# * * * * * * * * * * * * * * * * * * *
#

  # Query scan list, taking only OBSERVE_TARGET scans.
$N_scans = 0;
$good_target = 0;
@scan_targets = ();
@scan_numbers = ();
$target_list = "";
$scan_holder = "";
open(SCANS, "$asdm_file/Scan.xml");
while(<SCANS>){
    chomp;

    # Grab the scan number in case we need it later
    if (/scanNumber/){
	($tmp,@tmp_a) = split('\>',$_);
	($scan_holder,$tmp) = split('\<',$tmp_a[0]);
    }

    # If we're on a good target, read name into scan name array
    if ($good_target && /sourceName/){
	($tmp,@tmp_a) = split('\>',$_);
	($source_name,$tmp) = split('\<',$tmp_a[0]);
	push (@scan_targets,$source_name);
	push (@scan_numbers,$scan_holder);
	#print "In scan $scan_holder($N_scans) found source $scan_targets[$N_scans-1]\n";
	$good_target = 0;
	if ($target_list !~ /$source_name/){
	    $target_list .= "\t$source_name\n";
	}
    }
    if (/OBSERVE_TARGET/){$N_scans++;$good_target = 1;}
}
close(SCANS);

  # Query DM list to check number of DM submissions, populate submission names.
$N_dms = @dm_scans;
$N_dms--;
@dm_range = ();
$current_dm=0;
for ($idm=0;$idm<$N_dms;$idm++){
    push (@dm_range, "$dm_scans[$idm],$dm_scans[$idm+1]");
}

print "\n\tFound $N_scans scans, $N_dms DM ranges, and $N_nodes valid nodes. Target list:\n$target_list\n";


  # Set up node-tracking arrays
@node_Nscans = ();
@node_info = ();
@node_scan = ();
@node_scans = ();
@node_dm = ();
for ($inode = 0;$inode<$N_noces;$inode++){
    push(@node_Nscans,0);
    push(@node_info,"");
    push(@node_scan,"");
    push(@node_dm,"");
    push(@node_scans,"");
}


# * * * * * * * * * * * * * * * * * * *
# (3)  CREATE ms AND pickle FILES
# * * * * * * * * * * * * * * * * * * *
# 
  # Determine number of scans and scans per node                                                                                        

print "Submitting pickling scripts... please be patient.\n";

if ($N_scans < $N_nodes){
    $N_per_node = 1;
    $N_spares = 0;
} else {
    $N_per_node = int($N_scans/$N_nodes);
    $N_spares = $N_scans % $N_nodes;
}

#@node_Nscans = ();
#@node_info = ();
$current_scan = 0;
print "Submitting scan: ";
for ($inode=0;$inode<$N_nodes;$inode++){
    $dt = localtime();
    $comma_scans = "";
    for ($iscan=$current_scan;$iscan<$current_scan+$N_per_node;$iscan++){
	$cmd = "batch\_filler.py $asdm_file $ms_base $scan_numbers[$iscan]";
	print "$iscan (obs scan index $scan_numbers[$iscan])";
	system "ssh -l $user $node_list[$inode] cd $working_dir\\;nohup python $cmd `</dev/null` \\>nohup-$ms_base-$node_list[$inode]-scan$scan_numbers[$iscan].$user 2\\>\\&1 \\& >& /dev/null";
#	system "ssh -l $user $node_list[$inode] cd $working_dir\\;nohup casapy --nologfile --nogui -c \\\"$cmd\\\" `</dev/null` \\>nohup-$ms_base-$node_list[$inode]-scan$iscan.$user 2\\>\\&1 \\& >& /dev/null";
#with logger	system "ssh -l $user $node_list[$inode] cd $working_dir\\;nohup casapy --nogui --logfile casapy-$ms_base-$node_list[$inode]-scan$iscan.log -c \\\"$cmd\\\" `</dev/null` \\>nohup-$ms_base-$node_list[$inode]-scan$iscan.$user 2\\>\\&1 \\& >& /dev/null";
	$comma_scans .= "$iscan,";
	sleep 1;
    }
    $current_scan = $iscan;

    $node_Nscans[$inode] = $N_per_node;
    if ($N_spares > 0){
	$cmd = "batch\_filler.py $asdm_file $ms_base $scan_numbers[$current_scan]";
	print "$current_scan (obs scan index $scan_numbers[$current_scan])";
	system "ssh -l $user $node_list[$inode] cd $working_dir\\;nohup python $cmd `</dev/null` \\>nohup-$ms_base-$node_list[$inode]-scan$scan_numbers[$current_scan].$user 2\\>\\&1 \\& >& /dev/null";
#	system "ssh -l $user $node_list[$inode] cd $working_dir\\;nohup casapy --nologfile --nogui -c \\\"$cmd\\\" `</dev/null` \\>nohup-$ms_base-$node_list[$inode]-scan$current_scan.$user 2\\>\\&1 \\& >& /dev/null";
#with logger	system "ssh -l $user $node_list[$inode] cd $working_dir\\;nohup casapy --nogui --logfile casapy-$ms_base-$node_list[$inode]-scan$current_scan.log -c \\\"$cmd\\\" `</dev/null` \\>nohup-$ms_base-$node_list[$inode]-scan$current_scan.$user 2\\>\\&1 \\& >& /dev/null";
	$comma_scans .= "$current_scan,";
	$N_spares = $N_spares - 1;
	$node_Nscans[$inode]++;
	$current_scan++;
	sleep 1;
    }
    chop($comma_scans);
    $node_scans[$inode] = $comma_scans;
    $node_info[$inode] = "Started pickling scans $comma_scans on $dt. Still $node_Nscans[$inode] left to finish.";
}

print "\n\n\n* * * * * Current processing status * * * * *\n* \n";
for ($inode = 0;$inode<$N_nodes;$inode++){
    print "* $node_list[$inode]: $node_info[$inode]\n";
}
print "* \n* * * * * * * * * * * * * * * * * * * * * * * *\n";




# * * * * * * * * * * * * * * * * * * *
# (3)   SUBMIT FIRST ROUND OF SCANS
# * * * * * * * * * * * * * * * * * * *
# 
# !!! Script submits one scan per node sequentially across the
# !!! nodes. As those scans finish, we can submit another scan to the
# !!! completed node.
#
# !!! SSH submission/tracking might be better/easier using
# !!! Net::SSH::Perl but for now let's try the simpler version that
# !!! doesn't require installing extra perl extensions. Below is a
# !!! commented-out example on how to get error reporting and
# !!! stdout/err through ssh, followed by the current solution.
# 

  # Set up commands and submit them thru ssh
#@node_info = ();
#@node_scan = ();
#@node_dm = ();
#$current_scan = 0;
#for ($inode=0;$inode<$N_nodes;$inode++){
#
#    # Create command for this node and submit it, noting date.
#    $cmd = "batch\_leanpipedt.py $asdm_file $ms_base $scan_numbers[$current_scan] $dm_range[$current_dm]";
#
#    # Set up node submit info for later user review
#    $dt = localtime();
#    $node_info[$inode] = "Started scan $current_scan->DM$dm_range[$current_dm] \[$scan_targets[$current_scan]\] on $dt";
#    $node_scan[$inode] = $current_scan;
#    $node_dm[$inode] = $dm_range[$current_dm];
#
#    # Submit jobs through ssh
#    system "ssh -l $user $node_list[$inode] cd $working_dir\\;nohup casapy --nogui --logfile casapy-$ms_base-$node_list[$inode]-scan$current_scan-dmrange$dm_range[$current_dm].log -c \\\"$cmd\\\" `</dev/null` \\>nohup-$ms_base-$node_list[$inode]-scan$current_scan-dmrange$dm_range[$current_dm].$user 2\\>\\&1 \\& >& /dev/null";
#    
#    print "$node_list[$inode]: $node_info[$inode]\n";
#    
#    # Give the previous job a bit of time to set up in working_dir
#    sleep 5;
#
#    $current_dm++;
#    if ($current_dm >= $N_dms){
#	$current_dm = 0;
#	$current_scan++;
#    }
#}
#
#$scans_completed = "---";
#print "\n\n* * * * * Current processing status * * * * *\n* Scans completed:\n* $scans_completed\n";
#for ($inode = 0;$inode<$N_nodes;$inode++){
#    print "* $node_list[$inode]: $node_info[$inode]\n";    
#}
#print "* \n* Processing up to scan $current_scan/$N_scans.\n";
#print "* * * * * * * * * * * * * * * * * * * * * * * *\n";

# * * * * * * * * * * * * * * * * * * *
# (3) REPORT PROGRESS AND SCAN FOR COMPLETED JOBS
# * * * * * * * * * * * * * * * * * * *
# 
# !!! *** ALSO NEED TO INPUT TESTS HERE FOR FAILURES.
#
$N_finished = 0;
$scans_completed = "";
$current_scan=0;
#while ($current_scan < $N_scans){
while ($N_finished < $N_nodes){

    $are_ready = `find tracking\_dir/ -name "*.ready*" | wc -l`;chomp($are_ready);$are_ready += 0;

    if ($are_ready){
	open(DONELIST, "find tracking\_dir/ -name \"*.ready*\"|");
	while(<DONELIST>){
	    chomp;
	    
	    ($tmp1,$tmp2) = split('\\.',$_);
	    ($tmp3,$host) = split('/',$tmp1);

	    # Find node in nodelist, report completion.
	    for ($inode=0;$inode<$N_nodes;$inode++){
		if ($node_list[$inode] eq $host){
	    	    # Remove notification file
		    system "rm -rf $_";
		    last;
		}
	    }

	    # If node is invalid
	    if ($inode>=$N_nodes){
		print "Skipping invalid .ready file, $_.\n\n";
		next;
	    }

	    $node_Nscans[$inode]--;
	    #print "It is $node_Nscans[$inode]\n";
	    if ($node_Nscans[$inode] <= 0){
		print "*** $node_list[$inode] finished pickling its scans and will now commence processing.\n";
		system "touch tracking\_dir/$node_list[$inode].finished";

		$dt = localtime();
		$node_info[$inode] = "Idle since $dt";
		#!!! NEED TO PUT IN A CHECK TO MAKE SURE SCAN'S PICKLE AND MS FILE EXIST before submitting jobs!
	    } else {
		$node_info[$inode] = "Started pickling scans $node_scans[$inode] on $dt. Still $node_Nscans[$inode] left to finish.";
	    }
	}
    }

    $are_done = `find tracking\_dir/ -name "*.finished" | wc -l`;chomp($are_done);$are_done += 0;
    
    if ($are_done){
	open(DONELIST, "find tracking\_dir/ -name \"*.finished\"|");
	while(<DONELIST>){
	    chomp;

	    ($tmp1,$tmp2) = split('\\.',$_);
	    ($tmp3,$host) = split('/',$tmp1);

	    # Find node in nodelist, report completion.
	    for ($inode=0;$inode<$N_nodes;$inode++){
		if ($node_list[$inode] eq $host){
		    # Remove notification file
		    system "rm -rf $_";
		    last;
		}
	    }

	    # If node is invalid
	    if ($inode>=$N_nodes){
		print "Skipping invalid .finished file, $_.\n\n";
		next;
	    }

	    # Is the following stupidity check really necessary? Trust in Casey...
	    #if ($inode >= $N_nodes){throw_error(7);
	    if ($node_scan[$inode] ne ""){
		$completed_target = $scan_targets[$node_scan[$inode]];
		$scans_completed .= " $node_scan[$inode]->$node_dm[$inode]";
		print "*** $node_list[$inode] finished scan $node_scan[$inode] DM range $node_dm[$inode] on target $completed_target.";
		$N_finished++;
	    }
	    $dt = localtime();
	    $node_info[$inode] = "Idle since $dt";

	    if ($current_scan < $N_scans){
		# Set up command and node submit info
		$cmd = "batch\_leanpipedt.py $asdm_file $ms_base $scan_numbers[$current_scan] --dms=$dm_range[$current_dm]";
		$node_info[$inode] = "Started scan $current_scan->DM$dm_range[$current_dm] \[$scan_targets[$current_scan]\] on $dt";
		$node_scan[$inode] = $current_scan;
		$node_dm[$inode] = $dm_range[$current_dm];	
		

		# Submit job through ssh
#		system "ssh -l $user $node_list[$inode] cd $working_dir\\;nohup echo $cmd `</dev/null` \\>nohup-$ms_base-$node_list[$inode]-scan$scan_numbers[$current_scan]-dmrange$dm_range[$current_dm].$user 2\\>\\&1 \\& >& /dev/null"; # test
		system "ssh -l $user $node_list[$inode] cd $working_dir\\;nohup python $cmd `</dev/null` \\>nohup-$ms_base-$node_list[$inode]-scan$scan_numbers[$current_scan]-dmrange$dm_range[$current_dm].$user 2\\>\\&1 \\& >& /dev/null";
#		system "ssh -l $user $node_list[$inode] cd $working_dir\\;nohup casapy --nologfile --nogui -c \\\"$cmd\\\" `</dev/null` \\>nohup-$ms_base-$node_list[$inode]-scan$current_scan-dmrange$dm_range[$current_dm].$user 2\\>\\&1 \\& >& /dev/null";
#with logger		system "ssh -l $user $node_list[$inode] cd $working_dir\\;nohup casapy --nogui --logfile casapy-$ms_base-$node_list[$inode]-scan$current_scan-dmrange$dm_range[$current_dm].log -c \\\"$cmd\\\" `</dev/null` \\>nohup-$ms_base-$node_list[$inode]-scan$current_scan-dmrange$dm_range[$current_dm].$user 2\\>\\&1 \\& >& /dev/null";
		#SSS system "ssh -l $user $node_list[$inode] cd $working_dir\\;nohup casapy --nogui --logfile casapy-$ms_base-$node_list[$inode]-scan$current_scan.log -c \\\"$cmd\\\" `</dev/null` \\>nohup-$ms_base-$node_list[$inode]-scan$current_scan.$user 2\\>\\&1 \\& >& /dev/null" == 0 
		#SSS || remove_node();
		# !!! *** Note remove_node does not currently function, as I am
		# !!! *** not properly calling system to execute the ssh command
		# !!! *** so cannot test its success or failure. Need to fix this.
		print " $node_info[$inode]\n";
		$N_finished--;
		
		# Give the previous job a bit of time to set up in working_dir
		sleep 1;

		$current_dm++;
		if ($current_dm >= $N_dms){
		    $current_dm = 0;
		    $current_scan++;
		    if ($current_scan >= $N_scans){last;}
		}
	    }
	}
	close(DONELIST);
    }

    if ($are_done || $are_ready){
	print "\n\n* * * * * Current processing status * * * * *\n* Scans completed:\n* $scans_completed\n";
	for ($inode = 0;$inode<$N_nodes;$inode++){
	    print "* $node_list[$inode]: $node_info[$inode]\n";    
	}
	print "* \n* Processing up to scan $current_scan/$N_scans.\n";
	print "* * * * * * * * * * * * * * * * * * * * * * * *\n";
    }
    
    sleep 2;
}


# * * * * * * * * * * * * * * * * * * *
# (5)    FIND AND REPORT FAILURES
# * * * * * * * * * * * * * * * * * * *
# 
# * * * * * * * * * * * * * * * * * * *
# (6)    COLLECT/ORGANIZE RESULTS
# * * * * * * * * * * * * * * * * * * *
# 

print "\n\nPROCESSING HAS COMPLETED. Now producing plots for this batch job...\n";
system "python batch\_plot.py $ms_base";
print "\n\n\n\tTHIS BATCH JOB HAS COMPLETED.\n\n\tTime for a cup of tea.\n\n\n";


#
#
# * * * * * * * * * * * * * * * * * * *
#       S U B F U N C T I O N S
# * * * * * * * * * * * * * * * * * * *
#

# Reminder: http://www.asciitable.com/ for ord keys
sub get_password(){
    my $key = 0;
    my $password = "";  
    
    print "\nPlease input password for $user: ";   
    
    # Start reading the keys with 'no type' for control keys
    # Go until enter key is pressed [ ord(10) ]
    ReadMode(4);
    while(ord($key = ReadKey(0)) != 10){
	# Chop if del/backspace and replace * with ' '
	if(ord($key) == 127 || ord($key) == 8){
	    chop($password);
	    print "\b \b";
	} elsif(ord($key) >= 32) {  
	    $password = $password.$key;   
	}
    }
    ReadMode(0);

    return($password);
}


sub remove_node(){
    $node_info[$inode] = "*** NODE UNAVAILABLE; ssh submission failed on $dt";
    $N_failed_nodes++;
    $current_scan--;
    return;
}


sub throw_error(){
    my $n_in = @_;
    my $err = '';
    my $mode = '';
    if ($n_in == 1){
	$err = shift(@_);
    } elsif ($n_in == 2){
	($err,$mode) = @_;
    } else {
	# Do something if bad call to throw_error;
	print "\n\tERROR:\n\tWeird internal call to throw_error.\n";
	exit(999);
    }

    if ($err == 0){
	print STDERR "\n\tERROR:\n\tCommand line inputs should be:\n\tdistribute <ASDM file name> <ms output basename> <nodelist>\n\n";
    } elsif ($err == 1){
	print STDERR "\n\tERROR:\n\tCould not find the asdm file you provided ($asdm_file)\n\n";
    } elsif ($err == 2){
	print STDERR "\n\tERROR:\n\tFiles with basename $ms_base already exist.\n\tPlease remove these and try running again.\n\n";
    } elsif ($err == 3){
	print STDERR "\n\tERROR:\n\tEssential file nodelist.txt does not exist; please\n\tcreate this text file with one node name per line.\n\n";
    } elsif ($err == 4){
	print STDERR "\n\tERROR:\n\tCould not find batch\_leanpipedt.py; are you sure\n\tyou're running this script in the right directory?\n\n";
    } elsif ($err == 5){
	print STDERR "\n\tERROR:\n\tNo valid nodes in nodelist.txt; please create\n\tthis text file with one node name per line.\n\n";
    } elsif ($err == 6){
	print STDERR "\n\tERROR:\n\tCould not find Scan.xml file in the ASDM file\n\tcalled $asdm_file ; please ensure this exists\n\n";
    } elsif ($err == 7){
	print STDERR "\n\tERROR:\n\tCould not find batch\_filler.py; please make\n\tit exists in this directory.\n\t (or, casey, I can call the code dir version directly :)! )\n\n";
    } elsif ($err == 8){
	print STDERR "\n\tERROR:\n\tCould not find batch\_plot.py; please ensure\n\tit exists in this direcotry.\n\t (or, casey, I can call the code dir version directly :)! )\n\n";
    } elsif ($err == 9){
    } elsif ($err == 10){
    } elsif ($err == 11){
    }
    exit($err);
}


# !!! HERE IS A Net::SSH:Perl SNIPPET
#$host = "nmpost-01";
#my $password = "password";
#
##-- set up a new connection
#my $ssh = Net::SSH::Perl->new($host);
##-- authenticate
#$ssh->login($user, $pass);
##-- execute the command
#my($stdout, $stderr, $exit) = $ssh->cmd("$RUN_COMMAND");
#
