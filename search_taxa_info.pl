#! /usr/bin/perl

=pod

=head1 search_taxa_info.pl (written by Sishuo Wang at University of British Columbia)

searching for taxonomic information of species in multiple ways based on NCBI Taxanomy database (http://www.ncbi.nlm.nih.gov/taxonomy)

#	Please note that the taxonomic information retrieved from NCBI Taxanomy database can be incorrect, so manual check may be needed depending on your requests.

Also many thanks to the help from Fabien Burki from University of British Columbia who contributed greatly to the part of the script in terms of looking for taxonomic information via Bioperl.

Please do not hesitate to write e-mail to tomassonwss@gmail.com or sishuowang@hotmail.ca if you have any questions and/or suggestions of the script. Your help will be highly appreciated.

=head1 usage

perl search_taxa_info.pl <-taxa taxaname> [-db_file database file] [-case case insensitive search] [-fuzzy perform fuzzy search] [-rank taxonomic level] [-concise] [-Bio_DB_Taxonomy|-Bioperl] [-help]
 
Note:
 0.	You could search for multiple taxa meanwhile by specifiying '-taxa'. (e.g. "perl search_taxa_info.pl -taxa Homo -taxa Drosophila -taxa Arabidopsis").
 1.	C<-db_file> is recommended to be specified. If the database file is not given, then the script will look for the file taxa_info.db in the current folder.
 2.	If C<-case> is specified, you will be able to perform the search in a case-insensitive pattern (e.g., the taxa_name "homo sapiens" and Homo sapiens will give the same result if C<-case> is specified).
 3.	If you want to perform fuzzy search (e.g., use the keyword "Saccharomyces cerevi" to search for "Saccharomyces cerevisiae" and all other taxa whose names contain the keyword), please specify C<-fuzzy>.
 4.	If there is any space in the taxaname, please use quotation mark (either double or single is OK). For instance, it is recommended to type "perl search_taxa_info.pl -taxa "Homo sapiens" taxa_info.db" where the word "Homo sapiens" is surrounded by double quotation mark.
 5.	Legal taxonomic levels are superkindom, kingdom, phylum, order, class, family, genus, species. If you would like to include all taxanomic levels, it is recommended to use '-rank all'.
 6.	"-concise":	to output the result in a concise format (see examples below).
 	only active when "-rank" is specified.
 7.	"-Bio_DB_Taxonomy or -Bioperl":	If no results are found in the sqlite3 database, then perform searching via Bioperl (The modules Bio::SearchIO and Bio::DB:Taxonomy should be installed).

=head1 notes

	Please make sure that you put the file search_taxa_info.pl and taxa_info.db are ready before starting to perform the search.
	Please make sure that Perl 5.010 (or later version) and SQLite has been isntalled on you computer.

=head1 an example of the result

You may want to try the following commands to have a test:
 1.	B<perl search_taxa_info.pl -taxa C<Saccharomyces cerevisiae> -db_file taxa_info.db>
 2.	B<perl search_taxa_info.pl -taxa C<Saccharomyces cerevisiae> -db_file taxa_info.db -rank genus -rank family>
 3.	B<perl search_taxa_info.pl -taxa C<Saccharomyces cerevisiae> -db_file taxa_info.db -rank genus -rank family -concise>

You are expected to get the following results respectively.

1.
##########################################################################################
 Saccharomyces cerevisiae
       species	      Saccharomyces cerevisiae
         genus	                 Saccharomyces
        family	            Saccharomycetaceae
         order	             Saccharomycetales
         class	               Saccharomycetes
        phylum	                    Ascomycota
       kingdom	                         Fungi
  superkingdom	                     Eukaryota

2.
Saccharomyces cerevisiae
family:	Saccharomycetaceae
genus:	Saccharomyces

3.
Saccharomyces cerevisiae	Saccharomycetaceae	Saccharomyces

=cut


################################################################################
use strict;
use DBI;
use Getopt::Long;
use File::Basename;
use 5.010;


my ($search_id_href, $db_file, $fuzzy, $case, $rank_swi, $rank_select_href, $concise, $Bio_DB_Taxonomy_swi);
my ($db, %key_order, %lineage, $focus);
my @key_order = qw(species genus family order class phylum kingdom superkingdom);
my %legal_fuzzy_way = (
	'before'     => 1,
	'after'      => 1,
	'everywhere' => 1,
);

my @stdModules = ("Bio::DB::Taxonomy", "Bio::SearchIO");
@key_order{@key_order} = (1) x 8;

scalar @ARGV < 1 and &show_usage;
($search_id_href, $db_file, $fuzzy, $case, $rank_swi, $rank_select_href, $concise, $Bio_DB_Taxonomy_swi) = &read_and_check_param(@ARGV);
$db_file=dirname($0) . '/' . "taxa_info.db" if not $db_file;


###################################################################################
my ($sth, $rv);
my $dbargs = {AutoCommit => 0, PrintError => 1};
my $dbh = DBI->connect("dbi:SQLite:dbname=$db_file","","",{AutoCommit => 0});
$dbh->commit();

foreach my $search_id (keys %{$search_id_href}){
	&get_taxid('', '', $search_id, '');
	if (not exists $lineage{$search_id}){
		&look_for_syn($search_id);
	}
	if (not exists $lineage{$search_id}){
		if ($Bio_DB_Taxonomy_swi){
			&look_for_Bio_DB_Taxonomy($search_id, \@stdModules);
		}
	}
}

&output();



#####################################################################################
sub look_for_Bio_DB_Taxonomy{
	my $search_id = shift();
	my $stdModules_aref = shift();
	@stdModules = @$stdModules_aref;

	my $return_check_modules = &check_modules();
	(not $return_check_modules) and return();

	my $db = new Bio::DB::Taxonomy(-source => 'entrez');
	my $taxonid = $db->get_taxonid("$search_id");
	my @nodes= $db->get_Taxonomy_Node(-taxonid => $taxonid);
	foreach my $node (@nodes){
		#print OUT $node->scientific_name . "\n";
		#print $node->scientific_name ."\n";
		my $parent = $db->get_Taxonomy_Node($node->parent_id);
		# print $parent->node_name, "\t";
		while (defined $parent and $parent->node_name ne 'root'){
			my $rank = $parent->rank;
			if (exists $key_order{$rank}){
				my $parent_node_name = $parent->node_name;
				$lineage{$focus}{$parent_node_name} = 1;
			}
			$parent = $db->get_Taxonomy_Node($parent->parent_id);
		}
	}

#===============================================================================#
	sub check_modules{
		my %notFoundModules=();
		foreach my $module (@stdModules){
			eval("use $module;");
			if ($@){
				$notFoundModules{$module}=1;
			}
		}
		if (keys %notFoundModules){
			print "The following modules which are necessary if you want to perform searching in the NCBI taxonomy database using bioperl may not be corretcly installed:\n";
			map {print $_."\n"} keys %notFoundModules;
			print "\n";
			return(0);
		}
		return(1);
	}
}


sub look_for_syn{
	my $search_id = shift();
	my $dirname = dirname($0);
	open (IN,'<',"$dirname/additional/syn.tab");
	while(<IN>){
	chomp;
		my ($taxid, $name) =split /\t/;
		if ($name eq $search_id){
			my $prepare_cmd = "select * from taxa where taxid = \"$taxid\"";
			$sth = $dbh->prepare("$prepare_cmd");
			$rv = $sth->execute() or die $DBI::errstr;

			my @row_aref_array = &fetchrow_array($sth);
			foreach (@row_aref_array){
				my ($taxid, $ptaxid, $name, $rank) = @$_;
				if ($ptaxid){
				$focus=$search_id;
					&get_taxid('',$ptaxid);
				}
			}
			last;
		}
	}
	close IN;
}


sub output{
	#foreach my $key1 (keys %$search_id_href){
	if (not $fuzzy and not $case){
		map {print join("\t", $_, "NotFound!\n")} grep {not exists $lineage{$_}} keys %$search_id_href;
		}
	for my $key1 (keys %lineage){
		if (not $rank_swi){
			printf '###' x 30 . "\n";
			print "taxaname:\t$key1\n";
			foreach my $rank (@key_order){
				printf "%15s",$rank."\t";
				printf "%30s", $lineage{$key1}{$rank} if exists $lineage{$key1}{$rank};
				print "\n";
			}
			print "\n";
		}
		else{
			$concise ? do {print "$key1\t"} : do {print "$key1\n"};
			foreach my $rank (reverse @key_order){
				next if not exists $rank_select_href->{$rank};
				my $rank_2;
				if (exists $lineage{$key1}{$rank}){
					$rank_2 = $lineage{$key1}{$rank};
				}
				else{
					$rank_2='None';
				}
				if (not $concise){
					print "$rank:\t$rank_2\n";
				}
				else{
					print "$rank_2\t";
				}
			}
			print "\n";
		}
	}
}


sub get_taxid{
	my ($search_item, $search_column, @search_item, $case_cmd);
	my ($prepare_cmd, $execute_cmd, $row_aref);
	my @search_column = qw(taxid ptaxid name rank);
	my @argu = @_;
	foreach (0..$#search_column){
		if ($argu[$_]){
			given($_){
				when($_ == 1)	{$search_column='taxid'}
				when($_ == 2)	{$search_column='name'}
			}
			$search_item = $argu[$_];
			last;
		}
	}

	if ($case){
		$case_cmd='COLLATE NOCASE';
	}

	if ($fuzzy){
		if ($search_column eq 'name'){
			my $search_item_fuzzy = "$search_item";
			if ($fuzzy eq 'before'){
				$search_item_fuzzy = "\"$search_item%\"";
			}
			elsif ($fuzzy eq 'after'){
				$search_item_fuzzy = "\"%$search_item\"";
			}
			else{
				$search_item_fuzzy = "\"%$search_item%\"";
			}
			$prepare_cmd = "select * from taxa where $search_column like $search_item_fuzzy $case_cmd";
			#$prepare_cmd = "select * from taxa where $search_column like \"%$search_item%\" $case_cmd";
		}
		else{
			$prepare_cmd = "select * from taxa where $search_column == \"$search_item\" $case_cmd";
		}
	}
	else{
		$prepare_cmd = "select * from taxa where $search_column = \"$search_item\" $case_cmd";
	}
	
	#print "$prepare_cmd\n";
	$sth = $dbh->prepare("$prepare_cmd");
	$rv = $sth->execute() or die $DBI::errstr;

	my @row_aref_array = &fetchrow_array($sth);
	foreach (@row_aref_array){
		return() if not $_;
		my ($taxid, $ptaxid, $name, $rank) = @$_;
		$focus = $name if $search_column eq 'name';
		#next if $rank eq 'no rank';
		$lineage{$focus}{$rank} = $name;
		#print $name."\t".$rank."\n";
		if ($ptaxid){
			&get_taxid('',$ptaxid);
		}
		else{
			return();
		}
	}
	return(1);
}


####################################################################################
sub fetchrow_array{
	my @row_aref_array;
	my ($sth) = @_;
	while(my @row = $sth->fetchrow_array()){
		my ($taxid, $ptaxid, $name, $rank) = @row;
		return(0) if $ptaxid == 0;
		push @row_aref_array, [@row];
	}
	return(@row_aref_array);
}


sub show_usage{
	system ("perldoc $0");
	exit 1;
}


sub read_and_check_param(){
my (%search_id, $db_file, $fuzzy, $case, $rank_swi, %rank_select, $concise, $Bio_DB_Taxonomy_swi);
while(@_){
	$_ = shift;
	given ($_){
		when(/-?-taxa/)		{$search_id{shift()}=1}
		when(/-?-db_file/)	{$db_file=shift}
		when(/-?-rank/)		{$rank_swi=1;
                             #my @ranks = shift().split(",");
                             #map {$rank_select{$_}=1} @ranks;
                             $rank_select{shift()}=1;
                            }
		when(/-?-concise/)	{$concise=1}
		when(/-?-case/)		{$case=1}
		when(/-?-fuzzy/)	{$fuzzy=1;
					my $input_fuzzy_way = shift;
					if ($input_fuzzy_way =~ /\-/){unshift(@_, $input_fuzzy_way)}
					else{
						if (exists $legal_fuzzy_way{$input_fuzzy_way}) {$fuzzy = $input_fuzzy_way;}
						else {$fuzzy = 'everywhere';}
					}
		}
		when(/-?-(Bio_DB_Taxonomy|
			  Bioperl)
			/xi)		{$Bio_DB_Taxonomy_swi=1}
		when(/-?-(h|help)/)	{&show_usage()}
		default			{print "param\t$_\tis illegal! Exit ......\n"; exit 1}
	}
}
&check_rank(\%rank_select);
&check_taxa(\%search_id);
return(\%search_id, $db_file, $fuzzy, $case, $rank_swi, \%rank_select, $concise, $Bio_DB_Taxonomy_swi);


#===============================================================================#
sub check_taxa{
	my ($search_id_href) = @_;
	foreach (keys %$search_id_href){
		if (not $_){
			print "-taxa does not have a legal value!\n";
			&show_usage()
		}
	}
}


sub check_rank{
	my ($rank_select_href) = @_;
	foreach (keys %{$rank_select_href}){
		if ($_ eq 'all'){
			map {$rank_select_href->{$_}=1} @key_order;
			return ();
		}
		else{
			print "Note: $_ is not a legal taxonomic unit!\n\n" and &show_usage() if not exists $key_order{$_};
		}
	}
}
}


