#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Getopt::Long;

my $dagchainer_file;
my $genelist_file;

GetOptions("dagchainer=s" => \$dagchainer_file,
    "genelist=s" => \$genelist_file);

unless (defined $dagchainer_file && defined $genelist_file) {
print<<EOF;
Usage:
    kudo_build_syntelog_table
        --dagchainer    [-d]   dagchainer output file
        --genelist      [-g]   ath gene list

EOF
exit 1;
}



open(ATH,"<".$genelist_file) || die "Cannot open file\n";
open(DG,"<".$dagchainer_file) || die "Cannot open file\n";
open(F,">".$dagchainer_file.".failed.colinear") || die "Cannot write file\n";
open(C,">".$dagchainer_file.".success.colinear") || die "Cannot write file\n";
open(ACh,">".$dagchainer_file.".all.chains") || die "Cannot write file\n";
open(PCh,">".$dagchainer_file.".chains.passed") || die "Cannot write file\n";

my %dat;
my %chain_members;
my %genes_of_interest;
my $chain_name ='';
my $chain_score=0;
my $num_pairs  =0;
my $last_chr='';
my $gc = 0;
my @queries = ('C1','C1.r','C2','C2.r','C3','C3.r','C4','C4.r','C5','C5.r','C6','C6.r','C7','C7.r','C8','C8.r','C9','C9.r',
    'Un','Un.r','overlap');


while(<ATH>) {
    chomp;
    my ($g,$gm,$length,$chr_ord,$chr,$ord,$g_s,$g_e,$j,$desc)=(split("\t",$_));
    if($last_chr ne $chr) {
        $gc = 0;
    }
    $genes_of_interest{$chr}{$g}=$gc++;
    $last_chr = $chr;
}

while(<DG>) {
    chomp;
    if ($_ =~ /^#/) {  ##header lines
        my $tmp = $_;
        my $q ='';
        my $r ='';
        my $rank=0;
        my $rev='';

       if($_ =~ /reverse/) {
            $tmp =~ m/alignment\ (\S*)\ vs.\ (\S*)\ \(reverse\)\ Alignment\ #(\d*)\ \ score\ =\ (\d*.\d*)\ \(num\ aligned\ pairs:\ (\d*)\):/;
            ($q,$r,$rank,$chain_score,$num_pairs) = ($1,$2,$3,$4,$5);
            $rev='.r';
        }
        else {
            $tmp =~ m/alignment\ (\S*)\ vs.\ (\S*)\ Alignment\ #(\d*)\ \ score\ =\ (\d*.\d*)\ \(num\ aligned\ pairs:\ (\d*)\):/;
            ($q,$r,$rank,$chain_score,$num_pairs) = ($1,$2,$3,$4,$5);
        }
        if($r =~ /^Chr/) {
            my $tmp = $r;
            $r = $q;
            $q = $tmp;
        }
        $r = 'Un' if ($r =~ /^utg/ || $r =~ /^Scaffold/);
        $chain_name = join("_",$r.$rev,$q,$rank);
    }
    else {  ## align-pairs
        my ($chr1,$g1,$g1_s,$g1_e,$chr2,$g2,$g2_s,$g2_e,$eval,$score) = (split("\t",$_));
        if ($chr2 =~ /^Chr/) {
            $g2 =~ s/\.[0-9]*$//g;
            $dat{$chr2}{$chain_name}{$g2}=$g1;
            $chain_members{$g1}{$chain_name}=$chain_score;
        }
        else {
            $g1 =~ s/\.[0-9]*$//g;
            $dat{$chr1}{$chain_name}{$g1}=$g2;
            $chain_members{$g2}{$chain_name}=$chain_score;
        }
        print ACh join("\t",$chain_name,$_)."\n";
    }
}

warn "Preprocessing chains ...\n";
my %black_list;


## black list orthologous pairs with lower chain score
foreach my $bog (keys %chain_members) {
    my @memberships = (keys %{$chain_members{$bog}});
    my $num_memberships = scalar @memberships;
    next unless ($num_memberships > 1);
    @memberships = (sort {$chain_members{$bog}{$b} <=> $chain_members{$bog}{$a} } @memberships);
    for(my $i=0; $i<$num_memberships;$i++) {
        $black_list{$memberships[$i]}{$bog}=1 if ($i>0);
#        print join("\t",$bog,$memberships[$i],$chain_members{$bog}{$memberships[$i]})."\n";
    }
}

# drop black listed pairs
# if a chain drop below 6 or few members drop the chain
foreach my $chr (sort {$a cmp $b} keys %dat) {
    warn "    chromosome $chr ... \n";

    # sort chains by # of genes in descending order
    my @chains = (sort {scalar keys %{$dat{$chr}{$b}} <=> scalar keys %{$dat{$chr}{$a}} } keys %{$dat{$chr}});
    warn "  found ".scalar @chains." initial chains\n";
    foreach my $ch (@chains) {
        foreach my $g (keys %{$dat{$chr}{$ch}}) {
            my $ortho = $dat{$chr}{$ch}{$g};
            delete $dat{$chr}{$ch}{$g} if (exists($black_list{$ch}{$ortho}));
        }
        if(scalar (keys %{$dat{$chr}{$ch}}) <= 6) {          ## check chain size after filtering
            warn "  delete $ch\n";
            delete $dat{$chr}{$ch};
        }
        else {
            print PCh $ch."\n";
            
        }
    }
}

warn "Done Pre-processing\n";

my %tracks;
foreach my $chr (sort {$a cmp $b} keys %dat) {
    warn "Processing chromosome $chr ... \n";
    # sort chains by # of genes in descending order
    my @chains = (sort {scalar keys %{$dat{$chr}{$b}} <=> scalar keys %{$dat{$chr}{$a}} } keys %{$dat{$chr}});
    warn "    found ".scalar @chains." chains\n";

    # consider synteny in Ath genes of interest
    my @genes = (sort {$a cmp $b} keys %{$genes_of_interest{$chr}});

    foreach my $ch (@chains) {
        warn $chr."\t".$ch."\t".scalar (keys %{$dat{$chr}{$ch}})."\n";
        my @ch_gs = (sort {$a cmp $b} keys %{$dat{$chr}{$ch}});
        my ($q,$r,$rank) = (split("_",$ch));
        
        my @overlap = grep  { exists($tracks{$chr}{$q}{$_}) } @ch_gs;
        if(scalar @overlap) {
            my $overlap_count = 0;
            foreach my $gene (@ch_gs) {
                if (exists $tracks{$chr}{$q}{$gene}) {
                    $overlap_count++;
                }
            }
            if ($overlap_count / scalar(@ch_gs) > 0) {
                foreach my $g (@ch_gs) {
                    $tracks{$chr}{'overlap'}{$g} = $dat{$chr}{$ch}{$g};
                }
                my $s = $genes_of_interest{$chr}{$ch_gs[0]};
                my $e = $genes_of_interest{$chr}{$ch_gs[-1]};
                for(my $i=$s; $i<=$e; $i++) {
                    my $ortho = $dat{$chr}{$ch}{$genes[$i]};
                    $ortho = 'x' unless defined $ortho;
                    print F join("\t",$ch,$chr,$genes[$i],$ortho)."\n";
                }
            }
            else {
                foreach my $g (@ch_gs) {
                    $tracks{$chr}{$q}{$g} = $dat{$chr}{$ch}{$g};
                }
                foreach my $g (@ch_gs) {
                    $tracks{$chr}{'overlap'}{$g} = $dat{$chr}{$ch}{$g};
                }
            }
        }
        else {
            foreach my $g (@ch_gs) {
                $tracks{$chr}{$q}{$g} = $dat{$chr}{$ch}{$g};
            }
        }
    }
}


print C join("\t",'locus_id',@queries)."\n";
foreach my $chr (sort {$a cmp $b} keys %tracks) {
    my @genes = (sort {$a cmp $b} keys %{$genes_of_interest{$chr}});
    foreach my $g (@genes) {
        print C $g;
        foreach my $q (@queries) {
            my $ortho = $tracks{$chr}{$q}{$g};
            $ortho='x' unless defined $ortho;
            print C "\t".$ortho;
        }
        print C "\n";
    }

}

