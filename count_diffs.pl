#!/usr/bin/perl

$|++;

# Copyright (C) 2017, Joseph F. Ryan <joseph.ryan@whitney.ufl.edu>
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# for more info see <http://www.gnu.org/licenses/>

use strict;
use warnings;
use Fasta;
use Data::Dumper;

our $DIR = '/home/aohdera/akidata/CassGenome_Dec2016/Current/WhitneyAssemblies/03-Alatina/07-CNIDARIA-ORTHOFINDER/Results_Nov09/Sequences';
our $OUTFILE = 'output.txt';

our $VERSION = 0.01;

#Cassiopea_xamachana.v1_* = Cassiopea gene models
#Cxam.HS.Trinity.transdecoder= CxamHS transcriptome
#Cxam.Scyph.transdecoder= CxamScyph transcriptome
#Calvadosia_cruxmelitensis_Trinity_SeqPrep-sickle.fasta.transdecoder.pep.A_Ccrux.Trinity.* = Calvadosia transcriptome
#Calvadosia_cruxmelitensis.v3_* = Calvadosia gene models
#Alatina_alata.masurca.v3_* = Alatina gene models
#Alatina_alata.transdecoder_Alatina_* = Alatina transcriptome
#Hmag_Hm_XP_* = Hydra gene model
#Nvec_jgi|Nemve1|145823|e_gw.* = Nematostella gene model
#HumRef2017_ = human
#xam.HS_2

MAIN: {
    opendir DIR, $DIR or die "cannot opendir $DIR:$!";
    my @seq_files = grep { /fa$/ } readdir DIR;
    my %ogs = ();
    FA: foreach my $sf (@seq_files) {
        my $fp = JFR::Fasta->new("$DIR/$sf");
        my %spe = ();
        while (my $rec = $fp->get_record()) {
            $rec->{'def'} =~ m/^>(([^_]+).*)/;
            my $sp = get_sp($2);
            $spe{$sp}++;
        }
	#all        
        if (($spe{'cxg'} || $spe{'cxh'}
        || $spe{'cxs'}) && ($spe{'ccg'} || $spe{'cct'}) && ($spe{'aag'} || $spe{'aat'}) && $spe{'hmg'} && $spe{'nvg'} && $spe{'hum'}){
            push @{$ogs{'all'}}, $sf;
            #Cnidarian specific
        } elsif (($spe{'cxg'} || $spe{'cxh'}
        || $spe{'cxs'}) && ($spe{'ccg'} || $spe{'cct'}) && ($spe{'aag'} || $spe{'aat'}) && $spe{'hmg'} && $spe{'nvg'} && !$spe{'hum'}) {
            push @{$ogs{'Cnidaria:cxgy_cxhy_cxsy_ccgy_ccty_aagy_aaty_hmgy_nvgy_humn'}}, $sf;
            #Cassiopea specific (w/ human)
        } elsif (($spe{'cxg'} || $spe{'cxh'}
        || $spe{'cxs'}) && !($spe{'ccg'} || $spe{'cct'}) && !($spe{'aag'} || $spe{'aat'}) && !$spe{'hmg'} && !$spe{'nvg'} && $spe{'hum'}) {
            push @{$ogs{'Cxam_w/Hum:cxgy_cxhy_cxsy_ccgn_cctn_aagn_aatn_hmgn_nvgn_humy'}}, $sf;
            #Cassiopea specific (no human)
        } elsif (($spe{'cxg'} || $spe{'cxh'}
        || $spe{'cxs'}) && (!$spe{'ccg'} && !$spe{'cct'}) && (!$spe{'aag'} && !$spe{'aat'}) && !$spe{'hmg'} && !$spe{'nvg'} && !$spe{'hum'}) {
            push @{$ogs{'Cxam:cxgy_cxhy_cxsy_ccgn_cctn_aagn_aatn_hmgn_nvgn_humn'}}, $sf;
            #Calvadosia specific (w/ human)
        } elsif (!$spe{'cxg'} && !$spe{'cxh'}
        && !$spe{'cxs'} && ($spe{'ccg'} || $spe{'cct'}) && !$spe{'aag'} && !$spe{'aat'} && !$spe{'hmg'} && !$spe{'nvg'} && $spe{'hum'}) {
            push @{$ogs{'Ccrux_w/Hum:cxgn_cxhn_cxsn_ccgy_ccty_aagn_aatn_hmgn_nvgn_humy'}}, $sf;
            #Calvadosia specific (no human)
        } elsif (!$spe{'cxg'} && !$spe{'cxh'}
        && !$spe{'cxs'} && ($spe{'ccg'} || $spe{'cct'}) && !$spe{'aag'} && !$spe{'aat'} && !$spe{'hmg'} && !$spe{'nvg'} && !$spe{'hum'}) {
            push @{$ogs{'Ccrux:cxgn_cxhn_cxsn_ccgy_ccty_aagn_aatn_hmgn_nvgn_humn'}}, $sf;
            #Alatina specific (w/ human)
        } elsif (!$spe{'cxg'} && !$spe{'cxh'}
        && !$spe{'cxs'} && !$spe{'ccg'} && !$spe{'cct'} && ($spe{'aag'} || $spe{'aat'}) && !$spe{'hmg'} && !$spe{'nvg'} && $spe{'hum'}) {
            push @{$ogs{'Aala_w/Hum:cxgn_cxhn_cxsn_ccgn_cctn_aagy_aaty_hmgn_nvgn_humy'}}, $sf;
            #Alatina specific (no human)
        } elsif (!$spe{'cxg'} && !$spe{'cxh'}
        && !$spe{'cxs'} && !$spe{'ccg'} && !$spe{'cct'} && ($spe{'aag'} || $spe{'aat'}) && !$spe{'hmg'} && !$spe{'nvg'} && !$spe{'hum'}) {
            push @{$ogs{'Aala:cxgn_cxhn_cxsn_ccgn_cctn_aagy_aaty_hmgn_nvgn_humn'}}, $sf;
            #Medusozoa specific (w/ human)
        } elsif (($spe{'cxg'} || $spe{'cxh'}
        || $spe{'cxs'}) && ($spe{'ccg'} || $spe{'cct'}) && ($spe{'aag'} || $spe{'aat'}) && $spe{'hmg'} && !$spe{'nvg'} && $spe{'hum'}) {
            push @{$ogs{'Medusozoa_w/Hum:cxgy_cxhy_cxsy_ccgy_ccty_aagy_aaty_hmgy_nvgn_humy'}}, $sf;
            #Medusozoa specific (no human)
        } elsif (($spe{'cxg'} || $spe{'cxh'}
        || $spe{'cxs'}) && ($spe{'ccg'} || $spe{'cct'}) && ($spe{'aag'} || $spe{'aat'}) && $spe{'hmg'} && !$spe{'nvg'} && !$spe{'hum'}) {
            push @{$ogs{'Medusozoa:cxgy_cxhy_cxsy_ccgy_ccty_aagy_aaty_hmgy_nvgn_humn'}}, $sf;
            #Non-Medusozoa specific (w/ human)
        } elsif (!$spe{'cxg'} && !$spe{'cxh'}
        && !$spe{'cxs'} && !$spe{'ccg'} && !$spe{'cct'} && !$spe{'aag'} && !$spe{'aat'} && !$spe{'hmg'} && $spe{'nvg'} && $spe{'hum'}) {
            push @{$ogs{'non-Medusozoa:cxgn_cxhn_cxsn_ccgn_cctn_aagn_aatn_hmgn_nvgy_humy'}}, $sf;
       }
    }
    open OUT, ">outfile" or die "cannot open OUTFILE:$!";
    foreach my $key (sort keys %ogs) {
        my $count = scalar (@{$ogs{$key}});
        print "$key: $count\n";
        print OUT "$key ";
        foreach my $file (@{$ogs{$key}}) {
            print OUT "$file ";
        }
        print OUT "\n";
    }
}

sub get_sp {
    my $val = shift;
    if ($val eq 'Cassiopea') {
        return 'cxg';
    } elsif ($val eq 'Cxam.HS.Trinity.transdecoder') {
        return 'cxh';
    } elsif ($val eq 'Cxam.Scyph.transdecoder') {
        return 'cxs';
    } elsif ($val eq 'Ccrux.Trinity') {
        return 'cct';
    } elsif ($val eq 'Alatina') {
        return 'aag';
    } elsif ($val eq 'Calvadosia') {
        return 'ccg';
    } elsif ($val eq 'Aalata.transdecoder') {
        return 'aat';
    } elsif ($val eq 'hydra2.0') {
        return 'hmg';
    } elsif ($val eq 'Nvec') {
        return 'nvg';
    } elsif ($val eq 'HumRef2017.2') {
        return 'hum';
    } elsif ($val eq 'ToxProt'){
	return 'tox';
    } else {
        die "unexpected: $val";
    }
}



