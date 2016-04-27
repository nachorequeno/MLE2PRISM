#-------------------------------------------------------------------------------
#
#   MLE2PRISM  Copyright (C) 2016  Jose-Ignacio Requeno
#
#   This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
#   This is free software, and you are welcome to redistribute it under certain
#   conditions; type `show c' for details.
#
#-------------------------------------------------------------------------------
# File :  MLE2PRISM.perl
# Last version :  v1.0 ( 02/Feb/2016 )
# Description :  Translation of the maximum likelihood equations of a 
#       phylogenetic tree into the notation of the model checking tool PRISM.
#-------------------------------------------------------------------------------
# Historical report :
#
#   DATE :  02/Feb/2016
#   VERSION :  v1.0
#   AUTHOR(s) :  Jose-Ignacio Requeno
#
#-------------------------------------------------------------------------------

#!/bin/perl -w
  use strict;
  use FileHandle;
  use Bio::TreeIO;
  use Bio::SeqIO;
  use Bio::SearchIO; 
  use Bio::DB::Fasta;

  # Auxiliar functions #
  sub long_dna{
        my ($infile, @rest) = @_;
        my $seqout = new Bio::SeqIO(-file => $infile,
                                -format => 'fasta');
        my $sequence = $seqout->next_seq(); 
        $seqout->close();
        return (length($sequence->seq()));
  } 

  sub MLE_EQUATIONS{
        my ($subtree, $infile, $outfile, $pos, $alpha, $max_length, @rest) = @_;
        my (@alphabet) =  @$alpha;
        my $root = $subtree;
        my @descendents = $root->each_Descendent;
        # Initialize the likelihood of constant values (leafs)
        if (scalar @descendents == 0){
            my $root_sequence = lc($infile->seq($root->id, $pos, $pos));
            for my $element ( @alphabet ) {
                my $temp = undef;
                print $outfile "const int L_".$root->id."_".$pos."_".$element." = ";
                if ($root_sequence ne $element){
                    $temp .= "0";
                } else {
                    $temp .= "1";                    
                }
                $temp .= ";\n";
                print $outfile $temp;
            }
        }else{
            for my $child ( @descendents ) {
                MLE_EQUATIONS($child, $infile, $outfile, $pos, $alpha, $max_length);
            }
            for my $root_index ( @alphabet ) {
                my $temp = undef;
                print $outfile "\"L_".$root->id."_".$pos."_".$root_index."\" : ";
                for my $child ( @descendents ) {
                    # Branch length
                    my $length = sprintf("%f", $child->branch_length);
                                        $length = ($length == 0) ? 1 : $length;
                    $max_length = ($length > $max_length) ? $length : $max_length;
                    # MLE equation
                    $temp.="(";
                    for my $child_index ( @alphabet ) {
                                                my @descendents_child = $child->each_Descendent;
                                                if (scalar @descendents_child == 0){
                                                        #Lyl
                                                        $temp .= "L_".$child->id."_".$pos."_".$child_index."*";
                                                } else {
                                                        #Lyl
                                                        $temp .= "\"L_".$child->id."_".$pos."_".$child_index."\"*";
                                                }
                        #Pij(Kxy)
                                                $temp .= "filter(max, P=? [(x = ".$root_index.") U [".$length.",".$length."] (x = ".$child_index.")], x = ".$root_index.")+";
                    }
                    chop($temp);
                    $temp.=")*";
                }
                chop($temp);
                $temp.=";\n";
                print $outfile $temp;
            }
        }
        $_[5]=$max_length;
            return 0;
  }
  
  # MAIN #
  @ARGV == 4 or die "usage: MLE2PRISM INPUTFILE1.NEW INPUTFILE2.FASTA OUTFILE.PRISM\n";
  my ($infile1, $infile2, $outfile) = @ARGV;

  my $treeio = new Bio::TreeIO(-format => 'newick',
                                        -file => $infile1);
  my $fh = new FileHandle;
  $fh->open(">$outfile") or die "Could not open file\n";

  my $long_seq = long_dna($infile2);

  my $tree = $treeio->next_tree;
  my $root = $tree->get_root_node;
  my $db = Bio::DB::Fasta->new($infile2);
  
  # Alphabet (nucleic bases)
  my @alphabet = ("a", "c", "g", "t");
  
  # Likelihood equations for each position of the sequence
  my $max_length = 0;
  for(my $i=0; $i<$long_seq; $i++){
        my $root_sequence = lc($db->seq($root->id, $i, $i));
        my $temp = "\"L_".$root->id."_".$i."\":";
        for my $element ( @alphabet ) {
                $temp.= "p".$element."*";
                $temp.= "\"L_".$root->id."_".$i."_".$element."\"+";
        }
        chop($temp);
        $temp .= ";\n";
        print $fh $temp;
        MLE_EQUATIONS($root, $db, $fh, $i, \@alphabet, $max_length);
  }

  $fh->close();
  $treeio->close();
