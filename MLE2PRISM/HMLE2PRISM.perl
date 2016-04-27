#    HMLE2PRISM: Heuristic Maximum Likelihood Equations for PRISM Model Checking Tool

#-------------------------------------------------------------------------------
#
#   HMLE2PRISM  Copyright (C) 2016  Jose-Ignacio Requeno
#
#   This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
#   This is free software, and you are welcome to redistribute it under certain
#   conditions; type `show c' for details.
#
#-------------------------------------------------------------------------------
# File :  HMLE2PRISM.perl
# Last version :  v1.0 ( 02/Feb/2016 )
# Description :  Translation of a variant of maximum likelihood equations of a 
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

  sub DEPTH{
        my ($depth, $subtree, $leaf, @rest) = @_;
        my $root = $subtree;
        my $length = (defined $root->branch_length) ? ($depth + $root->branch_length) : $depth;
        my @descendents = $root->each_Descendent;
        if ((scalar @descendents == 0) && !($root->id eq $leaf)){
                return 0;
        }elsif ((scalar @descendents == 0) && ($root->id eq $leaf)){
                return $length;
        }else{
                #(scalar @descendents > 0)
                my $maxdepth = 0;
                my $temp_depth = 0;
                for my $child ( @descendents ) {
                        $temp_depth = DEPTH($length, $child, $leaf);
                        $maxdepth = ( $temp_depth > $maxdepth) ? $temp_depth : $maxdepth;
                        }
                return $maxdepth;
        }
  }
  
  # MAIN #
  @ARGV == 4 or die "usage: HMLE2PRISM INPUTFILE1.NEW INPUTFILE2.FASTA OUTFILE.PRISM\n";
  my ($infile1, $infile2, $outfile) = @ARGV;

  my $treeio = new Bio::TreeIO(-format => 'newick',
                                        -file => $infile1);
  my $fh = new FileHandle;
  $fh->open(">$outfile") or die "Could not open file\n";

  my $long_seq = long_dna($infile2);

  my $tree = $treeio->next_tree;
  my $root = $tree->get_root_node;
  my $db = Bio::DB::Fasta->new($infile2);
  my @descendents = $root->get_all_Descendents;
  my %hash_depth = ();
  
  # Alphabet (nucleic bases)
  my @alphabet = ("a", "c", "g", "t");
  
  # Initialize the likelihood of constant values (leafs)
  for my $child ( @descendents ) {
        if ($child->is_Leaf){
                $hash_depth{ $child->id } = DEPTH(0, $root, $child->id);  
        }
  }

  # Likelihood equations for each position of the sequence
  for(my $i=0; $i<$long_seq; $i++){
        my $temp = "\"L_".$root->id."_".$i."\":";
        for my $element ( @alphabet ) {
                $temp.= "p".$element."*";
                $temp.= "\"L_".$root->id."_".$i."_".$element."\"+";
        }
        chop($temp);
        $temp .= ";\n";
        print $fh $temp;
        for my $element ( @alphabet ) {
                $temp = "\"L_".$root->id."_".$i."_".$element."\":";
                for my $child ( @descendents ) {
                        if ($child->is_Leaf){
                                my $child_sequence = lc($db->seq($child->id, $i, $i));
                                $temp .= "filter(max, P=? [F=(".$hash_depth{ $child->id }.") x = ".$child_sequence."], x = ".$element.")*";
                        }
                }
                chop($temp);
                $temp .= ";\n";
                print $fh $temp;
        }
  }

  $fh->close();
  $treeio->close();
