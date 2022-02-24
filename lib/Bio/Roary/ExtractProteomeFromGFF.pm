package Bio::Roary::ExtractProteomeFromGFF;

# ABSTRACT: Take in a GFF file and create protein sequences in FASTA format

=head1 SYNOPSIS

Take in GFF files and create protein sequences in FASTA format
   use Bio::Roary::ExtractProteomeFromGFF;

   my $obj = Bio::Roary::ExtractProteomeFromGFF->new(
       gff_file        => $fasta_file,
     );
   $obj->fasta_file();

=cut

use Moose;
use Bio::SeqIO;
use Cwd;
use Bio::Roary::Exceptions;
use FileHandle;
use File::Basename;
use File::Temp;
use File::Copy;
use File::Path qw"make_path";
use Bio::Tools::GFF;
with 'Bio::Roary::JobRunner::Role';
with 'Bio::Roary::BedFromGFFRole';

has 'gff_file' => (is => 'ro', isa => 'Str', required => 1);
has 'apply_unknowns_filter' => (is => 'rw', isa => 'Bool', default => 1);
has 'maximum_percentage_of_unknowns' => (is => 'ro', isa => 'Num', default  => 5);
has 'output_filename' => (is => 'ro', isa => 'Str', lazy => 1,
                          builder => '_build_output_filename');
has 'fasta_file' => (is => 'ro', isa => 'Str', lazy => 1, builder => '_build_fasta_file');
has '_working_directory' => (is => 'ro', isa => 'File::Temp::Dir',
                             default => sub { File::Temp->newdir(DIR => getcwd, CLEANUP => 1); });
has '_working_directory_name' => (is => 'ro', isa => 'Str', lazy => 1,
                                  builder => '_build__working_directory_name');
has 'translation_table' => (is => 'rw', isa => 'Int', default => 11);

sub _build_fasta_file {
    my ($self) = @_;
    my $output_directory = $self->output_directory;
    my $output_filename = $self->output_filename;
    my $output_path = qq"${output_directory}/${output_filename}";
    my $input_gff = $self->gff_file;
    print "TESTME: Making directory: $output_directory\n";
    make_path($output_directory);
    my ($input_name, $input_path, $input_suffix) = fileparse($self->gff_file, qr/\.[^.]*/);
    print "TESTME: WHAT IS THE GFF FILE: $input_gff\n";
    if ($input_suffix eq '.gff') {
        $self->_extract_nucleotide_regions;
        $self->_convert_nucleotide_to_protein;
    } elsif ($input_suffix eq '.gb') {
        print "TESTME: Writing from gb file!\n";
        $self->_write_fasta_files_from_gb(
            input => $self->gff_file, output => $output_path,
            output_nt => $self->_nucleotide_fasta_file_from_gff_filename);
    }
    $self->_cleanup_fasta;
    $self->_cleanup_intermediate_files;
    $self->_filter_fasta_sequences($output_path);
    return($output_path);
}

sub _build__working_directory_name {
    my ($self) = @_;
    return $self->_working_directory->dirname();
}

sub _build_output_filename {
    my ($self) = @_;
    my ($filename, $directories, $suffix) = fileparse($self->gff_file, qr/\.[^.]*/);
    return join('/', ($self->_working_directory_name, $filename . '.faa'));
}

sub _cleanup_intermediate_files {
    my ($self) = @_;
    unlink($self->_unfiltered_output_filename);
    unlink($self->_fastatranslate_filename);
}

sub _nucleotide_fasta_file_from_gff_filename {
    my ($self) = @_;
    return join('/', ($self->output_directory,
                      join('.', ($self->output_filename, 'intermediate.fa'))));
}

sub _extracted_nucleotide_fasta_file_from_bed_filename {
    my ($self) = @_;
    return join('/', ($self->output_directory,
                      join('.', ($self->output_filename,'intermediate.extracted.fa'))));
}

sub _unfiltered_output_filename {
    my $self = shift;
    return join( '/', ($self->output_directory,
                       join('.', ($self->output_filename, 'unfiltered.fa'))));
}

sub _create_nucleotide_fasta_file_from_gff {
    my ($self, %args) = @_;
    my $output_filename = $self->_nucleotide_fasta_file_from_gff_filename;
    my $input_filename = $self->gff_file;
    my $output_fh = FileHandle->new(">${output_filename}");
    my $input_fh = FileHandle->new("<$input_filename");
    my $at_sequence = 0;
    my $lines_written = 0;
    while (my $line = <$input_fh>) {
        if ($line =~/^>/) {
            $at_sequence = 1;
        }
        if ($at_sequence == 1) {
            $lines_written++;
            print $output_fh $line;
        }
    }
    $output_fh->close();
    $input_fh->close();
    return($lines_written);
}

sub _write_fasta_files_from_gb {
    my ($self, %args) = @_;
    print "TESTME: INPUT: $args{input} OUTPUTNT: $args{output_nt} OUTPUTAA: $args{output}\n";
    my $input_fh = FileHandle->new("<$args{input}");
    my $output_nt_fh = FileHandle->new(">$args{output_nt}");
    my $output_aa_fh = FileHandle->new(">$args{output}");
    my $seqio = Bio::SeqIO->new(-format => 'genbank', -fh => $input_fh);
    my $output_nt = Bio::SeqIO->new(-format => 'Fasta', -fh => $output_nt_fh);
    my $output_aa = Bio::SeqIO->new(-format => 'Fasta', -fh => $output_aa_fh);
  SEQUENCES: while (my $seq = $seqio->next_seq) {
      my $seqid = $seq->id;
      my $contig_name = $seq->display_name;
      ## Write out the chromosome/contig
      my $contig_sequence = $seq->seq;
      my $nt_seq_obj = Bio::Seq->new(-display_id => $seqid, -seq => $contig_sequence);
     $output_nt->write_seq($nt_seq_obj);

      ## Now extract all the coding sequences, translate them, and write them out.
      my @feature_list = $seq->get_SeqFeatures();
    FEATURES: for my $feat (@feature_list) {
        my $type = $feat->primary_tag();
        next FEATURES unless ($type eq 'CDS');
        my $annot = $feat->annotation();
        my $name = $feat->display_name();
        my @ids = $feat->get_tag_values('protein_id');
        my $name = $ids[0];
        my $start = $feat->start();
        my $end = $feat->end();
        my $strand = $feat->strand();
        my $cds = $feat->seq->seq;
        my $aa = $feat->seq->translate->seq;
        my $aa_seq_obj = Bio::Seq->new(-display_id => $name, -seq => $aa);
        $output_aa->write_seq($aa_seq_obj);
    } ## End iterating over features of this contig/chromosome
  } ## End iterating over every chromosome/contig
    $input_fh->close();
    $output_nt_fh->close();
    $output_aa_fh->close();
}


sub _extract_nucleotide_regions {
    my ($self) = @_;

    $self->_create_nucleotide_fasta_file_from_gff;
    $self->_create_bed_file_from_gff;

    my $cmd =
        'bedtools getfasta -s -fi '
        . $self->_nucleotide_fasta_file_from_gff_filename
        . ' -bed '
        . $self->_bed_output_filename . ' -fo '
        . $self->_extracted_nucleotide_fasta_file_from_bed_filename
        . ' -name > /dev/null 2>&1';

    $self->logger->debug($cmd);
    system($cmd);
    unlink($self->_nucleotide_fasta_file_from_gff_filename);
    unlink($self->_bed_output_filename);
    unlink($self->_nucleotide_fasta_file_from_gff_filename . '.fai');
}

sub _cleanup_fasta {
    my $self = shift;
    my $infile = $self->_unfiltered_output_filename;
    my $outfile = join('/' ,($self->output_directory,$self->output_filename));
    return unless (-e $infile);

    open(my $in,  '<', $infile);
    open(my $out, '>', $outfile);
    while (my $line = <$in>) {
        chomp $line;
        if ($line =~ /^>/) {
            $line =~ s/"//g;
            # newer versions of Bedtools add (-) or (+) to the end of the sequence name, remove them
            $line =~ s!\([-+]\)!!;
        }

        if ($line =~ /^(>[^:]+)/) {
            $line = $1;
        }
        print $out "$line\n";
    }
    close $in;
    close $out;
}

sub _fastatranslate_filename {
    my ($self) = @_;
    return join('/', ($self->output_directory,
                      join( '.', ($self->output_filename, 'intermediate.translate.fa'))));
}

sub _fastatranslate {
    my ($self, $inputfile, $outputfile) = @_;

    my $input_fasta_file_obj = Bio::SeqIO->new(-file => $inputfile, -format => 'Fasta');
    my $output_protein_file_obj = Bio::SeqIO->new(-file => ">" . $outputfile,
                                                  -format => 'Fasta', -alphabet => 'protein');

    my %protein_sequence_objs;
    while (my $seq = $input_fasta_file_obj->next_seq) {
        $seq->desc(undef);
        my $protseq = $seq->translate(-codontable_id => $self->translation_table);
        $output_protein_file_obj->write_seq($protseq);
    }
    return 1;
}

sub _convert_nucleotide_to_protein {
    my ($self) = @_;
    $self->_fastatranslate($self->_extracted_nucleotide_fasta_file_from_bed_filename,
                           $self->_unfiltered_output_filename);
    unlink($self->_extracted_nucleotide_fasta_file_from_bed_filename);
}

sub _does_sequence_contain_too_many_unknowns {
    my ($self, $sequence_obj) = @_;
    my $maximum_number_of_Xs = int(($sequence_obj->length() * $self->maximum_percentage_of_unknowns) / 100);
    my $number_of_Xs_found = () = $sequence_obj->seq() =~ /X/g;
    if ($number_of_Xs_found > $maximum_number_of_Xs) {
        return 1;
    } else {
        return 0;
    }
}

sub _filter_fasta_sequences {
    my ($self, $filename) = @_;
    my $temp_output_file = $filename . '.tmp.filtered.fa';
    my $out_fasta_obj = Bio::SeqIO->new(-file => ">" . $temp_output_file, -format => 'Fasta');
    my $fasta_obj = Bio::SeqIO->new(-file => $filename, -format => 'Fasta');
    my $sequence_found = 0;
    while (my $seq = $fasta_obj->next_seq()) {
        if ($self->_does_sequence_contain_too_many_unknowns($seq)) {
            next;
        }
        $seq->desc(undef);
        $out_fasta_obj->write_seq($seq);
        $sequence_found = 1;
    }

    if ($sequence_found == 0) {
        $self->logger->error("Could not extract any protein sequences from "
                             . $self->gff_file
                             . ". Does the gff file contain the assembly as well as the annotation?" );
    }

    # Replace the original file.
    move($temp_output_file, $filename);
    return 1;
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;
