from modeller import *
from modeller.scripts import complete_pdb
from Bio import SeqIO

def initialize_modeller_environment(atom_files_directory):
    """
    Initializes the Modeller environment and sets the directory for atomic files.
    """
    env = Environ()
    env.io.atom_files_directory = atom_files_directory
    return env

def complete_protein_model(env, pdb_file, output_file):
    """
    Reads in a protein model from a PDB file, potentially completes it, and saves it.
    """
    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')
    mdl = complete_pdb(env, pdb_file)
    mdl.write(output_file)

def basic_model_read_write(env, pdb_file, output_file):
    """
    Reads and writes a protein model with minimal processing.
    """
    mdl = Model(env, file=pdb_file)
    mdl.write(file=output_file)

def convert_sequence_format(input_file, input_format, output_file, output_format):
    """
    Converts sequence files from one format to another using Biopython.
    """
    records = SeqIO.parse(input_file, input_format)
    count = SeqIO.write(records, output_file, output_format)
    print(f"Converted {count} records")

def prepare_sequence_database_and_profile(env, seq_database_file, target_alignment_file):
    """
    Prepares the sequence database, reads target sequences, and builds profiles for homology search.
    """
    log.verbose()
    
    # Prepare the input files and read in the sequence database
    sdb = SequenceDB(env)
    sdb.read(seq_database_file=seq_database_file, seq_database_format='PIR', chains_list='ALL', minmax_db_seq_len=(30, 4000), clean_sequences=True)
    sdb.write(seq_database_file='pdb_95.bin', seq_database_format='BINARY', chains_list='ALL')
    sdb.read(seq_database_file='pdb_95.bin', seq_database_format='BINARY', chains_list='ALL')

    # Read in the target sequence/alignment and convert to profile
    aln = Alignment(env)
    aln.append(file=target_alignment_file, alignment_format='PIR', align_codes='ALL')
    prf = aln.to_profile()
    prf.build(sdb, matrix_offset=-450, rr_file='${LIB}/blosum62.sim.mat', gap_penalties_1d=(-500, -50), n_prof_iterations=1, check_profile=False, max_aln_evalue=0.01)
    prf.write(file='build_profile.prf', profile_format='TEXT')
    aln = prf.to_alignment()
    aln.write(file='build_profile.ali', alignment_format='PIR')

def align_and_compare_structures(env, pdb_chains):
    """
    Aligns multiple protein structures, performs structural comparison, and visualizes relationships.
    """
    aln = Alignment(env)
    for pdb, chain in pdb_chains:
        m = Model(env, file=pdb, model_segment=('FIRST:'+chain, 'LAST:'+chain))
        aln.append_model(m, atom_files=pdb, align_codes=pdb+chain)
    aln.malign()
    aln.malign3d()
    aln.compare_structures()
    aln.id_table(matrix_file='family.mat')
    env.dendrogram(matrix_file='family.mat', cluster_cut=-1.0)

if __name__ == "__main__":
    # Example usage of the functions defined above
    
    # Initialize environment
    env = initialize_modeller_environment(['.', r'C:\Users\Owner\Downloads'])
    
    # Complete protein model
    complete_protein_model(env, '1bpv.pdb', 'your_protein_display.pdb')
    
    # Basic model read and write
    basic_model_read_write(env, '1bpv.pdb', 'output.pdb')
    
    # Convert sequence format
    convert_sequence_format('sequence.fasta', 'fasta', 'sequence.clustal', 'clustal')
    
    # Prepare sequence database and profile
    prepare_sequence_database_and_profile(env, 'pdb_95.pir', 'TvLDH.ali')
    
    # Align and compare structures
    pdb_chains = [('1b8p', 'A'), ('1bdm', 'A'), ('1civ', 'A'), ('5mdh', 'A'), ('7mdh', 'A'), ('1smk', 'A')]
    align_and_compare_structures(env, pdb_chains)
