from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from tqdm import tqdm
import os

human_fasta = "human.fa"
mouse_fasta = "mouse.fa"
output_file = "blast_results.txt"

# Count the number of sequences in the human FASTA file for tqdm
num_sequences = sum(1 for _ in SeqIO.parse(human_fasta, "fasta"))

with open(output_file, "w") as out_handle:

    # Parse the human sequences with tqdm for a progress bar
    for record in tqdm(SeqIO.parse(human_fasta, "fasta"), total=num_sequences, desc="Running BLAST"):
        human_seq_id = record.id
        human_seq = record.seq

        # Write the single sequence to a temporary FASTA file
        temp_query_file = "temp_query.fa"
        with open(temp_query_file, "w") as temp_handle:
            temp_handle.write(f">{human_seq_id}\n{human_seq}\n")

        # Run BLAST for this single sequence against the mouse database
        blastp_cline = NcbiblastpCommandline(query=temp_query_file, db="mouse_db", outfmt=5, out="temp_blast.xml", max_target_seqs=1)
        stdout, stderr = blastp_cline()

        # Parse BLAST results for this single sequence
        with open("temp_blast.xml") as result_handle:
            blast_record = NCBIXML.read(result_handle)

            if blast_record.alignments:
                alignment = blast_record.alignments[0]
                hsp = alignment.hsps[0]

                # Write results to the output file
                out_handle.write(f"Human ID: {human_seq_id}\n")
                out_handle.write(f"Mouse ID: {alignment.hit_id}\n")
                out_handle.write(f"Alignment:\n{hsp.query}\n{hsp.match}\n{hsp.sbjct}\n")
                out_handle.write(f"E-value: {hsp.expect}\n")
                out_handle.write(f"Bit Score: {hsp.bits}\n")
                out_handle.write("\n" + "="*50 + "\n\n")

        # Clean up temporary files
        os.remove(temp_query_file)
        os.remove("temp_blast.xml")

print(f"Results saved to {output_file}")
