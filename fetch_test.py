from Bio import Entrez

# NCBI requires an email address for API identification
Entrez.email = "aarav.jindal@example.com"  # Replace with your email address

def fetch_sequence(accession_id):
    """
    Fetches raw sequence data in FASTA format from the NCBI Nucleotide database.
    """
    try:
        handle = Entrez.efetch(
            db="nucleotide", 
            id=accession_id, 
            rettype="fasta", 
            retmode="text"
        )
        sequence_data = handle.read()
        handle.close()
        return sequence_data
    except Exception as e:
        return f"API Error: {str(e)}"

if __name__ == "__main__":
    # Test fetch using a known Thylacine mitochondrial locus (Accession: NC_004387)
    test_id = "NC_004387"
    print(f"Fetching sequence for accession ID: {test_id}...\n")
    fasta_result = fetch_sequence(test_id)
    
    # Print the first 250 characters of the retrieved sequence
    print("--- Output FASTA Sequence Preview ---")
    print(fasta_result[:250])
