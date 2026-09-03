import streamlit as st
from Bio import Entrez
import re

# --- PAGE CONFIG ---
st.set_page_config(page_title="Pamoros | Genomic Suite", page_icon="🧬")

# --- INITIALIZE SESSION STATE ---
if "user_authenticated" not in st.session_state:
    st.session_state.user_authenticated = False
if "user_email" not in st.session_state:
    st.session_state.user_email = ""

# ==========================================
# SECTION 1: CORE BACKEND ENGINES (LOCKED)
# ==========================================

def resolve_target_and_find_best_grna(organism_name, gene_keyword, user_email):
    Entrez.email = user_email
    search_query = f"{organism_name}[Organism] AND {gene_keyword}[Gene/Title]"
    try:
        handle = Entrez.esearch(db="nucleotide", term=search_query, retmax=1)
        search_results = Entrez.read(handle)
        handle.close()
    except Exception as e:
        return f"Error connecting to NCBI: {str(e)}"
    
    id_list = search_results.get("IdList", [])
    if not id_list:
        return f"Error: No sequence records found for '{organism_name}' with keyword '{gene_keyword}'."
    
    accession_id = id_list[0]
    handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text")
    lines = handle.read().strip().split("\n")
    handle.close()
    
    sequence = "".join(lines[1:])
    pattern = r'(?=([ATCGA-Z]{20}[ATCGA-Z]GG))'
    scored_candidates = []
    
    for match in re.finditer(pattern, sequence):
        full_site = match.group(1)
        grna = full_site[:20]
        pam = full_site[20:]
        gc_pct = ((grna.count('G') + grna.count('C')) / 20) * 100
        
        score = 100
        flags = []
        if gc_pct < 40 or gc_pct > 60:
            score -= 30
            flags.append(f"Suboptimal GC ({gc_pct:.1f}%)")
        if "TTTT" in grna:
            score -= 40
            flags.append("Contains TTTT termination signal")
            
        seed_region = grna[-8:]
        seed_gc = ((seed_region.count('G') + seed_region.count('C')) / 8) * 100
        if seed_gc < 37.5:
            score -= 15
            flags.append("Low seed region GC stability")
            
        scored_candidates.append({
            "accession": accession_id,
            "grna": grna,
            "pam": pam,
            "gc_content": gc_pct,
            "score": max(0, score),
            "flags": flags if flags else ["Optimal parameters"]
        })
        
    if not scored_candidates:
        return "No valid CRISPR-Cas9 target sites found in the retrieved locus."
        
    scored_candidates.sort(key=lambda x: x['score'], reverse=True)
    return scored_candidates[0]


def fetch_sequence(organism_name, gene_keyword, user_email):
    Entrez.email = user_email
    if gene_keyword.strip():
        search_query = f"{organism_name}[Organism] AND {gene_keyword}[Gene/Title]"
    else:
        search_query = f"{organism_name}[Organism] AND refseq[filter]"
        
    try:
        handle = Entrez.esearch(db="nucleotide", term=search_query, retmax=1)
        search_results = Entrez.read(handle)
        handle.close()
        
        id_list = search_results.get("IdList", [])
        if not id_list:
            return None, f"No sequences found for {organism_name}."
            
        accession_id = id_list[0]
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text")
        lines = handle.read().strip().split("\n")
        handle.close()
        
        sequence = "".join(lines[1:]).upper()
        return sequence, accession_id
    except Exception as e:
        return None, str(e)


def compare_dna_sequences(seq1, seq2):
    min_len = min(len(seq1), len(seq2))
    if min_len == 0:
        return 0.0, 0, 0
        
    matches = sum(1 for i in range(min_len) if seq1[i] == seq2[i])
    similarity_pct = (matches / min_len) * 100
    
    gc1 = ((seq1.count('G') + seq1.count('C')) / len(seq1)) * 100
    gc2 = ((seq2.count('G') + seq2.count('C')) / len(seq2)) * 100
    
    return similarity_pct, gc1, gc2

# -------------------------------------------------------------
# 👇 ADD NEW BACKEND FUNCTIONS BELOW THIS LINE (FOR FUTURE FEATURES)
# -------------------------------------------------------------




# ==========================================
# SECTION 2: AUTHENTICATION GATE
# ==========================================

if not st.session_state.user_authenticated:
    st.title("🧬 Welcome to Pamoros")
    st.write("Enter your email address to authenticate NCBI API requests before using the platform.")
    
    email_input = st.text_input("User Email Address", placeholder="researcher@university.edu")
    
    if st.button("Authenticate & Enter Platform"):
        if email_input and "@" in email_input and "." in email_input:
            st.session_state.user_authenticated = True
            st.session_state.user_email = email_input
            st.rerun()
        else:
            st.error("Please enter a valid email address.")


# ==========================================
# SECTION 3: FRONTEND INTERFACE & NAVIGATION
# ==========================================

else:
    st.sidebar.write(f"NCBI API User: **{st.session_state.user_email}**")
    
    # -------------------------------------------------------------
    # 👇 ADD NEW FEATURE NAMES TO THIS LIST FOR SIDEBAR SELECTION
    # -------------------------------------------------------------
    app_mode = st.sidebar.radio(
        "Select Feature:",
        [
            "CRISPR gRNA Discovery", 
            "Comparative DNA Alignment"
            # Add new feature tab titles here!
        ]
    )
    
    if st.sidebar.button("Change Email / Reset"):
        st.session_state.user_authenticated = False
        st.session_state.user_email = ""
        st.rerun()

    # --- FEATURE 1: CRISPR DISCOVERY ENGINE ---
    if app_mode == "CRISPR gRNA Discovery":
        st.title("🧬 CRISPR-Cas9 gRNA Discovery Engine")
        st.write("Automated sequence retrieval, PAM site scanning, and candidate scoring.")

        col1, col2 = st.columns(2)
        with col1:
            organism = st.text_input("Organism Name", value="Homo sapiens")
        with col2:
            gene = st.text_input("Gene or Locus Keyword", value="BRCA1")

        if st.button("Run gRNA Discovery Engine"):
            with st.spinner("Searching NCBI and scoring targets..."):
                result = resolve_target_and_find_best_grna(organism, gene, st.session_state.user_email)
                
                if isinstance(result, dict):
                    st.success("Target Successfully Identified!")
                    st.metric("Top Candidate Score", f"{result['score']}/100")
                    
                    st.subheader("Candidate Details")
                    st.code(f"gRNA Sequence (20bp): {result['grna']}\nPAM Locus:           {result['pam']}", language="text")
                    
                    st.write(f"**NCBI Accession:** `{result['accession']}`")
                    st.write(f"**GC Content:** `{result['gc_content']:.1f}%`")
                    st.write(f"**Biophysical Assessment:** {', '.join(result['flags'])}")
                else:
                    st.error(result)

    # --- FEATURE 2: COMPARATIVE DNA ALIGNMENT ---
    elif app_mode == "Comparative DNA Alignment":
        st.title("📊 Comparative DNA Alignment")
        st.write("Compare genetic sequence similarity across species. Leave the gene keyword blank to run a broad genome-level correlation.")

        col1, col2 = st.columns(2)
        with col1:
            org1 = st.text_input("Organism 1 Name", value="Homo sapiens")
        with col2:
            org2 = st.text_input("Organism 2 Name", value="Pan troglodytes")
            
        gene_target = st.text_input("Specific Gene Keyword (Optional)", placeholder="Leave blank for broad genome alignment (e.g., BRCA1)")

        if st.button("Run DNA Comparison"):
            with st.spinner("Fetching sequences from NCBI and calculating homology..."):
                seq1, acc1 = fetch_sequence(org1, gene_target, st.session_state.user_email)
                seq2, acc2 = fetch_sequence(org2, gene_target, st.session_state.user_email)
                
                if seq1 and seq2:
                    similarity, gc1, gc2 = compare_dna_sequences(seq1, seq2)
                    
                    st.success("Sequence Alignment Complete!")
                    st.metric("Base Pair Homology Match", f"{similarity:.2f}%")
                    
                    st.subheader("Analysis Metrics")
                    col_a, col_b = st.columns(2)
                    with col_a:
                        st.write(f"**{org1}**")
                        st.write(f"Accession: `{acc1}`")
                        st.write(f"GC Content: `{gc1:.1f}%`")
                    with col_b:
                        st.write(f"**{org2}**")
                        st.write(f"Accession: `{acc2}`")
                        st.write(f"GC Content: `{gc2:.1f}%`")
                    
                    st.progress(similarity / 100)
                else:
                    st.error(f"Error retrieving sequences: {acc1 if not seq1 else acc2}")

    # -------------------------------------------------------------
    # 👇 ADD NEW FRONTEND INTERFACES BELOW THIS LINE USING `elif app_mode == "..."`
    # -------------------------------------------------------------
