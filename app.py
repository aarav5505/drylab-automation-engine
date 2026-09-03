import streamlit as st
from Bio import Entrez
import re

# --- PAGE CONFIG ---
st.set_page_config(page_title="VectaIQ | CRISPR Engine", page_icon="🧬")

# --- INITIALIZE SESSION STATE ---
if "user_authenticated" not in st.session_state:
    st.session_state.user_authenticated = False
if "user_email" not in st.session_state:
    st.session_state.user_email = ""

# --- BACKEND LOGIC ---
def resolve_target_and_find_best_grna(organism_name, gene_keyword, user_email):
    # Set Entrez email dynamically to the user's input email
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

# --- APP FLOW CONTROL ---

# STEP 1: MANDATORY EMAIL ENTRY
if not st.session_state.user_authenticated:
    st.title("🧬 Welcome to VectaIQ")
    st.write("Enter your email address to authenticate NCBI API requests before using the engine.")
    
    email_input = st.text_input("User Email Address", placeholder="researcher@university.edu")
    
    if st.button("Authenticate & Enter Engine"):
        if email_input and "@" in email_input and "." in email_input:
            st.session_state.user_authenticated = True
            st.session_state.user_email = email_input
            st.rerun()
        else:
            st.error("Please enter a valid email address.")

# STEP 2: MAIN ENGINE INTERFACE
else:
    st.sidebar.write(f"NCBI API User: **{st.session_state.user_email}**")
    if st.sidebar.button("Change Email / Reset"):
        st.session_state.user_authenticated = False
        st.session_state.user_email = ""
        st.rerun()

    st.title("🧬 VectaIQ: CRISPR-Cas9 gRNA Discovery Engine")
    st.write("Automated sequence retrieval, PAM site scanning, and candidate scoring.")

    col1, col2 = st.columns(2)
    with col1:
        organism = st.text_input("Organism Name", value="Homo sapiens")
    with col2:
        gene = st.text_input("Gene or Locus Keyword", value="BRCA1")

    if st.button("Run gRNA Discovery Engine"):
        with st.spinner("Searching NCBI and scoring targets..."):
            # Pass user_email to function
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
