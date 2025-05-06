import os
import subprocess
import pandas as pd
from rdkit.Chem import PandasTools, Draw
from rdkit import Chem
import streamlit as st
from datetime import datetime
from rdkit.Chem import Recap
import pubchempy as pcp

# Streamlit app
st.title("REINVENTer 4 Drug Discovery")
st.text("- De Novo Molecular Design with AI by AZ -")
st.markdown("[Learn more about REINVENT4](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-024-00812-5)")
st.markdown("---")  # Add a horizontal rule (line)

# Set up directories relative to the script's location
script_dir = os.path.dirname(os.path.abspath(__file__))  # Get the script's directory
reinvent_dir = os.path.join(script_dir, "~/Documents/apps/REINVENT4")  # Adjusted REINVENT4 directory path
wd = os.path.join(reinvent_dir, "wd")  # Working directory path
os.makedirs(wd, exist_ok=True)  # Create directory if it doesn't exist
st.text(f"Working directory: {wd}")
# ...existing code...

# Update paths to be relative to the script's directory
toml_dir = os.path.join(wd, 'toml')
input = os.path.join(wd, 'input')
results_dir = os.path.join(wd, 'results')
sampling_log = os.path.join(results_dir, 'log')
model = os.path.join(wd, 'model')
model_log = os.path.join(model, 'log')

priors_dir = os.path.join(reinvent_dir, "priors")
reinvent = os.path.join(priors_dir, "reinvent.prior")
lib = os.path.join(priors_dir, "libinvent.prior")
link = os.path.join(priors_dir, "linkinvent.prior")
mol2mol_high = os.path.join(priors_dir, "mol2mol_high_similarity.prior")
mol2mol_med = os.path.join(priors_dir, "mol2mol_medium_similarity.prior")
mol2mol_mmp = os.path.join(priors_dir, "mol2mol_mmp.prior")
mol2mol_scaffold_generic = os.path.join(priors_dir, "mol2mol_scaffold_generic.prior")
mol2mol_scaffold = os.path.join(priors_dir, "mol2mol_scaffold.prior")
mol2mol_similarity = os.path.join(priors_dir, "mol2mol_similarity.prior")
pubchem = os.path.join(priors_dir, "pubchem_ecfp4_with_count_with_rank_reinvent4_dict_voc.prior")

# Sidebar for input parameters
st.sidebar.text("Current directory: ")
st.sidebar.text(os.getcwd())
st.sidebar.header("Sampling Parameters")
num_smiles = st.sidebar.number_input("Number of SMILES", min_value=1, value=155)
device = st.sidebar.selectbox("Device", ["mps", "cpu", "cuda"])

# Allow users to select a model file from priors_dir
model_files = {
    "Reinvent": reinvent,
    "LibInvent": lib,
    "LinkInvent": link,
    "Mol2Mol High Similarity": mol2mol_high,
    "Mol2Mol Medium Similarity": mol2mol_med,
    "Mol2Mol MMP": mol2mol_mmp,
    "Mol2Mol Scaffold Generic": mol2mol_scaffold_generic,
    "Mol2Mol Scaffold": mol2mol_scaffold,
    "Mol2Mol Similarity": mol2mol_similarity,
    "PubChem": pubchem,
}
model_file = st.sidebar.selectbox("Model File", options=model_files.keys(), format_func=lambda x: x)
selected_model_file = model_files[model_file]

unique_molecules = st.sidebar.checkbox("Unique Molecules", value=True)
randomize_smiles = st.sidebar.checkbox("Randomize SMILES", value=True)

st.sidebar.markdown('## only for mol2mol')
sample_stategy = st.sidebar.selectbox("Sampling Strategy", ["beamsearch", 'multinomial'])
temperature = st.sidebar.number_input("Temperature", min_value=0.0, max_value=1.0, value=1.0, step=0.1)

# Generate TOML file
def generate_toml(selected_model_file, num_smiles=155, device="mps", show=False, sample_stategy=None, temperature=1.0):
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = os.path.join(results_dir, selected_model_file+f"_sampling_{timestamp}.csv")
    toml_content = f"""
run_type = "sampling"
device = "{device}"
json_out_config = "_sampling.json"

[parameters]
model_file = "{selected_model_file}"
output_file = "{output_file}"
num_smiles = {num_smiles}
unique_molecules = true
randomize_smiles = true
"""
    toml_path = os.path.join(toml_dir, "sampling.toml")
    with open(toml_path, "w") as f:
        f.write(toml_content)
    return toml_path, output_file

# Run REINVENT4
def run_reinvent(toml_path):
    log_file = os.path.join(results_dir, "sampling.log")
    subprocess.call([f"{home_dir}/miniforge3/envs/r4/bin/reinvent", "-l", log_file, toml_path])

# Display molecules
def display_molecules(csv_file):
    df = pd.read_csv(csv_file)
    PandasTools.AddMoleculeColumnToFrame(df, smilesCol="SMILES", molCol="Mol")
    img = Draw.MolsToGridImage(df.sample(min(len(df), 30)).Mol.tolist(), molsPerRow=5, subImgSize=(300, 200))
    return img

# Main app logic
# de novo molecular sampling
st.header("De Novo Molecular Sampling")
if st.button("Run Sampling"):
    if not os.path.exists(model_file):
        st.error("Model file does not exist!")
    else:
        st.info("Generating TOML file...")
        toml_path, output_file = generate_toml(model_file, num_smiles, device)
        st.info("Running REINVENT4...")
        run_reinvent(toml_path)
        st.success("Sampling completed!")
        st.info("Displaying sampled molecules...")
        st.image(display_molecules(output_file))
st.markdown("---")  # Add a horizontal rule (line)

# Fetch SMILES from PubChem
st.header("Fetch SMILES from PubChem")
compound_name = st.text_input("Enter compound name", value="MMAE")
if st.button("Fetch SMILES"):
    try:
        smiles = pcp.get_compounds(compound_name, 'name')[0].isomeric_smiles
        st.success(f"SMILES for {compound_name}: {smiles}")
    except Exception as e:
        st.error(f"Error fetching SMILES for {compound_name}: {e}")
else:
    smiles = ""  # Default blank

# Display fetched SMILES
st.text_area("Fetched SMILES", value=smiles, height=68)

# Recap decomposition and molecule visualization
st.header("Recap Decomposition and Visualization")
if st.button("Run Recap Decomposition"):
    try:
        mol = Chem.MolFromSmiles(smiles)
        recap_tree = Recap.RecapDecompose(mol)
        modified_mols = []
        modified_smiles = []

        # Function to relabel dummy atoms for LibInvent
        def relabel_dummy_atoms(mol, map_num=1):
            rw_mol = Chem.RWMol(mol)
            for atom in rw_mol.GetAtoms():
                if atom.GetAtomicNum() == 0:  # Dummy atom (*)
                    atom.SetAtomMapNum(map_num)
            return rw_mol.GetMol()

        # Process Recap fragments
        for frag_smiles, node in recap_tree.children.items():
            frag_mol = node.mol
            largest_frag = max(frag_smiles.split('.'), key=lambda x: Chem.MolFromSmiles(x).GetNumAtoms())
            modified = relabel_dummy_atoms(Chem.MolFromSmiles(largest_frag), map_num=1)
            modified_mols.append(modified)
            modified_smiles.append(Chem.MolToSmiles(modified))

        # Add parent molecule
        modified_mols.insert(0, mol)
        modified_smiles.insert(0, smiles)

        # Display molecules
        st.image(Draw.MolsToGridImage(modified_mols, molsPerRow=3, subImgSize=(300, 300), legends=modified_smiles))
        st.success("Recap decomposition completed and molecules displayed.")
    except Exception as e:
        st.error(f"Error during Recap decomposition: {e}")