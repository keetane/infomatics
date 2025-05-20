import os
import subprocess
import pandas as pd
from rdkit.Chem import PandasTools, Draw
from rdkit import Chem
import streamlit as st
from datetime import datetime
from rdkit.Chem import Recap
import pubchempy as pcp
import datetime
# Get the current timestamp
now = datetime.datetime.now()
# Format the timestamp as a string
time = now.strftime("%Y%m%d_%H%M")

st.set_page_config(
    page_title="REINVENTer 4 Drug Discovery",
    page_icon=":pill:",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Streamlit app
st.title("REINVENTer 4 Drug Discovery")
st.text("- De Novo Molecular Design with AI by AZ -")
st.markdown("[Learn more about REINVENT4](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-024-00812-5)")
st.markdown("---")  # Add a horizontal rule (line)

# Set up directories relative to the script's location
script_dir = os.path.dirname(os.path.abspath(__file__))  # Get the script's directory
home_dir = os.path.expanduser("~")  # Get the user's home directory
reinvent_dir = os.path.join(home_dir, "Documents/apps/REINVENT4")  # Adjusted REINVENT4 directory path
wd = os.path.join(reinvent_dir, "wd")  # Working directory path
os.makedirs(wd, exist_ok=True)  # Create directory if it doesn't exist
input_dir = os.path.join(wd, 'input')  # Input directory path
os.makedirs(input_dir, exist_ok=True)  # Create input directory if it doesn't exist
model_dir = os.path.join(wd, 'model')  # Model directory path
os.makedirs(model_dir, exist_ok=True)  # Create model directory if it doesn't exist
results_dir = os.path.join(wd, 'results')  # Results directory path
os.makedirs(results_dir, exist_ok=True)  # Create results directory if it doesn't exist
toml_dir = os.path.join(wd, 'toml')  # TOML directory path
os.makedirs(toml_dir, exist_ok=True)  # Create TOML directory if it doesn't exist
sampling_log = os.path.join(results_dir, 'log')  # Sampling log directory path
os.makedirs(sampling_log, exist_ok=True)  # Create sampling log directory if it doesn't exist

# # Update paths to be relative to the script's directory
# toml_dir = os.path.join(wd, 'toml')
# input = os.path.join(wd, 'input')
# results_dir = os.path.join(wd, 'results')
# sampling_log = os.path.join(results_dir, 'log')
# model = os.path.join(wd, 'model')
# model_log = os.path.join(model, 'log')

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


# Main app logic
# # de novo molecular sampling
# st.header("De Novo Molecular Sampling")
# if st.button("Run Sampling"):
#     if not os.path.exists(model_file):
#         st.error("Model file does not exist!")
#     else:
#         st.info("Generating TOML file...")
#         toml_path, output_file = generate_toml(model_file, num_smiles, device)
#         st.info("Running REINVENT4...")
#         run_reinvent(toml_path)
#         st.success("Sampling completed!")
#         st.info("Displaying sampled molecules...")
#         st.image(display_molecules(output_file))
# st.markdown("---")  # Add a horizontal rule (line)

# Fetch SMILES from PubChem
st.header("Fetch SMILES from PubChem")
compound_name = st.text_input("Enter compound name", value="ruxolitinib")
if st.button("Fetch SMILES from PubChem into text_area"):
    try:
        smiles = pcp.get_compounds(compound_name, 'name')[0].isomeric_smiles
        # st.success(f"SMILES for {compound_name}: {smiles}")
    except Exception as e:
        st.error(f"Error fetching SMILES for {compound_name}: {e}")
else:
    smiles = "C1CCC(C1)[C@@H](CC#N)N2C=C(C=N2)C3=C4C=CNC4=NC=N3"  # Default blank

# Display fetched SMILES
smiles = st.text_area("enter SMILES", value=smiles, height=68)
st.image(Draw.MolToImage(Chem.MolFromSmiles(smiles), size=(300, 300)), caption="Parent Molecule")

# save SMILES of parent molecule from text_area
with open(f"{input_dir}/parent.smi", "w") as f:
    f.write(f"{smiles}\n")

# Recap decomposition and molecule visualization
st.header("Children by Recap Decomposition")

mol = Chem.MolFromSmiles(smiles)
recap_tree = Recap.RecapDecompose(mol)
warhead_mols = []
warhead_smiless = []
child_mols = []
child_smiless = []

# ダミー原子にatom map番号を付ける関数
def relabel_dummy_atoms(mol, map_num=1):
    """mol内のすべてのダミー原子(*)にatom map番号を付ける → [*:1]形式に"""
    rw_mol = Chem.RWMol(mol)
    for atom in rw_mol.GetAtoms():
        if atom.GetAtomicNum() == 0:  # dummy atom (*)
            atom.SetAtomMapNum(map_num)
    return rw_mol.GetMol()

# 分子を1stepだけ分解し、ラベルを変換
for smiles, node in recap_tree.children.items():
    mol = node.mol
    warhead_mols.append(mol)
    warhead_smiless.append(Chem.MolToSmiles(mol))
    modified = relabel_dummy_atoms(mol, map_num=1)
    child_mols.append(modified)
    child_smiless.append(Chem.MolToSmiles(modified))

BB_namelist = [f'BB{i+1}' for i in range(len(warhead_mols))]
st.image(Draw.MolsToGridImage(warhead_mols, molsPerRow=3, subImgSize=(600,300), legends=BB_namelist))



# Sidebar for input parameters
st.sidebar.text("Current directory: ")
st.sidebar.text(os.getcwd())
st.sidebar.text("Working directory: ")
st.sidebar.text(wd)
st.sidebar.header("Sampling Parameters")
num_mols = st.sidebar.number_input("Number of SMILES", min_value=1, value=155)
device = st.sidebar.selectbox("Device", ["mps", "cpu", "cuda"])

unique_molecules = st.sidebar.checkbox("Unique Molecules", value=True)
randomize_smiles = st.sidebar.checkbox("Randomize SMILES", value=True)
overwrite = st.sidebar.checkbox("Overwrite", value=True)

st.sidebar.markdown('## mol2mol sampling')
# Allow users to select a model file from priors_dir
model_files = {
    "Mol2Mol High Similarity": mol2mol_high,
    "Mol2Mol Medium Similarity": mol2mol_med,
    "Mol2Mol MMP": mol2mol_mmp,
    "Mol2Mol Scaffold Generic": mol2mol_scaffold_generic,
    "Mol2Mol Scaffold": mol2mol_scaffold,
    "Mol2Mol Similarity": mol2mol_similarity,
    "PubChem": pubchem,
    "Reinvent": reinvent,
}
model_file = st.sidebar.selectbox("Model File", options=model_files.keys(), format_func=lambda x: x)
selected_model_file = model_files[model_file]

# other options for mol2mol sampling
sample_stategy = st.sidebar.selectbox("Sampling Strategy", ["beamsearch", 'multinomial'])
temperature = st.sidebar.number_input("Temperature", min_value=0.0, max_value=1.0, value=1.0, step=0.1)

# Generate TOML file
def generate_toml(method=selected_model_file, num_smiles=num_mols, smiles_file=None, device=device, show=False, sample_strategy=None, temperature=1.0):
    if method == link:
        filename = 'LinkInvent'
    elif method == lib:
        filename = 'LibInvent'
    else:
        filename = method.split("/")[-1].replace(" ", "_")[:-6]
        
    if overwrite is not True:
        filename = method.split("/")[-1].replace(" ", "_") + '_' + time
    else:
        pass
    output_file = os.path.join(results_dir, f"{filename}.csv")
    toml_content = f"""
run_type = "sampling"
device = "{device}"
json_out_config = "{sampling_log}/_sampling.json"

[parameters]
model_file = "{method}"
output_file = "{output_file}"
num_smiles = {num_smiles}
"""
    if unique_molecules is True:
        toml_content += f'''
        unique_molecules = true
        '''
    else:
        pass
    if randomize_smiles is True:
        toml_content += f'''
        randomize_smiles = true
        '''
    else:
        pass

    if smiles_file is not None:
        toml_content += f'''
        smiles_file = "{smiles_file}"  # 1 compound per line
        '''
    if sample_strategy == 'beamsearch':
        toml += f'''
        sample_strategy = "beamsearch"  # multinomial or beamsearch (deterministic)
        '''
    elif sample_strategy == 'multinomial':
        toml += f'''
        sample_strategy = "multinomial"  # multinomial or beamsearch (deterministic)
        temperature = {temperature} # temperature in multinomial sampling
        '''
    else:
        pass

    toml_path = os.path.join(toml_dir, "sampling.toml")
    with open(toml_path, "w") as f:
        f.write(toml_content)
    return toml_path, output_file

# Run REINVENT4
def run_reinvent(toml_path):
    log_file = os.path.join(sampling_log, "sampling.log")
    subprocess.call([f"{home_dir}/miniforge3/envs/r4/bin/reinvent", "-l", log_file, toml_path])

# Display molecules
def display_molecules(csv_file):
    df = pd.read_csv(csv_file)
    PandasTools.AddMoleculeColumnToFrame(df, smilesCol="SMILES", molCol="Mol")
    img = Draw.MolsToGridImage(df.sample(min(len(df), 30)).Mol.tolist(), molsPerRow=5, subImgSize=(300, 200))
    return img

# mol2mol sampling
if st.sidebar.button('Mol2Mol Sampling'):
    # Generate TOML file for mol2mol sampling
    toml_path, output_file = generate_toml(method=selected_model_file, smiles_file=f'{input_dir}/parent.smi', num_smiles=num_mols, device=device)
    # Run REINVENT4
    run_reinvent(toml_path)
    st.success("Sampling completed!")

# LibInvent用の親分子を選択
st.sidebar.markdown('## LibInvent Sampling')
selected_child = st.sidebar.multiselect(
    "Select Child Molecule", options=BB_namelist, default=["BB2"], key="child_multiselect"
)
if selected_child:
    selected_child_index = BB_namelist.index(selected_child[0])
    selected_child_smiles = child_smiless[selected_child_index]
    with open(f"{input_dir}/child.smi", "w") as f:
        f.write(f"{selected_child_smiles}\n")
else:
    selected_child_smiles = ""
    st.sidebar.text("No parent selected.")

if st.sidebar.button('LibInvent Sampling'):
    # Generate TOML file for LibInvent
    toml_path, output_file = generate_toml(method=lib, smiles_file=f'{input_dir}/child.smi', num_smiles=num_mols, device=device)
    # Run REINVENT4
    run_reinvent(toml_path)
    st.success("Sampling completed!")

# LinkInvent用のwarheadを選択
st.sidebar.markdown('## LinkInvent')
selected_warhead = st.sidebar.multiselect(
    "Select Warhead", options=BB_namelist, default=["BB1", "BB3"], key="warhead_multiselect"
)
if selected_warhead:
    selected_warhead_indices = [BB_namelist.index(w) for w in selected_warhead]
    selected_warhead_smiles_list = [warhead_smiless[i] for i in selected_warhead_indices]
    # st.sidebar.text(f"Selected Warhead SMILES: {' | '.join(selected_warhead_smiles_list)}")
    with open(f"{input_dir}/warheads.smi", "w") as f:
        f.write('|'.join(selected_warhead_smiles_list))
else:
    selected_warhead_smiles = ""
    st.sidebar.text("No warhead selected.")

if st.sidebar.button('LinkInvnet'):
    # Generate TOML file for LinkInvent
    toml_path, output_file = generate_toml(method=link, smiles_file=f'{input_dir}/warheads.smi', num_smiles=num_mols, device=device)
    # Run REINVENT4
    run_reinvent(toml_path)
    st.success("Sampling completed!")
