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
st.sidebar.header("Re:Inventer 4 Drug Discovery")
st.sidebar.markdown("[Learn more about REINVENT4](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-024-00812-5)")


# Set up directories relative to the script's location
wd_input = st.sidebar.text_input('set your working directory', value='~/Documents/apps/REINVENT4/wd')
wd = os.path.abspath(os.path.expanduser(wd_input))
os.makedirs(wd, exist_ok=True)  # Create directory if it doesn't exist
# script_dir = os.path.dirname(os.path.abspath(__file__))  # Get the script's directory
home_dir = os.path.expanduser("~")  # Get the user's home directory
reinvent_dir = os.path.join(home_dir, "Documents/apps/REINVENT4")  # Adjusted REINVENT4 directory path
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

# Set up file paths for priors
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


# Fetch SMILES from PubChem
st.header("Parent Molecule Drawer")

# separate with columns
col1, col2 = st.columns([2, 1])
with col1:

    compound_name = st.text_input("Enter compound name", value="ruxolitinib")
    if st.button("Fetch SMILES from PubChem into text_area"):
        try:
            smiles = pcp.get_compounds(compound_name, 'name')[0].isomeric_smiles
        except Exception as e:
            st.error(f"Error fetching SMILES for {compound_name}: {e}")
    else:
        smiles = "C1CCC(C1)[C@@H](CC#N)N2C=C(C=N2)C3=C4C=CNC4=NC=N3"  # Default blank

    # Display fetched SMILES
    smiles = st.text_area("enter SMILES", value=smiles, height=68)

with col2:
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
st.sidebar.header("Sampling Parameters")
num_mols = st.sidebar.number_input("Number of SMILES", min_value=1, value=155)
device = st.sidebar.selectbox("Device", ["mps", "cpu", "cuda"])

unique_molecules = st.sidebar.checkbox("Unique Molecules", value=True)
randomize_smiles = st.sidebar.checkbox("Randomize SMILES", value=True)
overwrite = st.sidebar.checkbox("Overwrite", value=True)

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
model_file = [] #st.sidebar.selectbox("Model File", options=model_files.keys(), format_func=lambda x: x)
selected_model_file = model_files['Reinvent']

st.sidebar.header('mol2mol options')
# other options for mol2mol sampling
sample_stategies = st.sidebar.selectbox("Sampling Strategy", ["beamsearch", 'multinomial'])
temperatures = st.sidebar.number_input("Temperature", min_value=0.0, max_value=1.0, value=1.0, step=0.1)



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
    if randomize_smiles is True:
        toml_content += f'''
        randomize_smiles = true
        '''
    if smiles_file is not None:
        toml_content += f'''
        smiles_file = "{smiles_file}"  # 1 compound per line
        '''
    # Fix: use toml_content instead of undefined toml
    if sample_strategy == 'beamsearch':
        toml_content += f'''
        sample_strategy = "beamsearch"  # multinomial or beamsearch (deterministic)
        '''
    elif sample_strategy == 'multinomial':
        toml_content += f'''
        sample_strategy = "multinomial"  # multinomial or beamsearch (deterministic)
        temperature = {temperature} # temperature in multinomial sampling
        '''

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

# 3columns for sampling
st.header("ReInvent Sampling")
col1, col2, col3 = st.columns([1, 1, 1])

# mol2mol sampling
with col1:
    st.markdown('#### Mol2Mol Sampling')
    model_file = st.selectbox("Model File", options=model_files.keys(), format_func=lambda x: x)
    selected_model_file = model_files[model_file]
    if st.button('Mol2Mol Sampling from Parent'):
        # Generate TOML file for mol2mol sampling
        toml_path, output_file = generate_toml(method=selected_model_file, smiles_file=f'{input_dir}/parent.smi', num_smiles=num_mols, device=device, sample_strategy=sample_stategies, temperature=temperatures)
        # Run REINVENT4
        run_reinvent(toml_path)
        st.success("Sampling completed!")

# LibInvent用の親分子を選択
with col2:
    st.markdown('#### LibInvent Sampling')
    selected_child = st.multiselect(
        "Select Child Molecule", options=BB_namelist, default=["BB2"], key="child_multiselect"
    )
    if selected_child:
        selected_child_index = BB_namelist.index(selected_child[0])
        selected_child_smiles = child_smiless[selected_child_index]
        with open(f"{input_dir}/child.smi", "w") as f:
            f.write(f"{selected_child_smiles}\n")
    else:
        selected_child_smiles = ""
        st.text("No parent selected.")

    if st.button('LibInvent Sampling'):
        # Generate TOML file for LibInvent
        toml_path, output_file = generate_toml(method=lib, smiles_file=f'{input_dir}/child.smi', num_smiles=num_mols, device=device)
        # Run REINVENT4
        run_reinvent(toml_path)
        st.success("Sampling completed!")

# LinkInvent用のwarheadを選択
with col3:
    st.markdown('#### LinkInvent')
    selected_warhead = st.multiselect(
        "Select Warhead", options=BB_namelist, default=["BB1", "BB2"], key="warhead_multiselect"
    )
    if selected_warhead:
        selected_warhead_indices = [BB_namelist.index(w) for w in selected_warhead]
        selected_warhead_smiles_list = [warhead_smiless[i] for i in selected_warhead_indices]
        # st.sidebar.text(f"Selected Warhead SMILES: {' | '.join(selected_warhead_smiles_list)}")
        with open(f"{input_dir}/warheads.smi", "w") as f:
            f.write('|'.join(selected_warhead_smiles_list))
    else:
        selected_warhead_smiles = ""
        st.text("No warhead selected.")

    if st.button('LinkInvnet'):
        # Generate TOML file for LinkInvent
        toml_path, output_file = generate_toml(method=link, smiles_file=f'{input_dir}/warheads.smi', num_smiles=num_mols, device=device)
        # Run REINVENT4
        run_reinvent(toml_path)
        st.success("Sampling completed!")

st.divider()


# Transfer learning from ChEMBL or CSV
from chembl_webresource_client.new_client import new_client
import requests
from rdkit.Chem import Descriptors

# parameters for TL in sidebar
st.sidebar.header("Transfer Learning Parameters")
# Epochs
num_epochs = st.sidebar.slider("Number of Epochs", min_value=1, max_value=100, value=3)
save_every_n_epochs = st.sidebar.slider("Save Every N Epochs", min_value=1, max_value=num_epochs, value=3)

# Batch sizes
batch_size = st.sidebar.number_input("Batch Size", min_value=1, max_value=1024, value=50)
sample_batch_size = st.sidebar.number_input("Sample Batch Size", min_value=1, max_value=2048, value=100)

# Number of reference molecules
num_refs = st.sidebar.number_input("Number of Reference Molecules", min_value=0, max_value=2000, value=100)

# Similarity settings
similarity_types = ["tanimoto", "cosine", "dice", "euclidean"]
pairs_type = st.sidebar.selectbox("Similarity Type", options=similarity_types, index=0)

pairs_upper_threshold = st.sidebar.slider("Upper Similarity Threshold", min_value=0.0, max_value=1.0, value=1.0)
pairs_lower_threshold = st.sidebar.slider("Lower Similarity Threshold", min_value=0.0, max_value=1.0, value=0.7)

pairs_min_cardinality = st.sidebar.number_input("Min Cardinality", min_value=1, max_value=1000, value=1)
pairs_max_cardinality = st.sidebar.number_input("Max Cardinality", min_value=pairs_min_cardinality, max_value=2000, value=199)

# 構成データを辞書形式で定義
config = {
    "num_epochs": num_epochs,
    "save_every_n_epochs": save_every_n_epochs,
    "batch_size": batch_size,
    "sample_batch_size": sample_batch_size,
    "num_refs": num_refs,
    "pairs": {
        "type": pairs_type,
        "upper_threshold": pairs_upper_threshold,
        "lower_threshold": pairs_lower_threshold,
        "min_cardinality": pairs_min_cardinality,
        "max_cardinality": pairs_max_cardinality
    }
}



st.header("Transfer learning from ChEMBL or CSV")

# generate toml file for transfer learning
def TL(TL_input):
    TL_input_filename = TL_input.split("/")[-1].replace(".smi", "")
    TL_toml=f'''
    run_type = "transfer_learning"
    device = "{device}"  # set torch device e.g. "cpu". For macOS, use "mps"
    tb_logdir = "{sampling_log}/tb_TL"  # name of the TensorBoard logging directory
    json_out_config = "{sampling_log}/json_transfer_learning.json"  # write this TOML to JSON

    [parameters]
    num_epochs = {num_epochs}  # number of steps to run
    save_every_n_epochs = {save_every_n_epochs}  # save checkpoint model file very N steps
    batch_size = {batch_size}  # batch size for training
    num_refs = {num_refs}  # number of reference molecules randomly chosen for similarity
                    # set this to zero for large datasets (>200 molecules)!
    sample_batch_size = {sample_batch_size}  # number of sampled molecules to compute sample loss
    # Uncomment one of the comment blocks below.  Each generator needs a model
    # file and possibly a SMILES file with seed structures.

    input_model_file = "{reinvent}" 
    smiles_file = "{TL_input}"  # read 1st column
    output_model_file = '{model_dir}/TL_{TL_input_filename}.model'  # sampled SMILES and NLL in CSV format
    validation_smiles_file = "{TL_input}"  # read 1st column
    
    # Define the type of similarity and its parameters
    pairs.type = '{pairs_type}'  # e.g. "tanimoto", "cosine", "dice", "euclidean"
    pairs.upper_threshold = {pairs_upper_threshold}  # upper threshold for similarity
    pairs.lower_threshold = {pairs_lower_threshold}  # lower threshold for similarity
    pairs.min_cardinality = {pairs_min_cardinality}  # minimum number of similar molecules
    pairs.max_cardinality = {pairs_max_cardinality}  # maximum number of similar molecules

    '''

    toml_path = f"{toml_dir}/TL.toml"

    with open(toml_path, "w") as f:
        f.write(TL_toml)
    return toml_path

col1, col2 = st.columns([3, 1])
with col1:
    st.markdown('#### from ChEMBL')
    uniprot_id = st.text_input("Enter Uniprot ID", value="P00533") #EGFR

    # Always fetch ChEMBL targets for the current Uniprot ID (not only on button press)
    try:
        target_records = list(new_client.target.get(
            target_components__accession=uniprot_id
        ).only(
            'target_chembl_id', 
            'organism', 
            'pref_name', 
            'target_type'
        ))
        df = pd.DataFrame.from_records(target_records)
    except Exception as e:
        st.error(f"Error fetching ChEMBL data for {uniprot_id}: {e}")
        df = pd.DataFrame()

    if not df.empty and 'target_chembl_id' in df.columns:
        st.dataframe(df, use_container_width=True)
        selected_target = st.selectbox(
            "Select ChEMBL target", 
            options=df['target_chembl_id'].tolist()
        )
    else:
        st.warning("No valid ChEMBL targets found for the given Uniprot ID.")
        selected_target = None

    if st.button('Transfer learning from ChEMBL'):
        if selected_target is not None:
            # 1. 活性のある化合物のデータ取得
            activities = new_client.activity.filter(
                target_chembl_id__in=selected_target,
                pchembl_value__gte=5,
                assay_type='B',
            ).only([
                'molecule_chembl_id',
                'molecule_pref_name',
                'target_pref_name',
                'parent_molecule_chembl_id',
                'pchembl_value',
                'canonical_smiles',
                'document_chembl_id',
                'document_journal',
            ])
            df = pd.DataFrame.from_records(activities)
            # st.dataframe(df, use_container_width=True)
            # 2. 主成分（分子量最大の構成要素）を抽出する関数
            def extract_main_component(smiles):
                try:
                    fragments = smiles.split('.')
                    max_mol = None
                    max_weight = -1
                    for frag in fragments:
                        mol = Chem.MolFromSmiles(frag)
                        if mol is None:
                            continue
                        mw = Descriptors.MolWt(mol)
                        if mw > max_weight:
                            max_mol = mol
                            max_weight = mw
                    if max_mol:
                        return Chem.MolToSmiles(max_mol, isomericSmiles=False)
                    else:
                        return None
                except Exception as e:
                    print(f"Error splitting SMILES: {smiles} -> {e}")
                    return None

            # 3. 分子量最大部分のSMILES抽出 → 標準化処理
            df["smiles"] = df["canonical_smiles"].apply(extract_main_component)

            # 4. 無効なエントリを削除
            df_clean = df.dropna(subset=["smiles"])

            # 5. 保存
            df_clean[['smiles']].to_csv(input_dir + f'/{selected_target}.smi', sep='\t', index=False)
            df_clean.to_csv(input_dir + f'/{selected_target}.csv', sep='\t', index=False)

            st.text(f"Saved {len(df_clean)} ligands with largest fragment SMILES to input/{selected_target}.csv and input/{selected_target}.smi")

            # st.dataframe(df_clean, use_container_width=True)

            # 6. TOMLファイルを生成
            toml_path = TL(input_dir + f'/{selected_target}.smi')

            # 7. REINVENT4を実行
            run_reinvent(toml_path)
            st.success("Transfer learning from ChEMBL completed!")
with col2:
    st.markdown('#### from csv')
    uploaded_file = st.file_uploader("Upload CSV file", type=["csv"])
    if uploaded_file is not None:
        try:
            df = pd.read_csv(uploaded_file, sep='\t')
            if 'smiles' in df.columns:
                # overwrite_TL = st.checkbox("overwrite TL model", value=True)
                if st.button('Transfer learning from csv'):
                    # Save uploaded file to input_dir and use its path for TL
                    csv_save_path = os.path.join(input_dir, uploaded_file.name)
                    with open(csv_save_path, "wb") as f:
                        f.write(uploaded_file.getbuffer())
                    # Generate .smi file from the 'smiles' column
                    smi_save_path = csv_save_path.rsplit('.', 1)[0] + ".smi"
                    df[['smiles']].to_csv(smi_save_path, sep='\t', index=False)
                    toml_path = TL(smi_save_path)
                    # Run REINVENT4
                    run_reinvent(toml_path)
                    st.success("Transfer learning from CSV completed!")
            else:
                st.error("CSV file must contain a 'smiles' column.")
        except Exception as e:
            st.error(f"Error reading CSV file: {e}")