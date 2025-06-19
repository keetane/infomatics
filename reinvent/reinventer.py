import streamlit as st
import os
import datetime # datetimeのインポートを統一
from rdkit import Chem
from rdkit.Chem import Draw, PandasTools
import pandas as pd

st.set_page_config(
    page_title="REINVENTer 4 Drug Discovery",
    page_icon=":pill:",
    layout="wide",
    initial_sidebar_state="expanded",
)

# --- サイドバーの設定 ---
st.sidebar.header("Re:Inventer 4 Drug Discovery")
st.sidebar.markdown("[Learn more about REINVENT4](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-024-00812-5)")
st.sidebar.title("Working Directory")

# セッションステートに 'user_subfolder' がなければ初期値を 'all_users' に設定
if 'user_subfolder' not in st.session_state:
    st.session_state.user_subfolder = "all_users"

# サイドバーにサブフォルダ名を入力するテキストボックス
user_subfolder_input = st.sidebar.text_input(
    "~/Documents/apps/REINVENT4/wd/",
    value=st.session_state.user_subfolder
)

# 入力値が変更されたらセッションステートを更新し、再実行
if user_subfolder_input != st.session_state.user_subfolder:
    st.session_state.user_subfolder = user_subfolder_input
    st.rerun()

# --- ディレクトリパスの構築と作成 ---
# ユーザーのホームディレクトリを展開し、REINVENT4/wdまでのベースパスを構築
base_path = os.path.expanduser(os.path.join("~", "Documents", "apps", "REINVENT4", "wd"))

# 最終的な作業ディレクトリパスを決定
# もしユーザーがサブフォルダ名を入力していればそれを結合、そうでなければベースパスのみ
st.session_state.wd = os.path.join(base_path, st.session_state.user_subfolder)

# ディレクトリを作成 (既に存在すれば何もしない)
try:
    os.makedirs(st.session_state.wd, exist_ok=True)
except Exception as e:
    st.sidebar.error(f"ディレクトリの作成中にエラーが発生しました: {e}")

# --- CSVファイル選択のためのselectbox ---
st.sidebar.markdown("### Results CSV File Selection")

results_dir_to_search = os.path.join(st.session_state.wd, 'results')
csv_files = []
mols_path = None # Initialize mols_path

# Debugging: Show the directory being searched

if os.path.exists(results_dir_to_search) and os.path.isdir(results_dir_to_search):
    for root, _, files in os.walk(results_dir_to_search):
        for f in files:
            if f.endswith('.csv'):
                csv_files.append(os.path.join(root, f))
    
    csv_files.sort() # Sort files for better user experience

    if csv_files:
        # Create a list of relative paths for display in the selectbox
        display_files = [os.path.relpath(f, results_dir_to_search) for f in csv_files]
        display_files.insert(0, "--- ファイルを選択してください ---") # Add a default "Select a file" option
        
        selected_file = st.sidebar.selectbox(
            "Choose a CSV file:", 
            options=display_files,
            key="csv_selector" # Unique key for the widget
        )

        if selected_file != "--- ファイルを選択してください ---":
            # Get the absolute path of the selected file
            mols_path_index = display_files.index(selected_file) - 1 # Adjust index because of the inserted default option
            mols_path = csv_files[mols_path_index]
            # Debugging: Show the assigned mols_path
        else:
            st.sidebar.info("CSVファイルをサイドバーから選択してください。")
    else:
        st.sidebar.warning("resultsディレクトリ内にCSVファイルが見つかりません。")
else:
    st.sidebar.error("resultsディレクトリが見つからないか、アクセスできません。")

st.sidebar.text(mols_path)


# --- メインコンテンツの表示 ---
df = pd.read_csv(mols_path, sep=',', encoding='utf-8') if mols_path else pd.DataFrame()
df = df.sort_values(by='NLL', ascending=True)
df['ROMol'] = df['SMILES'].apply(lambda x : Chem.MolFromSmiles(x))
# st.dataframe(df) if mols_path else st.info("CSVファイルが選択されていません。")
img = Draw.MolsToGridImage(
    df['ROMol'],
    molsPerRow=3,
    subImgSize=(300, 200),
    legends=df['NLL'].astype(str).tolist() if 'NLL' in df.columns else None
)
# img = Draw.MolToImage(
#     df['ROMol'][0],
#     size=(300, 300),
#     kekulize=True,
#     wedgeBonds=True,
#     # fitImage=True
# )

st.text(
    type(df['NLL'][0])
)
st.image(img, caption="Molecules from the selected CSV file", use_column_width=True)
