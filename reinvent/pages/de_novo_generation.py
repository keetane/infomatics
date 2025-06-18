import streamlit as st
import os

# ここにサイドバー設定のコードを直接記述
st.sidebar.title("設定")

# セッションステートに 'current_directory' がなければ初期値を設定
# この行は不要な場合もありますが、念のため記述しておきます。
# メインアプリが最初に実行されるため、通常はすでに初期化されています。
if 'current_directory' not in st.session_state:
    st.session_state.current_directory = os.getcwd()

# サイドバーにtext_inputを設置
# keyを指定することで、ページ遷移しても値が保持されます
current_dir_input = st.sidebar.text_input(
    "カレントディレクトリを設定:",
    value=st.session_state.current_directory,
    key="sidebar_current_directory_input" # 全ページで共通のkeyを使用
)

# text_inputの値が変更されたらsession_stateを更新
if current_dir_input != st.session_state.current_directory:
    st.session_state.current_directory = current_dir_input
    st.rerun()


st.title("ページ1")

# session_stateからカレントディレクトリの値を取得
current_directory = st.session_state.get('current_directory', '設定されていません')

st.write(f"ページ1でのカレントディレクトリ: **{current_directory}**")

# ここでcurrent_directoryを使用して何らかの処理を行うことができます
if os.path.exists(current_directory) and os.path.isdir(current_directory):
    st.subheader(f"ページ1から見た '{current_directory}' の情報:")
    try:
        files_count = len(os.listdir(current_directory))
        st.info(f"含まれるファイル/ディレクトリの数: {files_count}")
    except Exception as e:
        st.error(f"ディレクトリ情報を取得できませんでした: {e}")
else:
    st.warning("設定されたディレクトリが見つからないか、有効なディレクトリではありません。")