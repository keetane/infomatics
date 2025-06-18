import streamlit as st
import os

# ここにサイドバー設定のコードを直接記述
st.sidebar.title("設定")

# セッションステートに 'current_directory' がなければ初期値を設定
if 'current_directory' not in st.session_state:
    st.session_state.current_directory = os.getcwd()

# サイドバーにtext_inputを設置
current_dir_input = st.sidebar.text_input(
    "カレントディレクトリを設定:",
    value=st.session_state.current_directory,
    key="sidebar_current_directory_input" # 全ページで共通のkeyを使用
)

# text_inputの値が変更されたらsession_stateを更新
if current_dir_input != st.session_state.current_directory:
    st.session_state.current_directory = current_dir_input
    st.rerun()


st.title("ページ2")

# session_stateからカレントディレクトリの値を取得
current_directory = st.session_state.get('current_directory', '設定されていません')

st.write(f"ページ2でのカレントディレクトリ: **{current_directory}**")

# ここでcurrent_directoryを使用して何らかの処理を行うことができます
st.write("このページでは、設定されたカレントディレクトリに基づいて追加の処理を実行できます。")

# 例: 設定されたディレクトリが空かどうか
if os.path.exists(current_directory) and os.path.isdir(current_directory):
    st.subheader(f"'{current_directory}' の状態:")
    try:
        if not os.listdir(current_directory):
            st.info("設定されたディレクトリは現在空です。")
        else:
            st.success("設定されたディレクトリにはファイル/ディレクトリが含まれています。")
    except Exception as e:
        st.error(f"ディレクトリの状態を確認できませんでした: {e}")
else:
    st.warning("設定されたディレクトリが見つからないか、有効なディレクトリではありません。")