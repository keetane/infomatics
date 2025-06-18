import streamlit as st
import os

# メインページでもサイドバーをセットアップ
st.sidebar.title("working directory 設定")

# セッションステートに 'current_directory' がなければ初期値を 'all_users' に設定
if 'current_directory' not in st.session_state:
    st.session_state.current_directory = "all_users" # デフォルト値を変更

# サイドバーにtext_inputを設置
# keyを指定することで、ページ遷移しても値が保持されます
current_dir_input = st.sidebar.text_input(
    "~/Documents/infomatics/reinvent/wd/",
    value=st.session_state.current_directory,
    key="sidebar_current_directory_input" # 全ページで共通のkeyを使用
)

# text_inputの値が変更されたらsession_stateを更新
if current_dir_input != st.session_state.current_directory:
    st.session_state.current_directory = current_dir_input
    st.rerun()

st.title("メインページ")
st.write(f"現在のカレントディレクトリ: **{st.session_state.current_directory}**")

st.write("サイドバーからカレントディレクトリを設定してください。")

# 設定されたカレントディレクトリが存在するか確認し、存在すれば表示
if st.session_state.current_directory == "all_users":
    st.info("カレントディレクトリは 'all_users' に設定されています。")
elif os.path.exists(st.session_state.current_directory) and os.path.isdir(st.session_state.current_directory):
    st.subheader(f"'{st.session_state.current_directory}' の内容:")
    try:
        files = os.listdir(st.session_state.current_directory)
        if files:
            for file in files:
                st.markdown(f"- `{file}`")
        else:
            st.info("このディレクトリは空です。")
    except Exception as e:
        st.error(f"ディレクトリの内容を読み込めませんでした: {e}")
else:
    st.warning("設定されたディレクトリが見つからないか、有効なディレクトリではありません。")