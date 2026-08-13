import streamlit as st
import networkx as nx

from core.centrality import (
    centrality_dataframe,
    centrality_table_html,
    CENTRALITY_ORDER,
)

try:
    from core.centrality import (
        ccc_dataframe,
        ccc_table_html,
        CCC_ORDER,
    )
except ImportError:
    # Versi core/centrality.py belum diperbarui (atau masih cache).
    ccc_dataframe = None
    ccc_table_html = None
    CCC_ORDER = ["ccdc", "cccc", "ccbc"]


def render_centrality(G, graph_label: str):
    """Render Centrality section (same behavior as before).

    Returns the computed centrality dataframe.
    """
    st.subheader("Centrality")

    cent_key = f"centrality_df::{graph_label}"
    if cent_key not in st.session_state:
        use_gc_for_cent = False
        if G.number_of_nodes() > 0:
            try:
                use_gc_for_cent = (nx.number_connected_components(G.to_undirected(as_view=True)) > 1)
            except Exception:
                use_gc_for_cent = False

        st.session_state[cent_key] = centrality_dataframe(G, use_giant_component=use_gc_for_cent)

    df_cent = st.session_state[cent_key]

    c1, c2, c3 = st.columns([2, 2, 2], vertical_alignment="center")
    with c1:
        sort_by = st.selectbox(
            "Urutkan berdasarkan",
            CENTRALITY_ORDER,
            index=2,
            format_func=lambda k: k,
        )
    with c2:
        topk = st.slider("Top-k node", 5, min(50, G.number_of_nodes()), 10, 1)
    with c3:
        show_ranks = st.checkbox("Tampilkan rank", value=True)

    table_html = centrality_table_html(
        df_cent,
        sort_by=sort_by,
        topk=int(topk),
        digits=4,
        show_ranks=bool(show_ranks),
    )
    st.markdown(table_html, unsafe_allow_html=True)

    csv_data = df_cent.reset_index().to_csv(index=False)
    st.download_button(
        "Download centrality (CSV)",
        data=csv_data,
        file_name=f"centrality_{graph_label.replace(' ', '_')}.csv",
        mime="text/csv",
    )

    st.caption("Arah ↑: makin besar umumnya makin sentral. Hover pada nama metrik untuk melihat definisi singkat.")
    return df_cent


def render_ccc_overlap(G, res: dict, graph_label: str):
    """Render overlap-aware Community-Consideration Centrality (CCC).

    - Membutuhkan hasil deteksi komunitas (res['communities']).
    - Menghasilkan tabel + CSV download seperti centrality biasa.
    """
    st.subheader("Community-Consideration Centrality (Overlap)")
    if ccc_dataframe is None or ccc_table_html is None:
        st.error(
            "Modul CCC overlap belum tersedia di core/centrality.py. "
            "Silakan ganti file core/centrality.py dengan versi yang sudah memuat ccc_dataframe, "
            "lalu restart Streamlit (Ctrl+C, jalankan ulang)."
        )
        return None


    if not res or not isinstance(res, dict):
        st.info("Jalankan deteksi komunitas terlebih dahulu.")
        return None

    comms = (res or {}).get("communities", []) or []
    if len(comms) == 0:
        st.info("Belum ada komunitas (hasil kosong). Jalankan deteksi komunitas terlebih dahulu.")
        return None

    meta = (res or {}).get("meta", {}) or {}
    run_id = meta.get("run_id", "RUN")

    # default: pakai GC jika graph tidak terhubung (konsisten dengan centrality biasa)
    use_gc_for_ccc = False
    if G.number_of_nodes() > 0:
        try:
            use_gc_for_ccc = (nx.number_connected_components(G.to_undirected(as_view=True)) > 1)
        except Exception:
            use_gc_for_ccc = False

    with st.expander("Pengaturan CCC (overlap)", expanded=True):
        alpha = st.slider(
            "alpha (trade-off komunitas vs global)",
            min_value=0.0,
            max_value=1.0,
            value=0.2,
            step=0.05,
            help="alpha=0 fokus komunitas; alpha=1 kembali ke metrik global.",
        )

        c1, c2, c3 = st.columns([1, 1, 1], vertical_alignment="center")
        with c1:
            do_ccdc = st.checkbox("Hitung CCDC (degree)", value=True)
        with c2:
            do_cccc = st.checkbox("Hitung CCCC (closeness)", value=True)
        with c3:
            do_ccbc = st.checkbox("Hitung CCBC (betweenness)", value=False)

        bet_k = None
        seed = None
        if do_ccbc:
            st.caption("CCBC (betweenness) dihitung exact (tanpa aproksimasi/sampling).")
            bet_k = None
            seed = None

    # cache key unik per run + alpha + opsi
    alpha_key = f"{float(alpha):.2f}"
    opt_key = f"ccdc{int(do_ccdc)}_cccc{int(do_cccc)}_ccbc{int(do_ccbc)}"
    ccc_key = f"ccc_df::{graph_label}::{run_id}::a{alpha_key}::{opt_key}"

    if ccc_key not in st.session_state:
        st.session_state[ccc_key] = ccc_dataframe(
            G,
            comms,
            alpha=float(alpha),
            use_giant_component=use_gc_for_ccc,
            compute_ccdc=bool(do_ccdc),
            compute_cccc=bool(do_cccc),
            compute_ccbc=bool(do_ccbc),
            betweenness_k=bet_k,
            seed=seed,
        )

    df_ccc = st.session_state[ccc_key]

    # tabel
    metric_cols = [c for c in CCC_ORDER if c in df_ccc.columns and df_ccc[c].notna().any()]
    if len(metric_cols) == 0:
        st.info("Tidak ada metrik CCC yang dihitung (cek pilihan checkbox).")
        return df_ccc

    t1, t2, t3 = st.columns([2, 2, 2], vertical_alignment="center")
    with t1:
        sort_by = st.selectbox(
            "Urutkan berdasarkan (CCC)",
            metric_cols,
            index=min(len(metric_cols) - 1, 0),
        )
    with t2:
        topk = st.slider("Top-k node (CCC)", 5, min(50, G.number_of_nodes()), 10, 1)
    with t3:
        show_ranks = st.checkbox("Tampilkan rank (CCC)", value=True, key=f"show_ranks_ccc::{ccc_key}")

    table_html = ccc_table_html(
        df_ccc,
        sort_by=sort_by,
        topk=int(topk),
        digits=4,
        show_ranks=bool(show_ranks),
    )
    st.markdown(table_html, unsafe_allow_html=True)

    # download
    csv_data = df_ccc.reset_index().to_csv(index=False)
    st.download_button(
        "Download CCC overlap (CSV)",
        data=csv_data,
        file_name=f"ccc_overlap_{graph_label.replace(' ', '_')}_{run_id}_a{alpha_key}.csv",
        mime="text/csv",
    )

    st.caption("CCC overlap: skor node = rata-rata berbobot (default 1/|M(i)|) atas komunitas yang diikuti node tersebut.")
    return df_ccc