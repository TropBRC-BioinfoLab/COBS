import streamlit as st

from core.metrics import build_membership_maps
from core.viz import (
    draw_overview_community_map,
    draw_subgraph_for_community,
    draw_community_interaction_graph,
)


def render_visualization(G, res):
    """Render Visualisasi section (same behavior as before)."""
    st.subheader("Visualisasi")

    viz_mode = st.selectbox(
        "Mode visualisasi",
        ["Overview (community map)", "Komunitas terpilih (subgraph)", "Interaksi antar komunitas (supernode)"],
        index=0,
    )

    use_1hop = False
    selected_cid = None

    comms = (res or {}).get("communities", []) or []
    if not comms:
        st.info("Belum ada komunitas untuk divisualisasikan. Jalankan deteksi komunitas terlebih dahulu.")
        return

    node_to_comms, primary_comm = build_membership_maps(G, comms)

    if viz_mode == "Komunitas terpilih (subgraph)":
        selected_cid = st.number_input(
            "Pilih community_id",
            min_value=1,
            max_value=max(1, len(comms)),
            value=1,
            step=1,
        )
        use_1hop = st.checkbox("Tambahkan 1-hop neighbors (konteks)", value=True)

    if viz_mode == "Overview (community map)":
        fig = draw_overview_community_map(
            G,
            primary_comm=primary_comm,
            node_to_comms=node_to_comms,
            seed=42,
        )
        st.pyplot(fig, clear_figure=True, use_container_width=True)

    elif viz_mode == "Komunitas terpilih (subgraph)":
        fig = draw_subgraph_for_community(
            G,
            comms=comms,
            selected_cid=int(selected_cid),
            node_to_comms=node_to_comms,
            add_1hop=use_1hop,
            seed=42,
        )
        st.pyplot(fig, clear_figure=True, use_container_width=True)

    else:
        fig = draw_community_interaction_graph(
            G,
            comms=comms,
            seed=42,
        )
        st.pyplot(fig, clear_figure=True, use_container_width=True)
