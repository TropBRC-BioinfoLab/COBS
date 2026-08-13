import math
import streamlit as st
import streamlit.components.v1 as components


def render_communities_table(rows, overlap_set):
    css = """
    <style>
      .com-table-wrap {width:100%;}
      table.com-table {width:100%; border-collapse:collapse; table-layout:fixed;}
      table.com-table th, table.com-table td{
        border:1px solid rgba(49,51,63,0.15);
        padding:10px;
        vertical-align:top;
        font-size:14px;
      }
      table.com-table th{
        background:rgba(49,51,63,0.05);
        font-weight:600;
        text-align:left;
      }

      /* Kolom 1 & 2 dibuat sempit */
      table.com-table col.c1{width:120px;}
      table.com-table col.c2{width:90px;}
      table.com-table col.c3{width:auto;}

      /* Badge nodes + wrapping ke baris berikutnya */
      .node-wrap{display:flex; flex-wrap:wrap; gap:6px; align-items:flex-start;}
      .node-badge{
        display:inline-block;
        padding:2px 10px;
        border-radius:999px;
        background:rgba(0,0,0,0.06);
        line-height:1.6;
        white-space:nowrap;
      }
      .node-badge.overlap{
        background: rgba(255,193,7,0.55);
        font-weight:700;
      }
    </style>
    """

    def badge(n):
        cls = "node-badge overlap" if n in overlap_set else "node-badge"
        return f"<span class='{cls}'>{n}</span>"

    body = []
    for r in rows:
        nodes_html = "<div class='node-wrap'>" + "".join(badge(n) for n in r["nodes"]) + "</div>"
        body.append(
            f"<tr><td>{r['community_id']}</td><td>{r['size']}</td><td>{nodes_html}</td></tr>"
        )

    html = f"""
    {css}
    <div class="com-table-wrap">
      <table class="com-table">
        <colgroup>
          <col class="c1"/><col class="c2"/><col class="c3"/>
        </colgroup>
        <thead>
          <tr><th>community_id</th><th>size</th><th>nodes</th></tr>
        </thead>
        <tbody>
          {''.join(body)}
        </tbody>
      </table>
    </div>
    """

    # height dinamis agar tidak kepotong
    height = min(800, 140 + 52 * len(rows))
    components.html(html, height=height, scrolling=True)


def render_communities_paginated(res):
    """Render 'Daftar komunitas' section with pagination (same behavior as before)."""
    st.subheader("Daftar komunitas")

    comms = (res or {}).get("communities", []) or []
    n = len(comms)
    if n == 0:
        st.info("Tidak ada komunitas untuk ditampilkan.")
        return

    if "page_size" not in st.session_state:
        st.session_state.page_size = 25
    if "page" not in st.session_state:
        st.session_state.page = 1

    page_size = st.session_state.page_size
    pages = max(1, math.ceil(n / page_size))

    st.session_state.page = max(1, min(st.session_state.page, pages))
    page = st.session_state.page

    start = (page - 1) * page_size
    end = min(start + page_size, n)

    memberships = (res or {}).get("memberships", {}) or {}
    overlap_set = {node for node, k in memberships.items() if k > 1}
    st.caption("Node overlap ditandai warna berbeda.")

    rows = []
    for i in range(start, end):
        nodes = sorted(list(comms[i]))
        rows.append({"community_id": i + 1, "size": len(nodes), "nodes": nodes})

    render_communities_table(rows, overlap_set)

    left, right = st.columns([3, 2], vertical_alignment="center")

    with left:
        st.caption(f"{start+1} to {end} of {n}")

    with right:
        c1, c2, c3, c4, c5 = st.columns([1, 1, 1, 1, 2], vertical_alignment="center")

        if c1.button("⏮", help="First", use_container_width=True):
            st.session_state.page = 1
            st.rerun()

        if c2.button("◀", help="Prev", use_container_width=True):
            st.session_state.page = max(1, st.session_state.page - 1)
            st.rerun()

        if c3.button("▶", help="Next", use_container_width=True):
            st.session_state.page = min(pages, st.session_state.page + 1)
            st.rerun()

        if c4.button("⏭", help="Last", use_container_width=True):
            st.session_state.page = pages
            st.rerun()

        new_page = c5.number_input(
            f"Page of {pages}",
            min_value=1,
            max_value=pages,
            value=st.session_state.page,
            step=1,
            label_visibility="collapsed",
        )
        if new_page != st.session_state.page:
            st.session_state.page = int(new_page)
            st.rerun()
