import math
import re
from urllib.parse import quote_plus

import pandas as pd
import requests
import streamlit as st

from core.enrichment import run_enrichr, ENRICHR_FIXED_LIBS


_GO_ID_RE = re.compile(r"(GO:\d{7})")
_HSA_ID_RE = re.compile(r"(hsa\d{5})", re.IGNORECASE)
_MAP_ID_RE = re.compile(r"(map\d{5})", re.IGNORECASE)
_KEGG_NUM_RE = re.compile(r"(?:KEGG:)?(\d{5})", re.IGNORECASE)


@st.cache_data(show_spinner=False, ttl=60 * 60 * 24)
def _kegg_find_pathway_id(term_name: str, organism_prefix: str = "hsa") -> str:
    """Resolve KEGG pathway id (e.g., hsa04931) from a pathway name using KEGG REST."""
    name = str(term_name or "").strip()
    if not name:
        return ""
    try:
        url = "https://rest.kegg.jp/find/pathway/" + quote_plus(name)
        r = requests.get(url, timeout=15)
        if r.status_code != 200:
            return ""
        lines = [ln for ln in r.text.splitlines() if ln.strip()]
        for ln in lines:
            left = ln.split("\t", 1)[0]
            if left.startswith("path:"):
                pid = left.split(":", 1)[1].strip()
                if pid.startswith(organism_prefix):
                    return pid
        for ln in lines:
            left = ln.split("\t", 1)[0]
            if left.startswith("path:"):
                return left.split(":", 1)[1].strip()
    except Exception:
        return ""
    return ""


def _parse_term_id_and_name(term: str, lib_label: str):
    s = str(term or "").strip()
    if not s:
        return "", ""

    if lib_label == "KEGG Pathway":
        m = _HSA_ID_RE.search(s)
        if m:
            tid = m.group(1).lower()
            name = re.sub(_HSA_ID_RE, "", s).strip()
        else:
            m2 = _MAP_ID_RE.search(s)
            m3 = _KEGG_NUM_RE.search(s)
            if m2:
                tid = m2.group(1).lower()
                name = re.sub(_MAP_ID_RE, "", s).strip()
            elif m3:
                code = m3.group(1)
                tid = "hsa" + code
                name = re.sub(_KEGG_NUM_RE, "", s).strip()
            else:
                return "", s

        name = name.replace("Homo sapiens", "").replace("homo sapiens", "")
        name = name.strip(" -()[]")
        name = re.sub(r"\s{2,}", " ", name).strip()
        return tid, name

    m = _GO_ID_RE.search(s)
    if m:
        tid = m.group(1)
        name = re.sub(r"\s*\(?"+re.escape(tid)+r"\)?\s*", " ", s).strip()
        name = name.strip(" -()[]")
        name = re.sub(r"\s{2,}", " ", name).strip()
        return tid, name

    return "", s


def _term_url(term_id: str, lib_label: str) -> str:
    if not term_id:
        return ""
    if lib_label == "KEGG Pathway":
        tid = str(term_id).strip().lower()
        if tid.startswith("kegg:"):
            tid = tid.split(":", 1)[1]
        if re.fullmatch(r"\d{5}", tid):
            tid = "hsa" + tid
        return f"https://www.genome.jp/entry/pathway+{tid}"
    return f"https://amigo.geneontology.org/amigo/term/{term_id}"


def _extract_genes(val):
    if val is None:
        return []
    try:
        import pandas as _pd
        if isinstance(val, float) and _pd.isna(val):
            return []
    except Exception:
        pass

    if isinstance(val, (list, tuple, set)):
        genes = [str(x).strip() for x in val if str(x).strip()]
    else:
        s = str(val).strip()
        if not s:
            return []
        if ";" in s:
            genes = [g.strip() for g in s.split(";") if g.strip()]
        elif "," in s:
            genes = [g.strip() for g in s.split(",") if g.strip()]
        else:
            genes = [g.strip() for g in s.split() if g.strip()]

    seen = set()
    out = []
    for g in genes:
        if g not in seen:
            seen.add(g)
            out.append(g)
    return out


def _render_enrichment_table(df_show: pd.DataFrame, lib_label: str, key_prefix: str) -> pd.DataFrame:
    if not isinstance(df_show, pd.DataFrame) or len(df_show) == 0:
        st.info("Tidak ada term yang lolos filter untuk library ini.")
        return pd.DataFrame(columns=["term_id", "term_name", "adj_p_value", "gene_count", "genes"])

    rows = []
    for _, r in df_show.iterrows():
        term_raw = (
            r.get("term")
            or r.get("Term")
            or r.get("term_name")
            or r.get("Term Name")
            or ""
        )
        term_id, term_name = _parse_term_id_and_name(term_raw, lib_label)
        if lib_label == "KEGG Pathway" and not term_id:
            pid = _kegg_find_pathway_id(term_name or term_raw, organism_prefix="hsa")
            if pid:
                term_id = pid

        genes_raw = (
            r.get("overlap_genes")
            or r.get("Genes")
            or r.get("genes")
            or r.get("overlapping_genes")
        )
        genes = _extract_genes(genes_raw)

        adj = (
            r.get("adj_p_value")
            or r.get("Adjusted.P.value")
            or r.get("adjusted_p_value")
            or r.get("adj_pvalue")
        )
        try:
            adj_f = float(adj)
        except Exception:
            adj_f = float("nan")

        rows.append(
            {"term_id": term_id, "term_name": term_name, "adj_p_value": adj_f, "gene_count": int(len(genes)), "genes": genes}
        )

    out = pd.DataFrame(rows)

    h = st.columns([1.5, 4.8, 1.5, 1.2, 1.6])
    h[0].markdown("**Term ID**")
    h[1].markdown("**Term Name**")
    h[2].markdown("**Adj p-value**")
    h[3].markdown("**Gene Count**")
    h[4].markdown("**Genes**")

    for i, row in out.iterrows():
        c = st.columns([1.5, 4.8, 1.5, 1.2, 1.6])
        url = _term_url(row["term_id"], lib_label)
        if row["term_id"] and url:
            c[0].markdown(f"[{row['term_id']}]({url})")
        else:
            c[0].write(row["term_id"] or "")

        c[1].write(row["term_name"] or "")
        if not (isinstance(row["adj_p_value"], float) and math.isnan(row["adj_p_value"])):
            c[2].write(f"{row['adj_p_value']:.3g}")
        else:
            c[2].write("")
        c[3].write(int(row["gene_count"]))

        state_key = f"{key_prefix}::open::{i}"
        btn_key = f"{key_prefix}::btn::{i}"
        if c[4].button("View genes", key=btn_key):
            st.session_state[state_key] = not st.session_state.get(state_key, False)

        if st.session_state.get(state_key, False):
            genes = row["genes"] or []
            if genes:
                st.markdown(" ".join([f"`{g}`" for g in genes]))
            else:
                st.caption("Tidak ada gene overlap yang tersedia.")
    return out


def render_enrichment(res, graph_label: str):
    st.subheader("Enrichment analysis (berdasarkan komunitas)")

    ENRICH_MAX_TERMS = 10
    ENRICH_ADJ_P_THR = 0.05
    ENRICH_SLEEP_S = 0.2

    if res is None:
        st.info("Jalankan deteksi komunitas dulu untuk melakukan enrichment.")
        return

    comms_all = res.get("communities", []) or []
    if len(comms_all) == 0:
        st.warning("Tidak ada komunitas untuk dianalisis.")
        return

    sizes = [(i + 1, len(c)) for i, c in enumerate(comms_all)]
    sizes_sorted = sorted(sizes, key=lambda x: x[1], reverse=True)

    st.caption(f"Total komunitas: **{len(comms_all)}** (diurutkan menurut ukuran).")

    with st.expander("Pengaturan enrichment", expanded=True):
        st.markdown(
            "Library tetap: **GO Biological Process**, **GO Molecular Function**, **GO Cellular Component**, **KEGG Pathway**."
        )
        st.caption(
            f"Parameter tetap: maksimal term = **{ENRICH_MAX_TERMS}**, ambang adj p-value = **{ENRICH_ADJ_P_THR}**, jeda request = **{ENRICH_SLEEP_S} detik**."
        )
        run_btn = st.button("Jalankan enrichment", type="primary")

    run_id = (res.get("meta") or {}).get("run_id", "run")

    if run_btn:
        st.session_state["ENRICH_LAST_RUN_ID"] = run_id

    show_results = (st.session_state.get("ENRICH_LAST_RUN_ID") == run_id)

    if not show_results:
        st.caption("Klik **Jalankan enrichment** untuk memanggil API dan menampilkan hasil.")
        return

    libs = list(ENRICHR_FIXED_LIBS.items())
    compute_if_missing = bool(run_btn)

    if compute_if_missing:
        total_jobs = max(1, len(sizes_sorted) * len(libs))
        progress = st.progress(0.0)
        job = 0

    comm_tab_labels = [f"C{cid} (n={sz})" for cid, sz in sizes_sorted]
    comm_tabs = st.tabs(comm_tab_labels)

    for (cid, sz), comm_tab in zip(sizes_sorted, comm_tabs):
        with comm_tab:
            genes = sorted(list(comms_all[cid - 1]))
            desc = f"{graph_label} | {res.get('algo_name','algo')} | C{cid} (n={len(genes)})"
            st.caption(desc)

            tab_labels = [x[0] for x in libs]
            tab_objs = st.tabs(tab_labels)

            community_frames = []

            for (lib_label, lib_id), tab in zip(libs, tab_objs):
                with tab:
                    cache_key = f"enrich::{run_id}::enrichr::{lib_id}::C{cid}"

                    if (cache_key not in st.session_state) and compute_if_missing:
                        try:
                            with st.spinner(f"Enrichment C{cid} • {lib_label} ..."):
                                df_raw = run_enrichr(
                                    genes,
                                    library=lib_id,
                                    description=desc,
                                    sleep_s=float(ENRICH_SLEEP_S),
                                )
                            st.session_state[cache_key] = df_raw
                        except Exception as e:
                            st.error(f"Gagal enrichment C{cid} • {lib_label}: {type(e).__name__}: {e}")
                            st.session_state[cache_key] = pd.DataFrame()

                        if compute_if_missing:
                            job += 1
                            progress.progress(min(job / total_jobs, 1.0))

                    df_raw = st.session_state.get(cache_key, pd.DataFrame())

                    if isinstance(df_raw, pd.DataFrame) and len(df_raw) > 0:
                        df_show = df_raw.copy()
                        if "adj_p_value" in df_show.columns:
                            df_show = df_show[df_show["adj_p_value"] <= float(ENRICH_ADJ_P_THR)]
                            df_show = df_show.sort_values(["adj_p_value", "p_value"], ascending=True)
                        elif "p_value" in df_show.columns:
                            df_show = df_show.sort_values(["p_value"], ascending=True)
                        df_show = df_show.head(int(ENRICH_MAX_TERMS))
                    else:
                        df_show = df_raw if isinstance(df_raw, pd.DataFrame) else pd.DataFrame(df_raw)

                    df_clean = _render_enrichment_table(df_show, lib_label=lib_label, key_prefix=cache_key)

                    if isinstance(df_clean, pd.DataFrame) and len(df_clean) > 0:
                        df_ctx = df_clean.copy()
                        df_ctx.insert(0, "community_id", cid)
                        df_ctx.insert(1, "community_size", sz)
                        df_ctx.insert(2, "provider", "Enrichr")
                        df_ctx.insert(3, "library", lib_label)
                        community_frames.append(df_ctx)

            if len(community_frames) > 0:
                df_comm = pd.concat(community_frames, ignore_index=True)
                if "genes" in df_comm.columns:
                    df_comm["genes"] = df_comm["genes"].apply(
                        lambda x: ";".join([str(g) for g in x]) if isinstance(x, (list, tuple, set)) else x
                    )
                csv = df_comm.to_csv(index=False).encode("utf-8")
                st.download_button(
                    f"Download CSV • C{cid}",
                    data=csv,
                    file_name=f"enrichment_{graph_label.replace(' ', '_')}_{run_id}_C{cid}.csv",
                    mime="text/csv",
                )

    if compute_if_missing:
        progress.empty()
