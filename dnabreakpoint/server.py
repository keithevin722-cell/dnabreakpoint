"""Web interface: paste a genome, get cut-site candidates and the newest research."""
import html
import os
import sys
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs

from .analysis import analyze, parse_input
from .research import KnowledgeBase, evidence, identify_sequence, research

MAX_BODY = 20_000_000
DISCLAIMER = ("Research tool only. Not medical advice; candidate guides are computational "
              "suggestions that need expert and clinical validation.")
PAGE = """<!doctype html><meta charset=utf-8><title>DNA Breakpoint</title>
<style>body{{font-family:sans-serif;max-width:900px;margin:2em auto;padding:0 1em}}
textarea{{width:100%;height:14em;font-family:monospace}}td,th{{padding:2px 8px;text-align:left}}</style>
<h1>DNA Breakpoint</h1><p><em>{disclaimer}</em></p>
<form method=post><p>Paste your genome / sequence (FASTA or raw). Add the gene, rsID or HGVS
variant (e.g. CFTR, rs113993960, NM_000492.4:c.1521_1523del) on header or extra lines, plus a description of the problem.</p>
<p><button type=button id=up>Upload file</button> <span id=fn>VCF, .vcf.gz, FASTA or text</span>
<input type=file id=f hidden accept=".vcf,.gz,.fa,.fasta,.txt"></p>
<textarea name=genome id=g required></textarea><p><label><input type=checkbox name=identify value=1>
Identify sequence with NCBI BLAST (sends DNA sequence to NCBI; 15-20,000 bases)</label></p><p><label><input type=checkbox name=refresh value=1>
Force fresh research</label> <button>Research treatments</button></p></form>
<p>Knowledge base: {papers} papers from {queries} searches.</p>{results}
<script>
document.getElementById("up").onclick=()=>document.getElementById("f").click();
document.getElementById("f").onchange=async e=>{{const f=e.target.files[0];if(!f)return;
document.getElementById("fn").textContent=f.name;
let b=f.stream();if(/\\.gz$/i.test(f.name))b=b.pipeThrough(new DecompressionStream("gzip"));
document.getElementById("g").value=await new Response(b).text()}};
</script>"""


def render_results(a, res, identification=None):
    e = html.escape
    t = a["terms"]
    out = [f"<h2>Analysis</h2><p>Sequence length: {a['length']} bp, GC: {a['gc']}</p>",
           f"<p>Detected: genes {e(', '.join(t['genes']) or '-')}; rsIDs {e(', '.join(t['rsids']) or '-')}; "
           f"variants {e(', '.join(t['variants']) or '-')}</p>"]
    if identification:
        statuses = {"invalid": "Enter at least 15 DNA bases using A, C, G, T or N.",
                    "too_long": "Sequence search is limited to 20,000 bases.",
                    "no_hits": "NCBI BLAST did not return a matching sequence.",
                    "timeout": "NCBI BLAST did not finish in time. Try again later.",
                    "unavailable": "NCBI BLAST is currently unavailable.",
                    "ok": "Matches are database similarities, not a definitive gene identification. Short sequences can match many regions; use a longer sequence for more confidence."}
        out.append("<h2>Sequence identification (NCBI BLAST)</h2><p><small>" + statuses[identification["status"]] + "</small></p>")
        if identification["status"] == "ok":
            if identification["hits"]:
                out.append("<ol>")
                for hit in identification["hits"]:
                    out.append(f"<li><a href='{e(hit['url'])}' rel=noopener>{e(hit['title'])}</a> "
                               f"<small>{hit['identity']}% identity; {hit['coverage']}% query coverage</small></li>")
                out.append("</ol>")
            else:
                out.append("<p>No sequence matches returned.</p>")
    if a.get("variants"):
        out.append(f"<h2>Variants from VCF ({len(a['variants'])})</h2><table>"
                   "<tr><th>Chrom<th>Pos<th>ID<th>Ref<th>Alt<th>Genes</tr>")
        for v in a["variants"][:50]:
            out.append(f"<tr><td>{e(v['chrom'])}<td>{e(v['pos'])}<td>{e(v['id'])}<td>{e(v['ref'][:30])}"
                       f"<td>{e(v['alt'][:30])}<td>{e(', '.join(v['genes']))}</tr>")
        out.append("</table>")
    if a["cut_sites"]:
        vp = a.get("variant_pos") is not None
        out.append("<h2>Candidate Cas9 cut sites</h2>"
                   + (f"<p><small>Ranked by distance to the variant at position {a['variant_pos'] + 1}, then GC.</small></p>" if vp
                      else "<p><small>No variant position given (add 'pos=N' or mark the variant in lowercase), so ranked by GC only.</small></p>")
                   + "<table><tr><th>Strand<th>Guide<th>PAM<th>Cut pos<th>GC" + ("<th>Dist" if vp else "") + "</tr>")
        for s in a["cut_sites"]:
            out.append(f"<tr><td>{s['strand']}<td><code>{s['guide']}</code><td>{s['pam']}"
                       f"<td>{s['cut_position']}<td>{s['gc']}" + (f"<td>{s['distance']}" if vp else "") + "</tr>")
        out.append("</table>")
    if not res:
        out.append("<p>No gene, rsID or variant detected; add one to search the literature.</p>")
    else:
        titles = {"delivery": "Delivery options in the research", "editor": "Editing approaches in the research",
                  "target": "Target-region themes in the research"}
        out.append("<h2>Research summary</h2><p><small>Ranked by how many of the newest papers mention each "
                   "item in the title or abstract. This shows what the literature is discussing, not a clinical recommendation.</small></p>")
        for group, rows in evidence(res).items():
            if not rows:
                continue
            out.append(f"<h3>{titles[group]}</h3><ol>")
            for r in rows:
                links = " ".join(f"<a href='{e(p['url'])}' rel=noopener title='{e(p['title'])}'>[{i + 1}]</a>"
                                 for i, p in enumerate(r["papers"]))
                out.append(f"<li>{e(r['name'])}: {r['count']} papers {links}</li>")
            out.append("</ol>")
    for name, papers in res.items():
        out.append(f"<h2>Latest research: {e(name)}</h2><ul>")
        for p in papers:
            out.append(f"<li><a href='{e(p['url'])}' rel=noopener>{e(p['title'])}</a> "
                       f"<small>{e(p['journal'])} {e(p['date'] or p['year'])}</small></li>")
        out.append("</ul>" if papers else "<li>No results yet.</li></ul>")
    return "".join(out)


def make_handler(kb):
    class Handler(BaseHTTPRequestHandler):
        def _send(self, body, code=200):
            data = body.encode()
            self.send_response(code)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)

        def _page(self, results=""):
            s = kb.stats()
            self._send(PAGE.format(disclaimer=DISCLAIMER, results=results, **s))

        def do_GET(self):
            self._page()

        def do_POST(self):
            n = int(self.headers.get("Content-Length") or 0)
            if n > MAX_BODY:
                return self._send("Too large", 413)
            form = parse_qs(self.rfile.read(n).decode("utf-8", "replace"))
            genome = form.get("genome", [""])[0]
            a = analyze(genome)
            res = research(kb, a["terms"], refresh=bool(form.get("refresh")))
            identification = None
            if form.get("identify") and not a.get("variants"):
                sequence, _ = parse_input(genome)
                identification = identify_sequence(sequence)
            self._page(render_results(a, res, identification))
    return Handler


def main():
    port = int(os.environ.get("PORT", sys.argv[1] if len(sys.argv) > 1 else 8000))
    kb = KnowledgeBase(os.environ.get("KB_PATH", "knowledge.db"))
    print(f"Serving on http://127.0.0.1:{port}")
    ThreadingHTTPServer(("127.0.0.1", port), make_handler(kb)).serve_forever()


if __name__ == "__main__":
    main()
