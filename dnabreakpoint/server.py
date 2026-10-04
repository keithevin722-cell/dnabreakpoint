"""Web interface: paste a genome, get cut-site candidates and the newest research."""
import html
import os
import sys
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs

from .analysis import analyze
from .research import KnowledgeBase, research

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
<textarea name=genome id=g required></textarea><p><label><input type=checkbox name=refresh value=1>
Force fresh research</label> <button>Research treatments</button></p></form>
<p>Knowledge base: {papers} papers from {queries} searches.</p>{results}
<script>
document.getElementById("up").onclick=()=>document.getElementById("f").click();
document.getElementById("f").onchange=async e=>{{const f=e.target.files[0];if(!f)return;
document.getElementById("fn").textContent=f.name;
let b=f.stream();if(/\\.gz$/i.test(f.name))b=b.pipeThrough(new DecompressionStream("gzip"));
document.getElementById("g").value=await new Response(b).text()}};
</script>"""


def render_results(a, res):
    e = html.escape
    t = a["terms"]
    out = [f"<h2>Analysis</h2><p>Sequence length: {a['length']} bp, GC: {a['gc']}</p>",
           f"<p>Detected: genes {e(', '.join(t['genes']) or '-')}; rsIDs {e(', '.join(t['rsids']) or '-')}; "
           f"variants {e(', '.join(t['variants']) or '-')}</p>"]
    if a.get("variants"):
        out.append(f"<h2>Variants from VCF ({len(a['variants'])})</h2><table>"
                   "<tr><th>Chrom<th>Pos<th>ID<th>Ref<th>Alt<th>Genes</tr>")
        for v in a["variants"][:50]:
            out.append(f"<tr><td>{e(v['chrom'])}<td>{e(v['pos'])}<td>{e(v['id'])}<td>{e(v['ref'][:30])}"
                       f"<td>{e(v['alt'][:30])}<td>{e(', '.join(v['genes']))}</tr>")
        out.append("</table>")
    if a["cut_sites"]:
        out.append("<h2>Candidate Cas9 cut sites</h2><table><tr><th>Strand<th>Guide<th>PAM<th>Cut pos<th>GC</tr>")
        for s in a["cut_sites"]:
            out.append(f"<tr><td>{s['strand']}<td><code>{s['guide']}</code><td>{s['pam']}"
                       f"<td>{s['cut_position']}<td>{s['gc']}</tr>")
        out.append("</table>")
    if not res:
        out.append("<p>No gene, rsID or variant detected; add one to search the literature.</p>")
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
            a = analyze(form.get("genome", [""])[0])
            res = research(kb, a["terms"], refresh=bool(form.get("refresh")))
            self._page(render_results(a, res))
    return Handler


def main():
    port = int(os.environ.get("PORT", sys.argv[1] if len(sys.argv) > 1 else 8000))
    kb = KnowledgeBase(os.environ.get("KB_PATH", "knowledge.db"))
    print(f"Serving on http://127.0.0.1:{port}")
    ThreadingHTTPServer(("127.0.0.1", port), make_handler(kb)).serve_forever()


if __name__ == "__main__":
    main()
