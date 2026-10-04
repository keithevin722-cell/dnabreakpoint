# dnabreakpoint
build an ai that can do the research from published material and give you the newest options as links through a web interface where you can paste your bad genome and it can see and understand and research where to break the dna that has the best outcomes based on your uploaded problem, should include virus delivery.

## Usage
`python -m dnabreakpoint.server [port]` then open http://127.0.0.1:8000, paste a sequence (FASTA/raw)
with the gene, rsID or HGVS variant. It suggests Cas9 cut-site candidates and searches Europe PMC
(newest first) for CRISPR/base/prime editing and AAV/lentiviral delivery research. Results are saved in a
SQLite knowledge base (`knowledge.db`, override with `KB_PATH`) and reused/refreshed weekly.
No-install version: just open `index.html` in any browser (works offline for cut sites; research needs internet; knowledge base saved in browser localStorage).
Sequence identification is available in the Python server: opt in to send the pasted sequence to NCBI BLAST; matches include identity and query coverage. The standalone `index.html` does not submit DNA to external search services.
In Codespaces, start the server and open the **Ports** tab to use its forwarded URL. Keep port visibility **Private**; the app has no authentication. Set `HOST=127.0.0.1` to disable network-interface binding.
Tests: `python -m unittest discover tests`. Research tool only; not medical advice.
