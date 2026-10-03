"""Parse pasted genome input and extract research topics and candidate cut sites."""
import re

GENE_RE = re.compile(r"\b[A-Z][A-Z0-9]{2,9}\b")
RSID_RE = re.compile(r"\brs\d{3,}\b", re.I)
HGVS_RE = re.compile(r"\b(?:NM_|NC_|NG_)\d+(?:\.\d+)?:[cgp]\.[^\s,;]+")
STOP = {"DNA", "RNA", "THE", "AND", "FOR", "WITH", "NOT", "GENE", "SNP", "PAM"}
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def parse_input(text):
    """Split pasted text into (nucleotide sequence, free-text annotations)."""
    seq, notes = [], []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith((">", "#", ";")):
            notes.append(line.lstrip(">#; "))
            continue
        cleaned = re.sub(r"[\s\d]", "", line).upper()
        if len(cleaned) >= 10 and set(cleaned) <= set("ACGTUN"):
            seq.append(cleaned.replace("U", "T"))
        else:
            notes.append(line)
    return "".join(seq), "\n".join(notes)


def extract_terms(notes):
    genes = sorted({g for g in GENE_RE.findall(notes) if g not in STOP and not g.startswith(("NM", "NC", "NG"))})
    return {
        "genes": genes[:5],
        "rsids": sorted({r.lower() for r in RSID_RE.findall(notes)})[:5],
        "variants": sorted(set(HGVS_RE.findall(notes)))[:5],
    }


def find_cut_sites(seq, guide_len=20, limit=10):
    """Find SpCas9 NGG PAM sites on both strands; rank by GC content nearest 50%."""
    sites = []
    rc = seq.translate(COMPLEMENT)[::-1]
    for strand, s in (("+", seq), ("-", rc)):
        for i in range(guide_len, len(s) - 2):
            if s[i + 1:i + 3] == "GG":
                guide = s[i - guide_len:i]
                if "N" in guide or "TTTT" in guide:
                    continue
                gc = (guide.count("G") + guide.count("C")) / guide_len
                pos = i - 3 if strand == "+" else len(seq) - i + 3
                sites.append({"strand": strand, "guide": guide, "pam": s[i:i + 3],
                              "cut_position": pos, "gc": round(gc, 2)})
    sites.sort(key=lambda x: abs(x["gc"] - 0.5))
    return sites[:limit]


def analyze(text):
    seq, notes = parse_input(text)
    terms = extract_terms(notes)
    gc = round((seq.count("G") + seq.count("C")) / len(seq), 3) if seq else None
    return {"length": len(seq), "gc": gc, "notes": notes, "terms": terms,
            "cut_sites": find_cut_sites(seq) if seq else []}
