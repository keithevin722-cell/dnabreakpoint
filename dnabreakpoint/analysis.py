"""Parse pasted genome input and extract research topics and candidate cut sites."""
import re

GENE_RE = re.compile(r"\b[A-Z][A-Z0-9]{2,9}\b")
RSID_RE = re.compile(r"\brs\d{3,}\b", re.I)
HGVS_RE = re.compile(r"\b(?:NM_|NC_|NG_)\d+(?:\.\d+)?:[cgp]\.[^\s,;]+")
STOP = {"DNA", "RNA", "THE", "AND", "FOR", "WITH", "NOT", "GENE", "SNP", "PAM"}
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _parse(text):
    seq, notes, mask, offset = [], [], [], 0
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith((">", "#", ";")):
            notes.append(line.lstrip(">#; "))
            continue
        raw = re.sub(r"[\s\d-]", "", line)
        cleaned = raw.upper()
        if len(cleaned) >= 10 and set(cleaned) <= set("ACGTUN"):
            mask += [offset + i for i, c in enumerate(raw) if c.islower()]
            offset += len(cleaned)
            seq.append(cleaned.replace("U", "T"))
        else:
            notes.append(line)
    return "".join(seq), "\n".join(notes), mask


def parse_input(text):
    """Split pasted text into (nucleotide sequence, free-text annotations)."""
    seq, notes, _ = _parse(text)
    return seq, notes


POS_RE = re.compile(r"\b(?:variant_)?pos(?:ition)?\s*[=:]\s*(\d+)", re.I)


def variant_position(notes, mask, seq_len):
    """0-based variant index from a 'pos=N' (1-based) note or a short lowercase-marked region."""
    m = POS_RE.search(notes)
    if m and 0 < int(m.group(1)) <= seq_len:
        return int(m.group(1)) - 1
    if mask and len(mask) <= 50 and len(mask) < seq_len:
        return mask[len(mask) // 2]
    return None


def extract_terms(notes):
    genes = sorted({g for g in GENE_RE.findall(notes) if g not in STOP and not g.startswith(("NM", "NC", "NG"))})
    return {
        "genes": genes[:5],
        "rsids": sorted({r.lower() for r in RSID_RE.findall(notes)})[:5],
        "variants": sorted(set(HGVS_RE.findall(notes)))[:5],
    }


def find_cut_sites(seq, guide_len=20, limit=10, variant_pos=None):
    """Find SpCas9 NGG PAM sites on both strands; rank by distance to the variant, then GC nearest 50%."""
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
                              "cut_position": pos, "gc": round(gc, 2),
                              "distance": None if variant_pos is None else abs(pos - variant_pos)})
    sites.sort(key=lambda x: (x["distance"] if variant_pos is not None else 0, abs(x["gc"] - 0.5)))
    return sites[:limit]


def is_vcf(text):
    head = text.lstrip()[:2000]
    return head.startswith("##fileformat=VCF") or "\n#CHROM" in "\n" + head


def _info_genes(info):
    genes = []
    for item in info.split(";"):
        key, _, val = item.partition("=")
        if key in ("GENE", "SYMBOL"):
            genes += val.split(",")
        elif key == "GENEINFO":
            genes += [g.split(":")[0] for g in val.split("|")]
        elif key == "ANN":  # SnpEff: allele|effect|impact|gene|...
            for ann in val.split(","):
                f = ann.split("|")
                if len(f) > 3:
                    genes.append(f[3])
    return genes


def parse_vcf(text):
    """Return (variants, terms) from VCF text; terms match extract_terms output."""
    variants, genes, rsids = [], set(), set()
    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        f = line.rstrip("\r\n").split("\t")
        if len(f) < 5:
            f = line.split()
        if len(f) < 5:
            continue
        chrom, pos, vid, ref, alt = f[:5]
        info = f[7] if len(f) > 7 else ""
        ids = [] if vid == "." else vid.split(";")
        rsids.update(i.lower() for i in ids if RSID_RE.fullmatch(i))
        vg = [g for g in _info_genes(info) if g and g != "."]
        genes.update(g for g in vg if GENE_RE.fullmatch(g) and g not in STOP)
        variants.append({"chrom": chrom, "pos": pos, "id": vid, "ref": ref, "alt": alt,
                         "genes": sorted(set(vg))})
    terms = {"genes": sorted(genes)[:5], "rsids": sorted(rsids)[:5], "variants": []}
    return variants, terms


def analyze(text):
    if is_vcf(text):
        variants, terms = parse_vcf(text)
        return {"length": 0, "gc": None, "notes": "", "terms": terms, "cut_sites": [],
                "variants": variants}
    seq, notes, mask = _parse(text)
    terms = extract_terms(notes)
    vpos = variant_position(notes, mask, len(seq))
    gc = round((seq.count("G") + seq.count("C")) / len(seq), 3) if seq else None
    return {"length": len(seq), "gc": gc, "notes": notes, "terms": terms, "variant_pos": vpos,
            "cut_sites": find_cut_sites(seq, variant_pos=vpos) if seq else [], "variants": []}
