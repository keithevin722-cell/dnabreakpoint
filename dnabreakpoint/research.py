"""Literature research (Europe PMC) with a SQLite knowledge base that grows over time."""
import json
import re
import sqlite3
import time
import urllib.parse
import urllib.request

API = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
TOPICS = {
    "editing": "(CRISPR OR \"base editing\" OR \"prime editing\") AND (therapy OR treatment)",
    "delivery": "(AAV OR lentiviral OR \"viral vector\" OR \"adeno-associated\") AND (delivery OR gene therapy)",
}
CACHE_SECONDS = 7 * 86400


def fetch_europepmc(query, page_size=15):
    params = urllib.parse.urlencode({"query": query + " sort_date:y", "format": "json",
                                     "pageSize": page_size, "resultType": "lite"})
    with urllib.request.urlopen(f"{API}?{params}", timeout=20) as r:
        data = json.load(r)
    out = []
    for p in data.get("resultList", {}).get("result", []):
        pid = p.get("pmid") or p.get("id")
        src = p.get("source", "MED")
        out.append({"id": f"{src}:{p.get('id')}", "title": p.get("title", ""),
                    "authors": p.get("authorString", ""), "journal": p.get("journalTitle", ""),
                    "date": p.get("firstPublicationDate", ""), "year": p.get("pubYear", ""),
                    "url": f"https://europepmc.org/article/{src}/{p.get('id')}" if pid else ""})
    return out


class KnowledgeBase:
    def __init__(self, path="knowledge.db", fetcher=fetch_europepmc):
        self.db = sqlite3.connect(path, check_same_thread=False)
        self.fetcher = fetcher
        self.db.executescript("""
        CREATE TABLE IF NOT EXISTS papers(id TEXT PRIMARY KEY, title TEXT, authors TEXT,
            journal TEXT, date TEXT, year TEXT, url TEXT);
        CREATE TABLE IF NOT EXISTS queries(query TEXT PRIMARY KEY, last_run REAL, runs INTEGER);
        CREATE TABLE IF NOT EXISTS query_papers(query TEXT, paper_id TEXT, PRIMARY KEY(query, paper_id));
        """)

    def search(self, query, refresh=False):
        """Return papers for a query newest first; learn new papers from the web when stale."""
        row = self.db.execute("SELECT last_run FROM queries WHERE query=?", (query,)).fetchone()
        stale = refresh or row is None or time.time() - row[0] > CACHE_SECONDS
        if stale:
            try:
                for p in self.fetcher(query):
                    self.db.execute("INSERT OR REPLACE INTO papers VALUES(?,?,?,?,?,?,?)",
                                    (p["id"], p["title"], p["authors"], p["journal"],
                                     p["date"], p["year"], p["url"]))
                    self.db.execute("INSERT OR IGNORE INTO query_papers VALUES(?,?)", (query, p["id"]))
                self.db.execute("INSERT INTO queries VALUES(?,?,1) ON CONFLICT(query) DO UPDATE "
                                "SET last_run=excluded.last_run, runs=runs+1", (query, time.time()))
                self.db.commit()
            except OSError:
                pass  # offline: fall back to what has been learned already
        cur = self.db.execute("SELECT p.id,p.title,p.authors,p.journal,p.date,p.year,p.url FROM papers p "
                              "JOIN query_papers q ON q.paper_id=p.id WHERE q.query=? "
                              "ORDER BY p.date DESC", (query,))
        keys = ("id", "title", "authors", "journal", "date", "year", "url")
        return [dict(zip(keys, r)) for r in cur]

    def stats(self):
        return {"papers": self.db.execute("SELECT COUNT(*) FROM papers").fetchone()[0],
                "queries": self.db.execute("SELECT COUNT(*) FROM queries").fetchone()[0]}


def build_queries(terms):
    subject = terms["genes"] + terms["rsids"] + [f'"{v}"' for v in terms["variants"]]
    if not subject:
        return {}
    subj = "(" + " OR ".join(subject) + ")"
    return {name: f"{subj} AND {q}" for name, q in TOPICS.items()}


def research(kb, terms, refresh=False):
    return {name: kb.search(q, refresh) for name, q in build_queries(terms).items()}


SIGNALS = {
    "delivery": {
        "AAV (adeno-associated virus)": r"\bAAV\d*|adeno-associated",
        "Engineered / synthetic capsid or vector": r"capsid|synthetic|engineered (?:AAV|vector|virus)|directed evolution",
        "Lentiviral vector": r"lentivir",
        "Adenoviral vector": r"adenovir",
        "Lipid nanoparticle (non-viral)": r"lipid nanoparticle|\bLNP",
        "Ex vivo / electroporation": r"ex vivo|electroporat",
    },
    "editor": {
        "Cas9 nuclease": r"\bCas9\b|CRISPR",
        "Base editing": r"base[- ]edit|\b[AC]BE\b",
        "Prime editing": r"prime[- ]edit",
        "Cas12 / Cas13": r"Cas12|Cas13",
    },
    "target": {
        "Exon": r"\bexon",
        "Exon skipping / splice site": r"exon skipping|splic",
        "Intron": r"intron",
        "Promoter / enhancer": r"promoter|enhancer",
    },
}


def evidence(res, examples=3):
    """Rank delivery/editor/target themes by how many retrieved paper titles mention them."""
    papers = {}
    for plist in res.values():
        for p in plist:
            papers[p["id"]] = p
    out = {}
    for group, pats in SIGNALS.items():
        rows = []
        for name, pat in pats.items():
            hits = [p for p in papers.values() if re.search(pat, p["title"], re.I)]
            if hits:
                rows.append({"name": name, "count": len(hits), "papers": hits[:examples]})
        out[group] = sorted(rows, key=lambda r: -r["count"])
    return out
