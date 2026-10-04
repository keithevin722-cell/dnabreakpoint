import unittest
from dnabreakpoint.analysis import analyze
from dnabreakpoint.research import KnowledgeBase, evidence, research

FAKE = [{"id": "MED:1", "title": "New AAV CRISPR", "authors": "A", "journal": "J",
         "date": "2026-01-01", "year": "2026", "url": "u"}]


class T(unittest.TestCase):
    def test_analyze(self):
        a = analyze(">CFTR rs113993960\nATGCGTACGTTAGCCTAGGCTAGCTAGGATCGATCGGATC\n")
        self.assertEqual(a["terms"]["genes"], ["CFTR"])
        self.assertEqual(a["terms"]["rsids"], ["rs113993960"])
        self.assertTrue(a["cut_sites"])

    def test_vcf(self):
        vcf = ("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
               "7\t117559590\trs113993960\tATCT\tA\t.\tPASS\tGENEINFO=CFTR:1080\n"
               "1\t100\t.\tG\tA\t.\tPASS\tANN=A|missense_variant|MODERATE|BRCA1|x\n")
        a = analyze(vcf)
        self.assertEqual(len(a["variants"]), 2)
        self.assertEqual(a["terms"]["rsids"], ["rs113993960"])
        self.assertEqual(a["terms"]["genes"], ["BRCA1", "CFTR"])

    def test_evidence(self):
        kb = KnowledgeBase(":memory:", lambda q: [
            {"id": "MED:1", "title": "AAV9 capsid delivery of CRISPR base editor to exon 5", "authors": "",
             "journal": "", "date": "2026-01-01", "year": "2026", "url": "u"}])
        ev = evidence(research(kb, analyze(">CFTR\nATGC")["terms"]))
        self.assertEqual(ev["delivery"][0]["name"], "AAV (adeno-associated virus)")
        self.assertEqual(ev["delivery"][0]["count"], 1)
        self.assertTrue(any(r["name"] == "Base editing" for r in ev["editor"]))
        self.assertEqual(ev["target"][0]["name"], "Exon")

    def test_variant_proximity(self):
        seq = "ATGCGTACGTTAGCCTAGGCTAGCTAGGATCGATCGGATC"
        a = analyze(f">CFTR pos=3\n{seq}")
        self.assertEqual(a["variant_pos"], 2)
        d = [s["distance"] for s in a["cut_sites"]]
        self.assertEqual(d, sorted(d))
        b = analyze(">CFTR\n" + seq[:10] + seq[10].lower() + seq[11:])
        self.assertEqual(b["variant_pos"], 10)
        self.assertIsNone(analyze(">CFTR\n" + seq)["variant_pos"])

    def test_abstract_evidence(self):
        kb = KnowledgeBase(":memory:", lambda q: [
            {"id": "MED:2", "title": "Gene therapy results", "abstract": "We used lentiviral vectors.",
             "authors": "", "journal": "", "date": "", "year": "", "url": "u"}])
        ev = evidence(research(kb, analyze(">CFTR\nATGC")["terms"]))
        self.assertEqual(ev["delivery"][0]["name"], "Lentiviral vector")

    def test_kb_learns(self):
        calls = []
        kb = KnowledgeBase(":memory:", lambda q: calls.append(q) or FAKE)
        terms = analyze(">CFTR\nATGC")["terms"]
        r = research(kb, terms)
        research(kb, terms)
        self.assertEqual(len(calls), 2)  # cached on second run
        self.assertEqual(r["editing"][0]["title"], "New AAV CRISPR")
        self.assertEqual(kb.stats()["papers"], 1)


if __name__ == "__main__":
    unittest.main()
