import unittest
from dnabreakpoint.analysis import analyze
from dnabreakpoint.research import KnowledgeBase, research

FAKE = [{"id": "MED:1", "title": "New AAV CRISPR", "authors": "A", "journal": "J",
         "date": "2026-01-01", "year": "2026", "url": "u"}]


class T(unittest.TestCase):
    def test_analyze(self):
        a = analyze(">CFTR rs113993960\nATGCGTACGTTAGCCTAGGCTAGCTAGGATCGATCGGATC\n")
        self.assertEqual(a["terms"]["genes"], ["CFTR"])
        self.assertEqual(a["terms"]["rsids"], ["rs113993960"])
        self.assertTrue(a["cut_sites"])

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
