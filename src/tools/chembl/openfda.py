# src/tools/chembl/openfda.py
# ===========================================================
# OpenFDA Black-Box Warning Fetcher for ChEMBL Molecules
# ===========================================================
#
# Author: Dmitri Kosenkov
# Created: 2025-10-09
# Updated: 2025-10-09 (refined caching, fallback, DailyMed URL)
#
# Provides a simple API and CLI tool for retrieving FDA
# boxed (black-box) warning text for a given preferred drug
# name (pref_name) using the openFDA drug labeling API.
#

import requests
from typing import Optional, Tuple


# ------------------------------
# Core client
# ------------------------------
class OpenFDAClient:
    """Query openFDA drug labeling API for boxed warning text."""

    BASE_URL = "https://api.fda.gov/drug/label.json"
    DAILyMED_BASE = "https://dailymed.nlm.nih.gov/dailymed/lookup.cfm?setid={}"
    _cache = {}  # local static cache to avoid redundant queries

    def __init__(self, api_base: Optional[str] = None):
        self.api_base = api_base or self.BASE_URL

    def get_boxed_warning(
        self, search_term: str, limit: int = 1
    ) -> Optional[Tuple[str, str]]:
        """
        Query openFDA for boxed warning text by active ingredient or brand name.
        Returns tuple (text, url) or None if not found.
        """
        if not search_term:
            return None

        # --- Cached lookup
        if search_term in self._cache:
            return self._cache[search_term]

        queries = [
            f'openfda.generic_name:"{search_term}"',
            f'openfda.brand_name:"{search_term}"'
        ]

        for q in queries:
            params = {"search": q, "limit": limit}
            try:
                resp = requests.get(self.api_base, params=params, timeout=10)
                if resp.status_code == 404:
                    continue  # no entry for this query
                resp.raise_for_status()
                js = resp.json()
                results = js.get("results", [])
                for rec in results:
                    bw = rec.get("boxed_warning")
                    if not bw:
                        continue

                    setid = rec.get("set_id") or rec.get("id")
                    url = self.DAILyMED_BASE.format(setid) if setid else ""
                    text = "\n".join(bw) if isinstance(bw, list) else str(bw)

                    self._cache[search_term] = (text, url)
                    return (text, url)
            except Exception as e:
                print(f"[OpenFDA] Error for '{search_term}': {e}")
                continue

        self._cache[search_term] = None
        return None


# ------------------------------
# Public API
# ------------------------------
def fetch_black_box_text_by_name(pref_name: str) -> Optional[Tuple[str, str]]:
    """
    Given a preferred drug name (pref_name), query openFDA for black-box warning text.
    Returns tuple (text, url) or None if not found or lookup fails.
    """
    if not pref_name or not pref_name.strip():
        print("[OpenFDA] Empty or invalid pref_name provided.")
        return None

    name = pref_name.strip()
    client = OpenFDAClient()
    result = client.get_boxed_warning(name)

    if result:
        text, url = result
        print(f"[OpenFDA] Black-box warning found for '{name}'.")
        if url:
            print(f"[OpenFDA] Label link: {url}")
        return result
    else:
        print(f"[OpenFDA] No warning found for '{name}'.")
        return None


# ------------------------------
# CLI self-test
# ------------------------------
if __name__ == "__main__":
    # Example known FDA black-box warning: Haloperidol
    test_name = "Haloperidol"
    print(f"Testing openFDA boxed warning fetch for '{test_name}' ...")
    result = fetch_black_box_text_by_name(test_name)
    if result:
        warning, url = result
        print("\n=== Boxed Warning Text ===")
        print(warning[:800])  # show first 800 chars
        if len(warning) > 800:
            print("... [truncated] ...")
        if url:
            print("\nFDA Label URL:", url)
    else:
        print("No boxed warning found.")