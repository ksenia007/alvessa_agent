from __future__ import annotations

from collections import defaultdict
from typing import Dict, List, Set, Tuple
import json
import time

import requests

from src.config import DISGENET_API_KEY, DEBUG

def get_entrez(gobj):
    eid = getattr(gobj, "entrez_id", None)
    if not eid and hasattr(gobj, "get_entrez_id"):
        try: eid = gobj.get_entrez_id()
        except Exception: eid = None
    if not eid:
        eid = None
    return str(eid) if eid else None

def _fetch_predictions_DisGeNet(gene):

    if DEBUG:
        print(f"[DisGeNet] Fetching predictions for gene symbol: {gene.symbol}")

    full_disease_list = []

    params = {}
    params["gene_ncbi_id"] = get_entrez(gene)

    if not params["gene_ncbi_id"]:
        return full_disease_list

    HTTPheadersDict = {}
    HTTPheadersDict['Authorization'] = DISGENET_API_KEY
    HTTPheadersDict['accept'] = 'application/json'

    page_count = 0
    total_disease_count = 0

    while not page_count or len(full_disease_list)<total_disease_count:
        params["page_number"] = str(page_count)

        response = requests.get("https://api.disgenet.com/api/v1/gda/summary",\
                            params=params, headers=HTTPheadersDict, verify=False)
        
        if not response.ok:
            if response.status_code == 429:
                while response.ok is False:
                    print("You have reached a query limit for your user. Please wait {} seconds until next query".format(\
                    response.headers['x-rate-limit-retry-after-seconds']))
                    time.sleep(int(response.headers['x-rate-limit-retry-after-seconds']))
                    print("Your rate limit is now restored")

                    response = requests.get("https://api.disgenet.com/api/v1/gda/summary",\
                                            params=params, headers=HTTPheadersDict, verify=False)
                    if response.ok is True:
                        break
                    else:
                        continue

        response_parsed = json.loads(response.text)

        if not page_count:
            total_disease_count = response_parsed["paging"]["totalElements"]
        
        disease_list = [d['diseaseName'] for d in response_parsed["payload"]]
        full_disease_list.extend(disease_list)

        page_count += 1

    return set(full_disease_list)


