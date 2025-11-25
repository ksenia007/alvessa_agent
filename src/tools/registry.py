from __future__ import annotations

import importlib
import pkgutil
from types import ModuleType
from typing import Dict, Iterable, List, Sequence

from src.alvessa.agents.entity_extraction import entity_extraction_node

from src.tools.base import Node

EXAMPLE_TOOL_SELECTION = """EXAMPLE PIPELINES (pay attention to dependencies):

1. Variant regulatory activity (e.g. SEI):
   ["extract_entities", "query_gwas_by_gene", "variant_annotations", "sei"]

2. Variant pathogenicity (e.g. AlphaMissense):
   ["extract_entities", "query_gwas_by_gene", "variant_annotations", "alphamissense"] 

3. HumanBase Expecto variant annotation:
   ["extract_entities", "humanbase_expecto", "query_gwas_by_gene", "variant_annotations", "humanbase_tissue_expecto_annotate_variants"]

4. Gene-level functional annotation:
   ["extract_entities", "gencode_gene_node", "humanbase_functions", "uniprot_base", "reactome", "BioGRID", "Summarize_bioGRID_GO", "uniprot_gwas"]

5. Protein structure, visualization and druggability:
   ["prot"]
   
Note these are only examples, and in real life you may need to run combinations of these tools **depending on the user intent and the entities extracted**.

Respond *ONLY* with a Python list of tool names. Example: ["humanbase_functions", "uniprot_base"] or ["humanbase_functions", "uniprot_base", "query_gwas_by_gene"] or ["query_gwas_by_gene", "BioGRID"]

Do not include any explanations or extra text outside the list.
"""

PACKAGE_NAME = __name__.rsplit(".", 1)[0]
EXCLUDED_MODULES = {"agent", "shared", "registry", "base"}


def _iter_tool_modules() -> Iterable[ModuleType]:
    package = importlib.import_module(PACKAGE_NAME)
    prefix = package.__name__ + "."
    for module_info in pkgutil.iter_modules(package.__path__):
        name = module_info.name
        if name in EXCLUDED_MODULES or name.startswith("_"):
            continue
        module_name = f"{prefix}{name}.node"
        try:
            yield importlib.import_module(module_name)
        except ModuleNotFoundError:
            continue


def _extract_nodes(module: ModuleType) -> Sequence[Node]:
    nodes = getattr(module, "NODES", None)
    if nodes is None:
        node = getattr(module, "NODE", None)
        if node is None:
            return ()
        nodes = (node,)
    for node in nodes:
        if not isinstance(node, Node):
            raise TypeError(f"Module {module.__name__} exports a non-Node entry: {node!r}")
    return tuple(nodes)


def _collect_nodes() -> List[Node]:
    dynamic_nodes: List[Node] = []
    for module in _iter_tool_modules():
        dynamic_nodes.extend(_extract_nodes(module))
    return dynamic_nodes


STATIC_NODES: Sequence[Node] = (
    Node(
        name="extract_entities",
        entry_point=entity_extraction_node,
        description=(
            "Extract genes and biomedical entities from the user question using Claude + GLiNER. "
            "This must run before any downstream tool to populate gene/trait context."
        ),
    ),
)


NODES: List[Node] = list(STATIC_NODES) + _collect_nodes()


def _build_catalog_and_map(nodes: Iterable[Node]) -> tuple[Dict[str, str], Dict[str, callable]]:
    catalog: Dict[str, str] = {}
    fn_map: Dict[str, callable] = {}

    for node in nodes:
        for key, desc in node.iter_catalog_entries():
            if key in catalog:
                raise ValueError(f"Duplicate catalog entry for node key '{key}'")
            catalog[key] = desc

        for key, fn in node.iter_bindings():
            if key in fn_map:
                raise ValueError(f"Duplicate function mapping for node key '{key}'")
            fn_map[key] = fn

    return catalog, fn_map


TOOL_CATALOG, TOOL_FN_MAP = _build_catalog_and_map(NODES)
