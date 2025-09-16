# Alvessa: agentic AI scientist


To run sample requests in terminal:
```
# install environment from requirements.txt
export ANTHROPIC_API_KEY=[..your key here..]
export BioGRID_API_KEY=[..your key here..]
python -m run
```

To run as a UI-enabled session:
```
# install environment from requirements.txt
export ANTHROPIC_API_KEY=[..your key here..]
export BioGRID_API_KEY=[..your key here..]
python -m serve_ui
```

To run tests:
```
pytest tests/test_entity_recognition.py
```

The structure of the code is as follows:

- run.py – main script to run the graph and save outputs.
- src/alvessa/workflow/graph_builder.py – builds the LangGraph pipeline and saves the png
- src/alvessa/clients/claude.py – handles Claude API calls with rate-limiting - should only catch extreme cases though
- src/config.py – stores constants like model names and flags. Also where API key is pulled in
- src/state.py – defines the shared LangGraph state structure (what is stored in the state)
- src/alvessa/domain/ – domain entities (`gene_class.py`, `variant_class.py`, `gene_components.py`)

API & local db - based tools (under `src/tools/<name>/node.py`):
- humanbase – fetches and filters HumanBase predictions
- uniprot – retrieves UniProt entries, extracts traits and performs additional processing
- biogrid – wraps BioGRID interaction queries
- dbsnp – resolves variant annotations via dbSNP
- gwas – fetches GWAS trait and variant associations
- sei – runs sequence regulatory predictions
- alphamissense – annotates variants with AlphaMissense pathogenicity classes
- prot – renders AlphaFold structures with FPocket summaries

- Graph agents:
- src/alvessa/agents/entity_extraction.py – extracts gene symbols from user input using Claude + GLiNER
- src/alvessa/agents/conditioned_claude.py – builds context and sends main Claude query. This is the main reasoning agent right now.
- src/alvessa/agents/verify.py – checks if Claude’s answer matches the evidence
- src/tools/go_summarization/node.py – summarizes GO terms per gene using embeddings + UniProt fallbacks

Output files:
- demo.txt – saved answers and evidence
- demo.log – print/debug log
- demo.json - dump of the state, used in the UI-enabled version
- graph_diagram.png – LangGraph visualization (optional)

UI files:
- ui.html
- serve_ui.py
