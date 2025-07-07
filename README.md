# Alvessa: agentic AI scientist

*v0 of the submission*

To run sample requests:
```
# install environment from requirements.txt
export ANTHROPIC_API_KEY=[..your key here..]
export BioGRID_API_KEY=[..your key here..]
python -m run
```

The structure of the code is as follows:

- run.py – main script to run the graph and save outputs.
- graph_builder.py – builds the LangGraph pipeline and saves the png
- claude_client.py – handles Claude API calls with rate-limiting - should only catch extreme cases though
- config.py – stores constants like model names and flags. Also where API key is pulled in
- state.py – defines the shared LangGraph state structure (what is stored in the state)

API - based tools:
- tool_humanbase.py – fetches and filters HumanBase predictions
- tool_uniprot.py – retrieves UniProt entries, extracts traits and also has additional processing

Graph nodes:
- entity_extraction.py – extracts gene symbols from user input using the model defines in config (GENE_EXTRACT_MODEL)
- conditioned_claude.py – builds context and sends main Claude query. This is the main reasoning agent right now.
- verify.py – checks if Claude’s answer matches the evidence

Output files:
- demo.txt – saved answers and evidence
- demo.log – print/debug log
- graph_diagram.png – LangGraph visualization (optional)