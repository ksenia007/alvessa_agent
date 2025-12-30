# Adding a New Tool to Alvessa

This tutorial demonstrates how to add a new tool to Alvessa.

## Overview

Tools in Alvessa are modular components that:
- Fetch data from external APIs or local databases
- Annotate entities (genes, variants, etc.) with additional information
- Add structured data and text summaries to entities
- Can have dependencies on other tools

## File Structure

Each tool lives in `src/tools/<tool_name>/` with:
```
tool_name/
├── node.py          # Main tool implementation (required)
├── query.py         # API/database query functions (optional)
└── __init__.py      # Package initialization (optional)
```

## Anatomy of a Tool

### 1. Header Docstring

Every tool starts with a header :

```python
"""
Description:

Tool to query [Database/API] for [what it does].
"""
```

### 2. Imports

Import dependencies + config:

```python
from src.state import State          # Alvessa internal
from src.tools.base import Node      # Alvessa internal
# + any required imports

DEBUG = True  # Enable debug logging
```

### 3. Main Tool Function

The core function that performs the annotation:

```python
def tool_name_agent(state: "State"):
    """
    Annotate [entities] with [what data].

    Parameters
    ----------
    state
        Current graph state.

    Returns
    -------
    None
        Updates entities in place.
    """
```

**Key Points:**
- Takes `state: "State"` as the first parameter
- Can take additional optional parameters

### 4. Tool Logic

The function follows this pattern:

For variants: 

```python
def tool_agent(state: "State"):
    # 1. Get entities from state
    variant_entities = state.get("variant_entities") or {}

    # 2. Iterate through entities
    for variant_name, variant in variant_entities.items():
        if not variant_name:
            continue

        try:
            # 3. Fetch/process data
            # ... API calls, database queries, etc.
            # use variant.get_location(build) to get variant locations

            # 4. Update entity with results
            variant.update_text_summaries(summary)
            variant.add_tool("tool_name")

        except Exception as e:
            if DEBUG:
                print(f"[ToolName] Error processing {variant_name}: {e}")

    # 5. Log completion
    if DEBUG:
        print(f"[ToolName] Processing complete")

    return
```

And for genes:

```python
def tool_agent(state: "State"):
    # 1. Get entities from state
    gene_entities = state.get("gene_entities") or {}

    # 2. Iterate through entities
    for gene_name, gene in gene_entities.items():
        if not gene_name:
            continue

        try:
            # 3. Fetch/process data
            # ... API calls, database queries, etc.
            # use gene.symbol to get the gene symbol

            # 4. Update entity with results
            gene.update_text_summaries(summary)
            gene.add_tool("tool_name")

        except Exception as e:
            if DEBUG:
                print(f"[ToolName] Error processing {gene_name}: {e}")

    # 5. Log completion
    if DEBUG:
        print(f"[ToolName] Processing complete")

    return
``` 


### 5. Adding information from the tools

Tools use these methods to update entities:

**For Variants:**
```python
variant.update_text_summaries(text)      # Add text summary
variant.add_tool(tool_name)              # Record tool was run
```

**For Genes:**
```python
gene.update_text_summaries(text)         # Add text summary
gene.add_tool(tool_name)                 # Record tool was run
```

### 7. Node Registration

Register the tool at the module level:

```python
NODES: tuple[Node, ...] = (
    Node(
        name="tool_name",                   # Tool identifier
        entry_point=tool_name_agent,        # Function to call
        description=(                        # LLM-facing description
            "Clear description of what this tool does."
        ),
        dependencies=("dependency_tool",),  # Optional: tools that must run first
    ),
)
```

**Node Parameters:**
- `name`: Unique identifier for the tool
- `entry_point`: The function to execute
- `description`: Description used by the LLM for tool selection

## Best Practices

### Error Handling

Use try-except blocks within the loop to prevent one failure from stopping all processing:

```python
for entity_name, entity in entities.items():
    try:
        # Process entity
        pass
    except Exception as e:
        if DEBUG:
            print(f"[Tool] Error processing {entity_name}: {e}")
```

### Rate Limiting

Add delays when calling external APIs:

```python
response = requests.get(url, params=params)
time.sleep(0.3)  # courteous pause
```

### Text Summaries

Format summaries consistently with a tool prefix:

```python
summary = f"*ToolName: {description of what was found}."
entity.update_text_summaries(summary)
```

### Tool Naming

Use consistent, descriptive names:
- Function: `tool_name_agent` (e.g., `dbsnp_variants_agent`)
- Node name: `tool_name` (e.g., `"variant_annotations"`)
- Tool tracking: Match node name (e.g., `add_tool("variant_annotations")`)


### Debug Logging

Use consistent debug logging:

```python
if DEBUG:
    print(f"[ToolName] Processing {len(entities)} entities")
    print(f"[ToolName] Error: {error_message}")
    print(f"[ToolName] Complete")
```

## Testing Your Tool

After implementing your tool we recommend to test with a simple query containing relevant entities

## Common Pitfalls

- ❌ Don't modify state structure unnecessarily
- ❌ Don't return state
- ❌ Don't use bare `except:` - be specific or use `except Exception:`
- ❌ Don't skip the header docstring
- ❌ Don't forget to add the tool name with `entity.add_tool()`
- ✅ Do iterate through entities, not modify the dict during iteration
- ✅ Do add text summaries with the `*ToolName:` prefix
- ✅ Do handle missing data gracefully with early `continue`
- ✅ Do use DEBUG logging for development and troubleshooting
- ✅ Do use detailed tool descriptions
