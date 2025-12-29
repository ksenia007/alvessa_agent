import pytest
from unittest.mock import patch, MagicMock

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.alvessa.agents import entity_extraction
# --- Test Variant Extraction Logic ---

def test_variant_extraction_node_rsids():
    """Tests the variant_extraction_node for various rsID formats."""
    input_text = "The study found associations with rs123, RS456, and rs123 (a duplicate)."
    initial_state = {"messages": [{"content": input_text}]}
    
    expected_output = {
        'variants': {
            'rs123': {'rsid': 'rs123'},
            'rs456': {'rsid': 'rs456'}
        },
        'chr_pos_variants': {},
        "ensembl_genes": [],
        "ensembl_proteins": [],
        "ensembl_transcripts": []
    }
    
    result = entity_extraction.variant_extraction_node(initial_state)
    assert sorted(result['variants'].keys()) == sorted(expected_output['variants'].keys())
    assert result['chr_pos_variants'] == expected_output['chr_pos_variants']

def test_variant_extraction_node_chr_pos():
    """Tests the variant_extraction_node for chr:pos:ref>alt format."""
    input_text = "Variants include chr7:55249071:C>T and chrX:1000:A>G. Also a duplicate chr7:55249071:C>T."
    initial_state = {"messages": [{"content": input_text}]}
    
    expected_output = {
        'variants': {},
        'chr_pos_variants': {
            'chr7:55249071:C>T': {
                'coordinates': {'chrom': '7', 'pos': 55249071, 'ref': 'C', 'alt': 'T', 'assembly': 'GRCh38.p14'}
            },
            'chrX:1000:A>G': {
                'coordinates': {'chrom': 'X', 'pos': 1000, 'ref': 'A', 'alt': 'G', 'assembly': 'GRCh38.p14'}
            }
        },
        "ensembl_genes": [],
        "ensembl_proteins": [],
        "ensembl_transcripts": []
    }
    
    result = entity_extraction.variant_extraction_node(initial_state)
    assert result == expected_output

def test_variant_extraction_node_mixed_and_invalid():
    """Tests for a mix of valid variants and invalid formats."""
    input_text = "Check rs789 and chr1:12345:A>T, but ignore chr1:bad:A>T and rsABC."
    initial_state = {"messages": [{"content": input_text}]}
    
    expected_output = {
        'variants': {
            'rs789': {'rsid': 'rs789'}
        },
        'chr_pos_variants': {
            'chr1:12345:A>T': {
                'coordinates': {'chrom': '1', 'pos': 12345, 'ref': 'A', 'alt': 'T', 'assembly': 'GRCh38.p14'}
            }
        },
        "ensembl_genes": [],
        "ensembl_proteins": [],
        "ensembl_transcripts": []
    }
    
    result = entity_extraction.variant_extraction_node(initial_state)
    assert result == expected_output

def test_variant_extraction_node_no_variants():
    """Tests behavior when no variants are present in the text."""
    input_text = "This text discusses genes like BRCA1 but has no genetic variants."
    initial_state = {"messages": [{"content": input_text}]}
    
    expected_output = {'variants': {}, 'chr_pos_variants': {}, "ensembl_genes": [], "ensembl_proteins": [], "ensembl_transcripts": []}
    
    result = entity_extraction.variant_extraction_node(initial_state)
    assert result == expected_output

def test_variant_extraction_node_ensembl_ids():
    """Tests the variant_extraction_node for various Ensembl ID formats."""
    input_text = "Gene ENSG00000157764, protein ENSP00000288602, transcript ENST00000288602. Duplicate ensg00000157764 and another gene ENSG00000243485."
    initial_state = {"messages": [{"content": input_text}]}
    
    expected_output = {
        'variants': {},
        'chr_pos_variants': {},
        'ensembl_genes': ['ENSG00000157764', 'ENSG00000243485'],
        'ensembl_proteins': ['ENSP00000288602'],
        'ensembl_transcripts': ['ENST00000288602']
    }
    
    result = entity_extraction.variant_extraction_node(initial_state)

    # Compare dictionaries and sorted lists to ensure order-agnostic comparison
    assert result['variants'] == expected_output['variants']
    assert result['chr_pos_variants'] == expected_output['chr_pos_variants']
    assert sorted(result['ensembl_genes']) == sorted(expected_output['ensembl_genes'])
    assert sorted(result['ensembl_proteins']) == sorted(expected_output['ensembl_proteins'])
    assert sorted(result['ensembl_transcripts']) == sorted(expected_output['ensembl_transcripts'])

# --- Test Post-Processing Logic ---

def test_post_processing():
    """Tests the _post_process_entities helper function."""
    raw_genes = ["BRCA1, TP53", "EGFR ", " short ", "a", "miRNA", "a_very_long_gene_name_that_is_invalid", "ENSG00000157764"]
    raw_traits = ["cancer, disease", " short ", "hf"]

    expected = {
        "genes": ["BRCA1", "TP53", "EGFR", "short", "ENSG00000157764"],
        "traits": ["cancer", "disease", "short"]
    }
    
    result = entity_extraction._post_process_entities(raw_genes, raw_traits)
    assert sorted(result["genes"]) == sorted(expected["genes"])
    assert sorted(result["traits"]) == sorted(expected["traits"])

# --- Test Individual Model Nodes (with Mocks) ---

@patch('src.alvessa.agents.entity_extraction.claude_call')
def test_claude_entity_extraction_node(mock_claude_call):
    """Tests the Claude node by mocking the claude_call function."""
    mock_response = MagicMock()
    mock_response.content = [MagicMock()]
    mock_response.content[0].text = "BRCA1, TP53, EGFR"
    mock_claude_call.return_value = mock_response

    initial_state = {"messages": [{"content": "Some text for Claude."}]}
    result = entity_extraction.claude_entity_extraction_node(initial_state)
    
    assert sorted(result["genes"]) == sorted(["BRCA1", "TP53", "EGFR"])

@patch('src.alvessa.agents.entity_extraction._get_gliner_model')
def test_gliner_entity_extraction_node(mock_get_gliner):
    """Tests the GLiNER node, checking for separate gene/protein extraction."""
    mock_model = MagicMock()
    mock_model.predict_entities.return_value = [
        {'text': 'BRCA1', 'label': 'gene'},
        {'text': 'breast cancer', 'label': 'disease'},
        {'text': 'TP53', 'label': 'protein'}
    ]
    mock_get_gliner.return_value = mock_model
    
    initial_state = {"messages": [{"content": "Some text for GLiNER."}]}
    result = entity_extraction.gliner_entity_extraction_node(initial_state)

    assert result["genes"] == ["BRCA1"]
    assert result["proteins"] == ["TP53"]
    assert result["traits"] == ["breast cancer"]

@patch('src.alvessa.agents.entity_extraction._get_flair_model')
def test_flair_entity_extraction_node(mock_get_flair):
    """Tests the Flair node, checking for separate gene/protein extraction."""
    mock_model = MagicMock()
    
    def mock_predict(sentence):
        # Mock labels for gene, protein, and disease
        mock_label_gene = MagicMock(value="Gene", data_point=MagicMock(text="CFTR"))
        mock_label_protein = MagicMock(value="Protein", data_point=MagicMock(text="AKT1"))
        mock_label_trait = MagicMock(value="Disease", data_point=MagicMock(text="cystic fibrosis"))
        
        sentence.get_labels.return_value = [mock_label_gene, mock_label_protein, mock_label_trait]

    mock_model.predict.side_effect = mock_predict
    mock_get_flair.return_value = mock_model
    
    with patch('src.alvessa.agents.entity_extraction.Sentence') as mock_sentence_class:
        mock_sentence_instance = MagicMock()
        mock_sentence_class.return_value = mock_sentence_instance

        initial_state = {"messages": [{"content": "Some text for Flair."}]}
        result = entity_extraction.flair_entity_extraction_node(initial_state)
        
        assert result["genes"] == ["CFTR"]
        assert result["proteins"] == ["AKT1"]
        assert result["traits"] == ["cystic fibrosis"]


@patch('src.alvessa.agents.entity_extraction.claude_call')
def test_verify_genes_in_query_parses_json(mock_claude):
    """Verifier should honor JSON lists from Claude and map back to canonical casing."""
    mock_response = MagicMock()
    mock_chunk = MagicMock()
    mock_chunk.text = '["BRCA1", "TP53"]'
    mock_response.content = [mock_chunk]
    mock_claude.return_value = mock_response

    query = "We studied BRCA1 and TP53 expression patterns."
    genes = ["BRCA1", "TP53", "EGFR"]

    result = entity_extraction._verify_genes_in_query(query, genes)
    assert result == ["BRCA1", "TP53"]


@patch('src.alvessa.agents.entity_extraction.claude_call', side_effect=RuntimeError("LLM offline"))
def test_verify_genes_in_query_fallback(mock_claude):
    """Verifier should fall back to text matching if Claude fails."""
    query = "Only TP53 was mentioned explicitly."
    genes = ["BRCA1", "TP53"]

    result = entity_extraction._verify_genes_in_query(query, genes)
    assert result == ["TP53"]

# --- Test Merged Node ---

@patch('src.alvessa.agents.entity_extraction._verify_genes_in_query')
@patch('src.alvessa.agents.entity_extraction._extract_entities_with_gliner')
@patch('src.alvessa.agents.entity_extraction._extract_entities_with_flair')
@patch('src.alvessa.agents.entity_extraction._extract_entities_with_claude')
def test_entity_extraction_node_merged(mock_claude, mock_flair, mock_gliner, mock_verify):
    """Tests the main merged node, ensuring results are combined and deduplicated correctly."""
    mock_claude.return_value = {"genes": ["BRCA1", "TP53"], "traits": []}
    mock_flair.return_value = {"genes": ["TP53"], "proteins": ["AKT1"], "traits": ["hereditary breast cancer", "cancer"]}
    mock_gliner.return_value = {"genes": ["BRCA2"], "proteins": ["AKT1", "P53"], "traits": ["cancer"]}
    mock_verify.return_value = ["BRCA1", "TP53", "BRCA2", "ENSG00000157764"]

    input_text = "Discussion on genes, proteins, cancer risk, rs12345, chr1:123:A>G, and ENSG00000157764."
    initial_state = {"messages": [{"content": input_text}]}

    result = entity_extraction.entity_extraction_node(initial_state)
    
    expected_variants = {'rs12345': {'rsid': 'rs12345'}}
    expected_chr_pos = {
        'chr1:123:A>G': {'coordinates': {'chrom': '1', 'pos': 123, 'ref': 'A', 'alt': 'G', 'assembly': 'GRCh38.p14'}}
    }

    assert sorted(result['genes']) == sorted(['BRCA1', 'TP53', 'BRCA2', 'ENSG00000157764'])
    assert sorted(result['traits']) == sorted(['hereditary breast cancer', 'cancer'])
    assert sorted(result['proteins']) == sorted(['AKT1', 'P53'])
    assert result['transcripts'] == [] # No transcripts in the input text
    assert result['variants'] == expected_variants
    assert result['chr_pos_variants'] == expected_chr_pos


@patch('src.alvessa.agents.entity_extraction._verify_genes_in_query')
@patch('src.alvessa.agents.entity_extraction._extract_entities_with_gliner')
@patch('src.alvessa.agents.entity_extraction._extract_entities_with_flair')
@patch('src.alvessa.agents.entity_extraction._extract_entities_with_claude')
def test_entity_extraction_detects_drugs(mock_claude, mock_flair, mock_gliner, mock_verify):
    """Ensures small-molecule parsing captures product names, synonyms, and metadata."""
    mock_claude.return_value = {"genes": [], "traits": []}
    mock_flair.return_value = {"genes": [], "proteins": [], "traits": []}
    mock_gliner.return_value = {"genes": [], "proteins": [], "traits": []}
    mock_verify.return_value = []

    input_text = "Vonoprazan (also called TAK-438 free base) was combined with ML162 during screening."
    initial_state = {"messages": [{"content": input_text}]}

    result = entity_extraction.entity_extraction_node(initial_state)

    assert {"Vonoprazan", "ML162"}.issubset(set(result["drugs"]))

    von_key = entity_extraction._canon_drug_key("Vonoprazan")
    ml_key = entity_extraction._canon_drug_key("ML162")

    assert von_key in result["drug_entities"]
    assert ml_key in result["drug_entities"]

    von_entry = result["drug_entities"][von_key]
    mention_terms = {m["text"] for m in von_entry["mentions"]}
    assert any("TAK-438" in term for term in mention_terms)
    assert von_entry["catalog_number"] == "HY-100007"
    assert "TAK-438 (free base)" in von_entry["synonyms"]

# --- Test Edge Condition Helpers using @pytest.mark.parametrize ---

state_with_all = {"genes": ["BRCA1"], "traits": ["cancer"], "variants": {"rs123": {}}, "proteins": ["P53"], "drugs": ["Vonoprazan"], "transcripts": ["ENST123"]}
state_with_genes_only = {"genes": ["BRCA1"], "traits": [], "variants": {}, "proteins": [], "drugs": [], "transcripts": []}
state_with_nothing = {"genes": [], "traits": [], "variants": {}, "proteins": [], "drugs": [], "transcripts": []}

@pytest.mark.parametrize("func, state, expected_result", [
    (entity_extraction.has_genes, state_with_all, True),
    (entity_extraction.has_genes, state_with_genes_only, True),
    (entity_extraction.has_genes, state_with_nothing, False),
    
    (entity_extraction.has_traits, state_with_all, True),
    (entity_extraction.has_traits, state_with_genes_only, False),
    (entity_extraction.has_traits, state_with_nothing, False),

    (entity_extraction.has_variants, state_with_all, True),
    (entity_extraction.has_variants, state_with_genes_only, False),
    (entity_extraction.has_variants, state_with_nothing, False),
    
    (entity_extraction.has_proteins, state_with_all, True),
    (entity_extraction.has_proteins, state_with_genes_only, False),
    (entity_extraction.has_proteins, state_with_nothing, False),

    (entity_extraction.has_transcripts, state_with_all, True),
    (entity_extraction.has_transcripts, state_with_genes_only, False),
    (entity_extraction.has_transcripts, state_with_nothing, False),

    (entity_extraction.has_drugs, state_with_all, True),
    (entity_extraction.has_drugs, state_with_genes_only, False),
    (entity_extraction.has_drugs, state_with_nothing, False),
])
def test_has_entities_conditions(func, state, expected_result):
    """Tests the edge condition helpers for all entity types."""
    assert func(state) is expected_result
