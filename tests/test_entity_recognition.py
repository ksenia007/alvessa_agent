import pytest
from unittest.mock import patch, MagicMock

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import entity_extraction
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
        'chr_pos_variants': {}
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
        }
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
        }
    }
    
    result = entity_extraction.variant_extraction_node(initial_state)
    assert result == expected_output

def test_variant_extraction_node_no_variants():
    """Tests behavior when no variants are present in the text."""
    input_text = "This text discusses genes like BRCA1 but has no genetic variants."
    initial_state = {"messages": [{"content": input_text}]}
    
    expected_output = {'variants': {}, 'chr_pos_variants': {}}
    
    result = entity_extraction.variant_extraction_node(initial_state)
    assert result == expected_output

# --- Test Post-Processing Logic ---

def test_post_processing():
    """Tests the _post_process_entities helper function."""
    raw_genes = ["BRCA1, TP53", "EGFR ", " short ", "a", "miRNA", "a_very_long_gene_name_that_is_invalid"]
    raw_traits = ["cancer, disease", " short ", "hf"]

    expected = {
        "genes": ["BRCA1", "TP53", "EGFR", "short"],
        "traits": ["cancer", "disease", "short"]
    }
    
    result = entity_extraction._post_process_entities(raw_genes, raw_traits)
    # Using sorted() is a standard way to compare lists when order doesn't matter
    assert sorted(result["genes"]) == sorted(expected["genes"])
    assert sorted(result["traits"]) == sorted(expected["traits"])

# --- Test Individual Model Nodes (with Mocks) ---

@patch('entity_extraction.claude_call')
def test_claude_entity_extraction_node(mock_claude_call):
    """Tests the Claude node by mocking the claude_call function."""
    mock_response = MagicMock()
    mock_response.content = [MagicMock()]
    mock_response.content[0].text = "BRCA1, TP53, EGFR"
    mock_claude_call.return_value = mock_response

    initial_state = {"messages": [{"content": "Some text for Claude."}]}
    result = entity_extraction.claude_entity_extraction_node(initial_state)
    
    assert sorted(result["genes"]) == sorted(["BRCA1", "TP53", "EGFR"])

@patch('entity_extraction._get_gliner_model')
def test_gliner_entity_extraction_node(mock_get_gliner):
    """Tests the GLiNER node by mocking the model loader."""
    mock_model = MagicMock()
    mock_model.predict_entities.return_value = [
        {'text': 'BRCA1', 'label': 'gene'},
        {'text': 'breast cancer', 'label': 'disease'},
        {'text': 'TP53', 'label': 'protein'}
    ]
    mock_get_gliner.return_value = mock_model
    
    initial_state = {"messages": [{"content": "Some text for GLiNER."}]}
    result = entity_extraction.gliner_entity_extraction_node(initial_state)

    assert sorted(result["genes"]) == sorted(["BRCA1", "TP53"])
    assert sorted(result["traits"]) == sorted(["breast cancer"])

@patch('entity_extraction._get_flair_model')
def test_flair_entity_extraction_node(mock_get_flair):
    """Tests the Flair node by mocking the model loader."""
    mock_model = MagicMock()
    
    def mock_predict(sentence):
        mock_label_gene = MagicMock()
        mock_label_gene.data_point.text = "CFTR"
        mock_label_gene.value = "Gene"
        
        mock_label_trait = MagicMock()
        mock_label_trait.data_point.text = "cystic fibrosis"
        mock_label_trait.value = "Disease"
        
        sentence.get_labels.return_value = [mock_label_gene, mock_label_trait]

    mock_model.predict.side_effect = mock_predict
    mock_get_flair.return_value = mock_model
    
    with patch('entity_extraction.Sentence') as mock_sentence_class:
        mock_sentence_instance = MagicMock()
        mock_sentence_class.return_value = mock_sentence_instance

        initial_state = {"messages": [{"content": "Some text for Flair."}]}
        result = entity_extraction.flair_entity_extraction_node(initial_state)
        
        assert sorted(result["genes"]) == sorted(["CFTR"])
        assert sorted(result["traits"]) == sorted(["cystic fibrosis"])

# --- Test Merged Node ---

@patch('entity_extraction._extract_entities_with_gliner')
@patch('entity_extraction._extract_entities_with_flair')
@patch('entity_extraction._extract_entities_with_claude')
def test_entity_extraction_node_merged(mock_claude, mock_flair, mock_gliner):
    """Tests the main merged node, ensuring results are combined and deduplicated."""
    mock_claude.return_value = {"genes": ["BRCA1", "TP53"], "traits": []}
    mock_flair.return_value = {"genes": ["TP53"], "traits": ["hereditary breast cancer", "cancer"]}
    mock_gliner.return_value = {"genes": ["BRCA2"], "traits": ["cancer"]}

    input_text = "Discussion on BRCA1, TP53, BRCA2 and cancer risk, mentioning rs12345 and chr1:123:A>G."
    initial_state = {"messages": [{"content": input_text}]}

    result = entity_extraction.entity_extraction_node(initial_state)
    
    expected_variants = {'rs12345': {'rsid': 'rs12345'}}
    expected_chr_pos = {
        'chr1:123:A>G': {'coordinates': {'chrom': '1', 'pos': 123, 'ref': 'A', 'alt': 'G', 'assembly': 'GRCh38.p14'}}
    }

    assert sorted(result['genes']) == sorted(['BRCA1', 'TP53', 'BRCA2'])
    assert sorted(result['traits']) == sorted(['hereditary breast cancer', 'cancer'])
    assert result['variants'] == expected_variants
    assert result['chr_pos_variants'] == expected_chr_pos

# --- Test Edge Condition Helpers using @pytest.mark.parametrize ---

# Define states once to be reused in the parameterized test
state_with_all = {"genes": ["BRCA1"], "traits": ["cancer"], "variants": {"rs123": {}}}
state_with_genes_only = {"genes": ["BRCA1"], "traits": [], "variants": {}}
state_with_nothing = {"genes": [], "traits": [], "variants": {}}

@pytest.mark.parametrize("func, state, expected_result", [
    # Test cases for has_genes
    (entity_extraction.has_genes, state_with_all, True),
    (entity_extraction.has_genes, state_with_genes_only, True),
    (entity_extraction.has_genes, state_with_nothing, False),
    # Test cases for has_traits
    (entity_extraction.has_traits, state_with_all, True),
    (entity_extraction.has_traits, state_with_genes_only, False),
    (entity_extraction.has_traits, state_with_nothing, False),
    # Test cases for has_variants
    (entity_extraction.has_variants, state_with_all, True),
    (entity_extraction.has_variants, state_with_genes_only, False),
    (entity_extraction.has_variants, state_with_nothing, False),
])
def test_has_entities_conditions(func, state, expected_result):
    """Tests the edge condition helpers (has_genes, has_traits, has_variants)."""
    assert func(state) is expected_result
