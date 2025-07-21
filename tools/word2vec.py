from gensim.models import Word2Vec
from sklearn.metrics import pairwise_distances
import numpy as np
from collections import defaultdict

def fps_word2vec(sentences, k, vector_size=100, window=5, min_count=1, separate_sampling=False):
    """
    Farthest‐Point Sampling on Word2Vec.
    
    Args:
      sentences   (List[str]): your corpus
      k           (int):       how many to sample
    
    Returns:
      List[int]: indices of the selected sentences
    """
    
    tokenized = [sentence.lower().split() for sentence in sentences]
    w2v_model = Word2Vec(sentences=tokenized, vector_size=vector_size, window=window, min_count=min_count, workers=1)
    def sentence_vector(words):
        valid_words = [w for w in words if w in w2v_model.wv]
        if not valid_words:
            return np.zeros(vector_size)
        return np.mean(w2v_model.wv[valid_words], axis=0)
    
    X = np.array([sentence_vector(s) for s in tokenized])
    
    def fps_on_subset(indices_subset, k_subset):
        if len(indices_subset) <= k_subset:
            return indices_subset
        
        subset_vectors = X[indices_subset]
        dist = pairwise_distances(subset_vectors, metric='cosine')
        
        centroid = subset_vectors.mean(axis=0)
        centroid_dists = pairwise_distances(subset_vectors, centroid[None, :], metric='cosine').reshape(-1)
        first = int(np.argmax(centroid_dists))
        selected_local = [first]
        
        min_dist = dist[first, :].copy()
        
        for _ in range(1, k_subset):
            next_idx = int(np.argmax(min_dist))
            selected_local.append(next_idx)
            min_dist = np.minimum(min_dist, dist[next_idx, :])
        
        return [indices_subset[i] for i in selected_local]
    
    if not separate_sampling:
        return fps_on_subset(list(range(len(sentences))), k)
    
    prefix_to_indices = defaultdict(list)
    for i, sent in enumerate(sentences):
        prefix = sent.split(':', 1)[0] 
        prefix_to_indices[prefix].append(i)
    
    sampled_indices_by_prefix = {}
    for prefix, indices_subset in prefix_to_indices.items():
        sampled_indices_by_prefix[prefix] = fps_on_subset(indices_subset, k)
    
    combined = []
    for prefix in sorted(sampled_indices_by_prefix.keys()):
        combined.extend(sampled_indices_by_prefix[prefix])
    return combined