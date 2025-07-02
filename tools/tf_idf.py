from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics import pairwise_distances
import numpy as np


def fps_tfidf(sentences, k, max_features=5000):
    """
    Farthest‐Point Sampling on TF–IDF (or raw counts).
    
    Args:
      sentences   (List[str]): your corpus
      k           (int):       how many to sample
      max_features(int):       max vocab size (for speed)
    
    Returns:
      List[int]: indices of the selected sentences
    """
    # 1) Vectorize
    vec = TfidfVectorizer(max_features=max_features, stop_words='english')
    X = np.array(vec.fit_transform(sentences).todense())
    
    # 2) Precompute cosine distances
    dist = pairwise_distances(X, metric='cosine')
    
    n = len(sentences)
    if k >= n:
        return list(range(n))
    
    # 3) Seed: pick the point farthest from the centroid
    centroid = X.mean(axis=0)
    centroid_dists = pairwise_distances(X, centroid[None, :], metric='cosine').reshape(-1)
    first = int(np.argmax(centroid_dists))
    selected = [first]
    
    min_dist = dist[first, :].copy()
    
    for _ in range(1, k):
        # Add the next best sentence to the selected list
        next_idx = int(np.argmax(min_dist))
        selected.append(next_idx)

        # Minimum here because we pick the least similar sentence to the selected sentences
        min_dist = np.minimum(min_dist, dist[next_idx, :])
    
    return selected


if __name__ == "__main__":
    sentences = [
        "The cat sat on the mat.",
        "A quick brown fox jumps over the lazy dog.",
        "Machine learning models can learn embeddings.",
        "Sentence diversity is useful for summarization.",
        "Natural language processing includes tokenization, parsing, and more.",
        "The lazy dog was too slow in the race."
    ]
    picks = fps_tfidf(sentences, k=3, use_tfidf=True, max_features=1000)
    print("Selected sentences:")
    for idx in picks:
        print(f" - [{idx}] {sentences[idx]}")
