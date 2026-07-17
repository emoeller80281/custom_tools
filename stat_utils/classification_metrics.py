import numpy as np
from sklearn.metrics import (
    f1_score,
    roc_auc_score,
    average_precision_score,
    accuracy_score,
    precision_score,
    recall_score
)

def compute_binary_classification_metrics(
    labels,
    scores,
    score_threshold: float = 0.5,
    random_state: int = 42,
):
    """
    Compute binary classification metrics given true labels and predicted scores.
    
    labels: array-like of 0/1 labels
    scores: array-like of predicted probabilities after sigmoid
    """

    labels = np.asarray(labels).astype(int).ravel()
    scores = np.asarray(scores).astype(float).ravel()

    preds = (scores >= score_threshold).astype(int)

    accuracy = accuracy_score(labels, preds)
    precision = precision_score(labels, preds, zero_division=0)
    recall = recall_score(labels, preds, zero_division=0)
    f1 = f1_score(labels, preds, zero_division=0)
    
    preds_sorted_indices = np.argsort(scores)[::-1]
    preds_sorted = preds[preds_sorted_indices]
    labels_sorted = labels[preds_sorted_indices]
    
    early_precision = precision_score(labels_sorted[:10_000], preds_sorted[:10_000], zero_division=0)

    if len(np.unique(labels)) < 2:
        auroc = np.nan
        auprc = np.nan
        rand_auroc = np.nan
        rand_auprc = np.nan
    else:
        auroc = roc_auc_score(labels, scores)
        auprc = average_precision_score(labels, scores)

        rng = np.random.default_rng(random_state)
        rand_scores = rng.permutation(scores)

        rand_auroc = roc_auc_score(labels, rand_scores)
        rand_auprc = average_precision_score(labels, rand_scores)

    return {
        "auroc": auroc,
        "auprc": auprc,
        "rand_auroc": rand_auroc,
        "rand_auprc": rand_auprc,
        "accuracy": accuracy,
        "precision": precision,
        "early_precision": early_precision,
        "recall": recall,
        "f1": f1,
        "n_edges": len(labels),
        "n_pos": int(labels.sum()),
        "n_neg": int((labels == 0).sum()),
        "score_threshold": score_threshold,
    }