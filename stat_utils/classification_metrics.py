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
    labels: np.ndarray | list,
    scores: np.ndarray | list,
    score_threshold: float = 0.5,
    random_state: int = 42,
) -> dict:
    """
    Compute binary classification metrics given true labels and predicted scores.
    
    Parameters
    ----------
    labels : array-like
        True binary labels (0 or 1).
    scores : array-like
        Predicted scores or probabilities.
    score_threshold : float
        Threshold to convert predicted scores into binary predictions (default: 0.5).
    random_state : int
        Random seed for reproducibility (default: 42).
        
    Returns
    -------
    dict
        Dictionary containing computed metrics:
        - "auroc": Area Under the Receiver Operating Characteristic Curve.
        - "auprc": Area Under the Precision-Recall Curve.
        - "rand_auroc": AUROC for random predictions.
        - "rand_auprc": AUPRC for random predictions.
        - "accuracy": Accuracy of predictions.
        - "precision": Precision of predictions.
        - "early_precision": Precision for the top 10,000 predictions.
        - "recall": Recall of predictions.
        - "f1": F1 score of predictions.
        - "n_edges": Total number of samples.
        - "n_pos": Number of positive samples.
        - "n_neg": Number of negative samples.
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