import numpy as np
import pandas as pd

from clustering.post_procces_clusters_by_mice import compute_confusion_matrix


def test_compute_confusion_matrix():
    s1 = pd.Series([1, 1, 0], index=["a", "b", "c"])
    s2 = pd.Series([2, 3, 3], index=["a", "b", "c"])

    cm = compute_confusion_matrix(s1, s2)

    assert set(cm.index) == {2, 3}
    assert set(cm.columns) == {1, 0}

    assert cm[1][2] == 1
    assert cm[0][2] == 0
    assert cm[1][3] == 1
    assert cm[0][3] == 1


def test_compute_confusion_matrix_more_than_2_classes_eye():
    s1 = pd.Series(["0", "1", "2"], index=["a", "b", "c"])
    s2 = pd.Series(["0", "1", "2"], index=["a", "b", "c"])

    cm = compute_confusion_matrix(s1, s2)

    assert (cm.values == np.eye(3)).all()


def test_compute_confusion_matrix_more_than_2_classes():
    s1 = pd.Series(["0", "1", "2"], index=["b", "a", "c"])
    s2 = pd.Series(["0", "1", "2"], index=["a", "b", "c"])

    cm = compute_confusion_matrix(s1, s2)
    expected_result = np.eye(3)
    expected_result[[0, 1]] = expected_result[[1, 0]]

    assert (cm.values == expected_result).all()
