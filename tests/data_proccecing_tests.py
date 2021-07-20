import numpy as np
import pandas as pd

from clustering.post_procces_clusters_by_mice import compute_confusion_matrix


def test_compute_confusion_matrix():
    s1 = pd.Series([1, 1, 0], index=["a", "b", "c"])
    s2 = pd.Series([2, 3, 3], index=["a", "b", "c"])

    cm = compute_confusion_matrix(s1, s2)

    assert set(cm.index) == {2, 3}
    assert set(cm.columns) == {1, 0}
    assert (cm.values == np.array([[0., 1.], [1., 1.]])).all()
