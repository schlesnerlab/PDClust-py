# PDClust-py
This is a highly-performant Python implementation of PDclust for clustering single cells with highly sparse CpG coverage.
The original paper introducing the algorithm can be found here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6093082/. The implementation in R was made available at https://github.com/hui-tony-zk/PDclust/tree/master.

Further information for each indiviual method is included as a Docstring comment.

This repository was created as part of a bachelor's thesis at the [Chair of Biomedical Informatics, Data Mining and Data Analytics](https://www.uni-augsburg.de/de/fakultaet/fai/informatik/prof/bioinf/mitarbeiter/matthias-schlesner/), led by Prof. Dr. Matthias Schlesner.


--- 
## Quickstart
Install the package by executing the following commands:
``` bash
git clone https://github.com/theo-felix-poeschl/PDClust-py.git
pip install PDClust-py/pdclust_py/
```
and include it with ```import pdclust_py.core as pdclust```

---
## Simple example for using the API
```python
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
import pdclust_py.core as pdclust
import seaborn as sns

path_to_bed_files = ##### include path to folder containg relevant bed files

# for use with Pandas version
dist_mat = pdclust.trivial_pdclust(path_to_bed_files)

# for use with Numba version
dist_mat = pdclust.numba_pdclust(path_to_bed_files)

square_form_dist_mat = squareform(distance_matrix)
link_mat = linkage(square_form_dist_mat, method='ward', metric='euclidean')

# Optional: Create cluster map from linkage matrix
sns.clustermap(link_mat)
plt.show()

```

#### Credits
- [@stephenkraemer](https://github.com/stephenkraemer) for the (amazing) supervision
- [@hui-tony-zk](https://github.com/hui-tony-zk) for the reference implementation in R
