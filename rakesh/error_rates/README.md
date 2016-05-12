Format
==

`err_rates.delim` Error rates for various configurations of nnets (includes logistic regr, svm, their regualrized versions etc.)

`knn20a.txt`, ..., `knn20e.txt` - Accuracies of KNN (for various k) for five subsets of 20 classes.

`knn400.txt` - Accuracies of KNN (for various 20*k) with all 400 classes.

`knn_probs/` - Probabilities for KNNs (times number of neighbours) for 20 subclasses

eg:- `knn_415154_probs_15`
- 415154 is the seed indicating a choice of sub-classes
- 15 is the number of neighbours.
- This file contains a 20 x 1000 matrix with entires in [0, 15]
  - 20 is number of subclasses
  - 1000 = 20 x 50 (50 is the number of samples per class)

`nnet_probs/` - Log probabilites for various models

eg:- `logistic_noise_678328.logprobs`
- Logistic regression with input noise regularization
- This file also contains a 20 x 1000 matrix (with log probabilities)
