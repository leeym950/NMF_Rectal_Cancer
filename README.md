# NMF_Rectal_Cancer

This file explains about this repository.

1. Find a new classification for rectal cancer(not including colon cancer), by using NMF method.

2. Analyze if the classification is correlated with clinical data. (pCR, survival, etc)

3. Validate the classification system, using Korean rectal cancer data, from Severance Hospital, Rep. of Korea.

4. Only R source codes are provided by this repository.
   Data are not included, due to large size > 25MB.
   
# Results

1. Found NMF rank = 4 is most promisible, but rank = 2 and rank =3 is also analyzed.

2. When rank = 2 or rank =3, survival analysis showed no significance.

3. When rank = 4, one group showed superior survival relatively to other three groups.

4. Geneset templates were made using PAMR package, distinguishing previously shown 'distinguish' group from other three.

5. Validation with Korean rectal cancer data showed no significant results.

Summary> Negative results!
