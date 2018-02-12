# rSVD-single-pass
Randomized SVD with single pass over data matrix.

This package includes the codes and experimental data regarding paper: 

Wenjian Yu, Yu Gu, Jian Li, Shenghua Liu, and Yaohang Li, Single-Pass PCA of Large High-Dimensional Data,
to appear in Proc. IJCAI 2017, pp. 3350-3356. Arxiv.org report <a href="http://arxiv.org/abs/1704.07669"> #1704.07669</a>.

<b>If you used this package, please cite the above paper, or send email to me letting me know. Thank you!</b>

1.This package includes the codes in Matlab for generating the test data, and drawing the figures in the paper, and the codes in C implementing the algorithms tested in the paper.

2.The Matlab programs for generating test matrices.

-PCAtestmatrix.m: generate 5 types of test matrices.

-genLargeMatrix.m: a script generating large matrix stored on hard disk.

-genFeretMatrix.m: generate the 150GB matrix from the FERET database.

(To obtain the FERET database, please follow the instructions on https://www.nist.gov/itl/iad/image-group/color-feret-database)

3.The Matlab programs for the algorithms.

-rSVDbasic.m: the basic randomized algorithm for computing SVD. (Algorithm 1)

-rSVD_exSP.m: the existing single-pass algorithm for randomized SVD. (Algorithm 2)

-rSVDsp.m: the proposed single-pass algorithm for computing SVD. (Algorithm 4)

4.The Matlab programs for drawing figures.

-showSigmaDist.m: a script drawing the singular value distribution (Fig. 2).

-plotfig3.m: draw singular values obtained from the algorithms for two matrices (Fig. 3)

-acc_valid.m: detailed comparison of the obtained singular values (for Fig. 4).

-acc_singvec.m: detailed comparison of the obained singular vectors (for Fig. 5).

-drawFeretSigma: draw the singular value of Feret matrix (Fig. 6a).

-drawEigface.m: draw the four eigenfaces (Fig. 6b).

5.The C program for generating large test data and implementing the algorithms.

Please unpack the tar file rSVDsp-C.tar, and then read the README there.
The program has been tested on a Ubuntu Linux machine, with Intel ICC with MKL library.

For comment/question/suggestion, please send email to yu-wj at tsinghua dot edu dot cn (Dr. Wenjian Yu).

 
