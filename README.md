# Missense Classification Using AdaBoost
This repository contains datasets and trained models of my Masters thesis work on - Supervised Classification of Missense Mutations as Pathogenic or Tolerated using Ensemble Learning Methods.

datasets/ folder contains 
- compiled Training dataset that was used to train our models
- specially compiled Benchmark dataset that contains minimal overlap (2-15%) with the 13 existing state-of-the-art classifiers (SIFT, PANTHER, PROVEAN, PhD-SNP, Mutation Assessor, FATHMM, SNPs&GO, SNPs&GO^3D, nsSNPAnalyzer, PolyPhen2, PON-P2, MutPred, SNAP, CONDEL and MetaSNP.
- Mutations in Seen and Unseen proteins of Benchmark dataset in our Training dataset
- Mutations in pure and mixed proteins (containing only 1 type of mutation class vs. containing both) of our Benchmark dataset.

models/ folder contains
- model trained using long range amino acid dependency preservation only
- model trained using 16 other sequence conservation scores
- model trained using sequence conservation and 21 physico-chemical properties.

evaluation/ folder contains graphs indicating performance metrics of the 3 models on the different datasets.
