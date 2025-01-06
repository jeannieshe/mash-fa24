## Transcriptomics-based histological scoring for metabolic-associated steatohepatitis using machine learning.

### About
A collection of notes and scripts from my fall 24 UROP in Lauffenburger Lab (MIT Department of Biological Engineering) on Metabolic Dysfunction-Associated Steatohepatitis (MASH). Thank you to Nikos Meimetis and Professor Doug Lauffenburger for the continued mentorship and guidance throughout this research.

### Project
MASH, the advanced stage of metabolic-associated steatosis liver disease, is characterized by severe accumulation of fat in the liver. Approximately 25% of adults in the US have metabolic dysfunction-associated steatotic liver disease (MASLD); of those, 20% have MASH [(Cleveland Clinic)](https://my.clevelandclinic.org/health/diseases/22988-nonalcoholic-steatohepatitis). However, we currently observe a lack of clinical therapies directly targeting liver inflammation [(Alexopoulos et al)](https://pubmed.ncbi.nlm.nih.gov/38030391/). Patients often receive anti-obesity or anti-diabetes treatments instead.

A patient's disease severity is determined by a doctor's histological scoring of a liver biopsy, which produces two clinical scores:
- Fibrosis stage score (discrete scale from 0 to 4, with 4 being most fibrotic)
- NAS (NASH Activity Score) (discrete scale from 0 to 8, with 8 being most severe)

The effort towards developing accurate *in vitro* liver models is crucial for better understanding disease progression and therapy development.

**My goal for this semester was to create a machine learning model to histologically score transcriptomic data with fibrosis stage score and NAS, bringing clinical meaning to in vitro liver model conditions.**

*Please [click here](NEET_Poster_Dec_2024.pdf) for a complete project poster that was presented in December 2024 at the MIT NEET Living Machines poster presentation.*

### Skills gained
Through the 14-week semester, I challenged myself to train and save machine learning models from scratch. 
- I started with a template of PyTorch models from my mentor, but performed benchmarking of these models from a blank document to practice the ML workflow on my own. 
- I trained models on the RNA-sequencing data from *in vivo* samples from [Pantano et al.](https://www.nature.com/articles/s41598-021-96966-5), [Govaere et al.](https://pubmed.ncbi.nlm.nih.gov/33268509/), and [Hoang et al.](https://pubmed.ncbi.nlm.nih.gov/31467298/). I later tested models on *in vitro* samples from [Kostrzewski et al.](https://doi.org/10.1038/s42003-021-02616-x).
- I familiarized myself with the cross-validation workflow for KNN, PLSR, SVM, and ElasticNet models. 
- I learned to compare these models using the Pearson's R score. 
- While PLSR was the simplest or one of the least computationally expensive models, it performed ultimately the best and was the one I selected to build my final model.
- I perused documentation to identify the best methods from packages that I should use (selecting sklearn's StratifiedKFold method for cross-validation, for example).
- I learned to save models with .pkl files.
- I spent many hours debugging my models until my models were reproducible.
- I became comfortable producing files to switch between analysis in Python and R.
- I learned to use ChatGPT to my advantage, especially to learn code to plot my data for presentable use.

### Script specifics
- In 0 and 1, I am utilizing GSEA and VIPER packages from R to try and predict an alternative score for *in vivo* samples. The models may potentially better predict these alternative scores rather than the NAS and fibrosis clinical scores.
- In 2, I am training my PLSR models.
- In 3, I am plotting the VIPER calculated scores vs. ground truth.
- In 4, 5, 6, 7, and 8, I train many different kinds of models to evaluate which might be best performing.
- In 9, I choose the PLSR model trained on the Govaere et al. data to move forward.
- In 10, 11, and 12, I predict NAS and fibrosis scores for the *in vitro* data and plot them against each other as a two-way verification.