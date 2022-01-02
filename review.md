Reviewer #1 (Comments to the Author ):

The manuscript by Wu T et al. describes a method for the quantification of neoantigen-mediated negative selection in cancer evolution.

It would be interesting to identify the presence of clusters of neoantigens with different ability to regulate the negative selection and escape in cancer evolution. Some antigens may have immunostimulatory function whereas other may have immunosuppressive function. The authors should provide this analysis.

In figure 5, the authors should consider adding an analysis where elimination and escape (yes and no) are evaluated by using the ratio CTL CD8+/Treg and NK/Treg cells within the same tumor sample rather then the absolute count of CTL CD8+ and NK cells alone.


Reviewer #2 (Comments to the Author ):

The work by Wu and colleagues addresses an important and controversial aspect of cancer evolution, the extent of immune mediated negative selection. The authors developed a quantitative metric for detecting two phases of immunoediting: elimination and escape. Specifically, they estimate two enrichment scores, ESccf based on the allele frequency distribution of neoantigenic and non-neoantigenic mutations, and ESrna based on the observed expression of neoantigen and non-neoantigenic mutations. They hypothesize that ESccf represents a metric of elimination, while ESrna represents a metric of escape. By permuting the neoantigen labels, the authors create a background distribution that can be used to estimate a P-value. The authors then apply this methodology to TCGA datasets and determine the proportion of tumors under elimination and escape. They further determine that the signals of elimination and escape are anticorrelated based on the median ES. Next, the authors used the model developed by Lakatos et al to demonstrate the association between selective coefficients and their ESccf score. The authors then attempt to correlate metrics of immune infiltration with those tumor with a strong elimination signal, although they report a slight increase of CD8 T cells and NK cells in those tumors with immune selection, they found this value not significant. Finally, the authors leveraged three datasets with immunotherapy information and found a significant association of elimination with increased survival.

While the idea of developing a non-parametric unsupervised statistic to derive an enrichment score for the allele frequency distribution and the expression of neoantigenic and non-neoantigenic mutations is very interesting, the description of the methodology is poor and lacks scientific detail making it hard to assess and several parts of the method section seems to be copy-pasted from other referenced manuscripts. It is also unclear what metrics/parameters were used in all analyses, i.e when they used driver mutations and the RECIST classification.

The novelty of their results comes as incremental evidence in favor of negative selection during cancer evolution. However, their strategy led to conclusions that have already been described: negative selection depends on the strength of selective coefficients (Lakatos et al 2020) and the emergence of escape is a confounding factor for negative selection (Lakatos et al 2020, Zapata et al 2020).

OVerall, the manuscript is well written, the strategy for detecting selection is sound, but on top of the issues that I have mentioned earlier, the depth of their analysis falls short for such an important debate in the field of cancer evolution and immunogenomics.

Some major specific concerns:

1. Figure 1 describes the conceptual framework of their stud-y. During the elimination phase, CCF of neoantigenic mutations will be lower than the rest of mutations as these are under strong selective pressure. I do believe in this idea that has been explored previously, however, it is hard to convey that escape can only be measured by the downregulation of expression of neoantigens. As the authors state, there are multiple mechanisms of immune evasion, but they are totally neglected in all their analysis. For example, they could easily obtain PD1 or PDL1 expression and determine their association with the elimination or escape signal.

2. The authors claim that mutation types do not influence the CCF distribution but I found it hard to believe that missense or truncating mutations, or just mutations in oncogenic and tumor suppressor genes, have effectively the same CCF distributions. The authors do not provide a comparison between them. They only show the aggregated CCF distribution for all cancer types and all mutations. In addition, the IC50 considered in this study considers the lowest possible value for a peptide with the mutation, but it is unclear if it is the lowest among all possible 9-mer surrounding the mutation or something different.

3. By using ESccf the authors validate the existence of an immunoediting-elimination signal in the cancer genome. However, this is only observed as significant in 2/30 (when including driver genes) or 4/30 (when excluding driver genes) which represent a mild signal overall and can reflect a high proportion of immune escape across multiple tumor types. By using ESrna, the authors claim a stronger immunoediting-escape signal than elimination, which indeed seems to be the case as 16/30 tumors show a significant difference. The authors then claim that both signals are anticorrelated indicating that escape masks the strength of immuno-editing elimination in tumors. It is then important that the authors perform all their analysis by excluding samples with a strong immunoediting escape signal for assessing immunoediting elimination. How different would the results be when performing this exclusion?

4. The authors claim that the minimum to quantify the immunoediting signal is one neoantigen. This seems very unlikely, as 1) you can not compute a distribution based signal with one value, 2) one neoantigen can be easily mispredicted as shown by studies that attempt to validate neoantigens in vitro, the success rate is very poor, 3) the neoantigen derived ES signal for elimination or escape can be overlapping, what happen when both scenarios are co-occurring low ESccf and low ESrna?. Also, How are these measures (both ESccf and ESrna) affected by the number of mutations / neoantigens present in a tumour? One would think the computation is not particularly reliable for low numbers?

6. How did they (if they did) account for the CCF>1 cluster (Lines 153)? Mutations in regions with high copy number can give rise to strange CCFs, which could lead to CCF under or overestimate, depending on if the mutation is on major or minor allele (e.g. CN=3:1, mutation is present in 66% but on the n=3 chromosome, it will be over-estimated).

7. How is this CCF distribution affected by sequencing noise? E.g. is it affected by the resolution of sequencing (that partly shapes the tail)?
   Also, why put more weight on the two tails?

9. The authors only used an immune infiltration classification from the ImmuneCellAI study. A proper comparison with other immune cell quantification methods should be conducted (Danaher et al, Thorsson et al, CIBERSORT, etc)

10. Lines 532. I don't quite get this logic, ESrna means the sample is downregulating neoantigens via RNA downreg, but Treg would be immune escape, which could actually mean no downregulation is needed. Also the same with other mechanisms of escape. If an antigen presenting associated gene is mutated (B2M), downregulation of neoantigenic expression will not be under selection anymore. Thus evidence of low ESrna does not mean escape.

11. Also, it is puzzling that tumors with the strongest ESccf would have the best response to immunotherapy as these patients will by definition have the lowest amount of neoantigens available due to the elimination process of immunoediting. Is it possible that by only accounting for the observed neoantigens the signal is actually reversed, and it reflects the efficacy of the immune cells on clearing neoantigens.

12. The authors claim that methods that use the ratio of dnds "inside" to "outside" HLA-binding regions is not a good metric to estimate immune selection as the majority of neoantigens will emerge in "outside" regions. However, it is my understanding that the idea of having an inside versus outside dnds estimation relies on the power given by the large portion of the genome "outside" that is under other selective forces, or even under neutral evolution (passenger). Therefore, despite the amount of neoantigens, these outside regions will have a tendency to be neutral (dn/ds~1), except in cases where selection is extremely strong or is acting on all the exome (which would mean that every possible change is immunogenic, which is very unlikely). You can similarly test the same using the ESccf of mutation HLA binding regions versus non-binding regions and compare.

Minor

1. In lines 128-132, the authors copy paste the description of Bailey et al for predicting driver mutations, they should paraphrase the section and explain why they choose a cut-off of more than one approach for their estimate. What happens if they use one or 3 approaches as a cutoff. Does this distinction change the results? Specially because in Bailey et al they described that CTAT cancer score outperforms the rest of categories.

2. Also, specific driver mutations will be enriched in oncogenes as the work from Bayleys suggest. The authors should separate missense from truncating or frameshifts events.

3. In some sections, It is unclear if the authors performed the calculation of the KS statistic between all mutations, or only mutations classified as drivers.

4. Immunogenic mutations are not the same as antigenic mutations, or even mutations that are predicted to bind MHC complexes. Authors should rigorously separate this distinction and clarify their terminology throughout all their manuscript, For example in lines 183 they refer to "immunogenic", and in lines 204 they refer to "neoantigen" mutations.

5. Supplementary Fig S2 conveys the CCF distribution of all TCGA mutations. However, it would be better to show the CCF distribution of all tumor types separated, plotting as well the CCF of neoantigens and non-neoantigens to understand the distribution.

6. One of their key metrics is ESrna, the authors mention that z denotes mRNA expression. In the manuscript is TPM but this is not described in the methods.

7. In their equations (should be numbered), the cumulative distributions of neo and noneo is based on subscript j, which is not described.Same for subscript i. I had to guess that i is the interval,

8. Several typos/orthographic errors: Lines 374, "A mutation was considered neoantigenic if there was at least one peptide derived from the mutated sequence is predicted as neoantigen ', the verb "is' should be removed ' (This is only one example)

9. Authors should number the equations.

10. Lines 376 This observation is not novel, please specify that this was already observed in Lakatos et al. Also place the references to the model where they correspond.

11. It is unclear whether the neoantigen prediction is done using the personalised 6-HLA alleles of each individual, a cohort HLA, a proto-HLA.

12. How robust is to use the peptide with the strongest binding affinity. What if the authors instead use an score that leverages all possible peptides such as the one used in MArty et al (2017) PMID: 29107334


Reviewer #3 (Comments to the Author ):

The current paper seeks to demonstrate a novel method for demonstrating neoantigen mediated negative selection in cancer evolution. This method is unique in that it does seem to demonstrate the existence of neoantigen mediated negative selection which could further provide a translationally relevant output if this method is able to be clinically validated. The methodology used herein appear appropriate based on other previously published literature and the authors use a well-established database for this endeavor, allowing them to achieve an appropriate sample size. Overall, this manuscript seeks to provide some important information; however, there are some areas the authors might consider addressing, as follows:

Major Points:

1. Though the introduction seems to present an unbiased take on the previously published results, the experimental rationale is vague and no clear hypothesis is stated.
2. There is an over-interpretation of the data for the biomarker use of this technique. The authors make the statement (starting at line 570) that the Hosmer-Lemeshow test used to compare ESCCF with TMB and Neoantigen Burden suggests that the ESCCF is more suitable for predicting prognosis than the other tests based on the p-value of the H-L test; however, this is not the case. While a significant p-value resulting from an H-L test would indicate a lack of model fit, the comparison of p-values of this test, do not conversely indicate "goodness of fit" such that one could say that one test is superior to another based on the p-value obtained. Additionally, the statement was made that the p-value for TMB is "close to significant"; however, I would caution the authors on making such a statement because the data based on use of appropriate statistics either is or is not significant, and therefore, the p-value references in line 570, based on the level of significance set, is not significant. Interpreting this value as near-significant thus represents an overinterpretation of the data. Though the Kaplan-Meier curve shows differences in survival probability between high and low ESCCF values which could suggest that this is a novel predictive biomarker for ICI therapy, that this is a better test than others remains unclear. Additionally, multivariate regression or other analysis looking at such factors which may have otherwise played a role in survival should have been assessed, such as ICI type (Nivolumab vs Pembrolizumab), patient gender, etc, to assess if there are other factors which may have been inherently different between these groups when interpreting these survival results.
3. While this study does use a large and well-studied pan-tumor database, to create its new quantification system, in vitro or in vivo would have been helpful in validating these in silico studies.
4. There should be consideration in moving supplemental figure 9 to the main body of the text. There are 12 lines in the text (502-514) spent as a pre-amble/explanation of the significance of this figure, and if it deserves such explanation, this suggests the text cannot stand alone without this and it should be potentially be moved.
5. There are many grammatical challenges with this paper which in some parts hinder the flow of the document. I would encourage the authors to seek proof-reading and revisions to this manuscript to ensure appropriate grammar with particular attention to subject-verb agreement, pluralization and sentence fragmentation.


Minor Points:

1. Not all of the data is described in the text of the results. There is no mention of the data from Figure 4B within the body of the text.
2. There are minor labeling inconsistencies in the text and the figures. For example, on line 352 (Fig 1a, b) should be denoted (Fig 1A, B) for consistency.
3. There are areas of redundancy within the text that could be streamlined for readability. For example: Line 581 and 598 are essentially the same.
4. There are areas of data interpretation in the results section which should be transferred to the discussion where it is more appropriate (ex, line 452, clause starting with "probably").
5. The sentence starting at line 473 is unclear and confusing.
6. In line 612 "n/s" is used, but not previously defined, is this mean to be "dN/dS"?