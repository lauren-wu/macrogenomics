# macrogenomics
The python codes used in the study of macrogenomics engineering
## Codes

1. mcExpressionOut.py: Output the result from Monte Carlo simulation and Brownian Dynamics simulation.
	
2. macrogenomics.py: The class required to do any chromatin packing macromolecular crowding model calculation
	
3. CP_MC_Model.py: Calculate the model predicted gene expression sensitivity from chromatin packing macromolecular crowding model and the sensitivity measured from microarray experiment.

	Requirement: Need the macrogenomics.py, full_genes.xlsx and the output from mcExpressionOut.py 
	
4. percentVariance.py: Calculate the variance expression by the model in experimental data.
	
	Requirement: Need the macrogenomics.py, full_genes.xlsx and the output from mcExpressionOut.py 
	
5. heterogeneity.py: Calculate the model predicted intercellular gene expression heterogeneity.
	
	Requirement: Need the macrogenomics.py, full_genes.xlsx and the output from mcExpressionOut.py 
	
6. cov.py: Calculate the model predicted intercellular coefficient of expression variation.
	
	Requirement: Need the macrogenomics.py, full_genes.xlsx and the output from mcExpressionOut.py 
	


## Input data:

full_genes.xlsx: the microarray measurement result
	
There are 4 groups of microarray measurements in total. Each group contains 4 replicated measurements. In the full_genes.xlsx file, the first columns correspond to the control measurement and 5-8, 9-12 and 13-16 columns correspond to the treated groups. The relative sigma values for them are: 1.001, 1, 0.9933 and 0.9151 


## Output data:
	
1. mcExpressionOut.py: max_mRNA_initial.csv, phi_initial.csv, second_derivative_TF_norm.csv, tot_con.csv
	
2. CP_MC_Model.py: sensitivity_model.csv, sensitivity_model_g.csv sensitivity_experiment.csv, g_function.csv
	
3. percentVariance.py: percentOfVariance.csv
	
4. heterogeneity.py: heterogeneity_model.csv, heterogeneity_experiment.csv
	
5. cov.py: cov_model.csv, cov_experiment.csv
