Supplementary material (code for result reproducibility) for article:
## Visual Analytics Unveils a Two-miRNA Signature for High-Accuracy Discrimination Between Pancreatic Adenocarcinoma and Well-Differentiated Pancreatic Neuroendocrine Tumours

### Authors
Nuria Valdés, MD-PhD (6); José M Enguita, PhD (1); Ignacio Díaz, PhD (1); Tamara Cubiella, Msc (2,3,4); Diego García, PhD (1); Raúl Rodriguez, MD-PhD (2,5); María Poch, MD (2,5) and María D. Chiara, MD-PhD (2,3,4)
1.	Department of Electrical Engineering, University of Oviedo, Gijón 33204, Spain 
2.	Health Research Institute of the Principado de Asturias (ISPA), Hospital Universitario Central de Asturias, Oviedo 33011, Spain
3.	Instituto Universitario de Oncología del Principado de Asturias (IUOPA), University of Oviedo, Spain
4.	CIBERONC (Network of Biomedical Research in Cancer), Madrid 28029, Spain
5.	Department of Pathology, Hospital Universitario de Cabueñes, Gijón 33394, Spain
6.	Department of Endocrinology and Nutrition, Hospital Universitario Cruces, Bilbao, Bizkaia. Biobizkaia, CIBERER, CIBERDEM, EndoERN 


Corresponding authors
J.M. Enguita jmenguita@uniovi.es and M.D. Chiara mdchiara.uo@uniovi.es 

### Dependencies
The packages needed to reproduce the results of the article are listed in `requirements.txt`. 

Run `pip install -r requirements.txt` to install all packages using pip.

### Instructions

run `download-gse.py` and `download-tcga.py` (in any order) to download and preprocess data.

run `python -i make_plots.py` to create figures and get the results

