Supplementary material (code for result reproducibility) for article:
## Visual Analytics Unveils a Two-miRNA Signature for High-Accuracy Discrimination Between Pancreatic Adenocarcinoma and Well-Differentiated Pancreatic Neuroendocrine Tumors

Submitted to:
### Journal of Gastroenterology

### Authors
José M Enguita (1), Ignacio Díaz (1), Tamara Cubiella (2,3,4), Diego García (1), Nuria Valdés (5) and María D. Chiara (2,3,4)

(1)	Department of Electrical Engineering, University of Oviedo, Gijón 33204, Spain 
(2)	Health Research Institute of the Principado de Asturias (ISPA), Hospital Universitario Central de Asturias, Oviedo 33011, Spain
(3)	Instituto Universitario de Oncología del Principado de Asturias (IUOPA), University of Oviedo, Spain
(4)	CIBERONC (Network of Biomedical Research in Cancer), Madrid 28029, Spain
(5)	Hospital Universitario Cruces, Bilbao, Bizkaia. Biobizkaia, CIBERER, CIBERDEM, EndoERN 

Corresponding authors
N. Valdés nvalga@gmail.com and M.D. Chiara mdchiara.uo@uniovi.es 

### Dependencies
The packages needed to reproduce the results of the article are listed in `requirements.txt`. Run `pip install -r requirements.txt` to install all packages using pip.

### Instructions

run `download-gse.py` and `download-tcga.py` (in any order) to download and preprocess data.

run `python -i make_plots.py` to create figures and get the results

