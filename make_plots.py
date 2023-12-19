
import pandas as pd
import numpy as np
from scipy import stats


def print_auc_ci(labels,scores,title):    
    '''
        Uses the Fast DeLong algorithm to estimate the variance of the AUC
        and estimates the 95% confidence interval (CI) taking +/- 1.96 
        standard deviation

        Original de DeLong paper

            DeLong, E. R., DeLong, D. M., & Clarke-Pearson, D. L. (1988). 
            Comparing the areas under two or more correlated receiver operating 
            characteristic curves: a nonparametric approach. Biometrics, 837-845.

        Fast DeLong algorithmm paper 

            Sun, Xu, and Weichao Xu. "Fast implementation of DeLong’s algorithm 
            for comparing the areas under correlated receiver operating characteristic 
            curves." IEEE Signal Processing Letters 21.11 (2014): 1389-1393.

        Library implementing the Fast Delong algorithm
        
            https://github.com/yandexdataschool/roc_comparison

    '''
    import compare_auc_delong_xu

    # estimation of variance using Fast DeLong method
    auc, var = compare_auc_delong_xu.delong_roc_variance(labels, scores)
    
    # lower and upper bounds
    low = auc - np.sqrt(var)*1.96   # lower bound for 95% confidence level
    sup = auc + np.sqrt(var)*1.96   # upper bound for 95% confidence level

    # clip the bounds to the [0,1] interval
    low = np.clip(low,0,1)
    sup = np.clip(sup,0,1)

    # print the results
    print('')
    print('----------------------------------')
    print(title + ':')
    print(f'AUCvar (DeLong): {var}')
    print(f'AUC = {auc} - CI = {low,sup}')
    print('----------------------------------')
    print('')




# Samples from the TCGA

keep_primary_only = True
remove_pseudo = True
remove_extra = True

print ('loading TCGA data ...')
df = pd.read_hdf("datasets/df_TCGA.hdf")

print ('preparing TCGA data ...')
df = df[df.cancer=='PAAD']

# separate sane pancreatic tissue
df_TCGA_Normal = df[df.type == 'sane']

# Keep only miRNA data
df_TCGA_Normal = df_TCGA_Normal.iloc[:,df_TCGA_Normal.columns.str.contains('MIMAT')]

# Keep only primary tumors
if(keep_primary_only):
	df=df[df.type=='primary']

# Remove pseudotumors (PERAN et al) if wanted
if (remove_pseudo):
     pseudo=["TCGA-F2-6880", "TCGA-HZ-7924", "TCGA-F2-7273", "TCGA-IB-AAUV","TCGA-F2-7276", 
             "TCGA-IB-AAUW", "TCGA-H8-A6C1", "TCGA-RL-AAAS", "TCGA-HZ-7920", "TCGA-US-A77J", "TCGA-HZ-7923",
             "TCGA-L1-A7W4", "TCGA-FB-A7DR",  "TCGA-HZ-8638"," TCGA-HZ-7289", "TCGA-HZ-7289", "TCGA-HZ-7918-01","TCGA-FB-AAPP",
             "TCGA-HV-A7OP", "TCGA-2J-AABP", "TCGA-IB-7654"]
     df = df[~df.index.isin(pseudo)]


# Remove some samples that were not present in the Morphing Projections app after data curation
if (remove_extra):
     extra=["TCGA-Q3-A5QY","TCGA-FB-A4P6","TCGA-3A-A9IL"]
     df = df[~df.index.isin(extra)]

# Keep only miRNA data
df = df.iloc[:,df.columns.str.contains('MIMAT')]

# Separate PANNETS
pannet_list = ['TCGA-3A-A9IN','TCGA-2L-AAQM','TCGA-3A-A9IS','TCGA-3A-A9IV','TCGA-3A-A9IJ','TCGA-3A-A9IR','TCGA-3A-A9IO']  
df_TCGA_PANNET = df[df.index.isin(pannet_list)]
df_TCGA_PDAC = df[~df.index.isin(pannet_list)]

# Normalize miRNA 129-5p and 203a data
print ('normalizing TCGA data ...')
list_TCGA_PANNET_129=[stats.percentileofscore(df.values.ravel(),x) for x in df_TCGA_PANNET['MIMAT0000242']]
list_TCGA_PDAC_129=[stats.percentileofscore(df.values.ravel(),x) for x in df_TCGA_PDAC['MIMAT0000242']]
list_TCGA_Normal_129=[stats.percentileofscore(df.values.ravel(),x) for x in df_TCGA_Normal['MIMAT0000242']]

list_TCGA_PANNET_203=[stats.percentileofscore(df.values.ravel(),x) for x in df_TCGA_PANNET['MIMAT0000264']]
list_TCGA_PDAC_203=[stats.percentileofscore(df.values.ravel(),x) for x in df_TCGA_PDAC['MIMAT0000264']]
list_TCGA_Normal_203=[stats.percentileofscore(df.values.ravel(),x) for x in df_TCGA_Normal['MIMAT0000264']]


# Read database GSE163031
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163031
print ('loading GSE163031 data ...')
df_GSE163031 = pd.read_hdf('datasets/df_GSE163031.hdf')

print('cleaning and normalizing ...')
lista_mimats = list(np.load('mimats.npy'))
mm = [x for x in df_GSE163031.columns if x in lista_mimats]
df_GSE163031 = df_GSE163031[mm + ['Type']]

# Get a dataframe to the miRNA values for normalization
df=df_GSE163031.loc[:,df_GSE163031.columns != 'Type']

# Normalize miRNA 129-5p and 203a data
df_GSE163031['MIMAT0000242']=[stats.percentileofscore(df.values.ravel(),x) for x in df_GSE163031['MIMAT0000242']]
list_GSE163031_Normal_129 = list(df_GSE163031[df_GSE163031.Type=='Normal']['MIMAT0000242'])
list_GSE163031_PDAC_129 = list(df_GSE163031[df_GSE163031.Type=='Pancre']['MIMAT0000242'])

df_GSE163031['MIMAT0000264']=[stats.percentileofscore(df.values.ravel(),x) for x in df_GSE163031['MIMAT0000264']]
list_GSE163031_Normal_203 = list(df_GSE163031[df_GSE163031.Type=='Normal']['MIMAT0000264'])
list_GSE163031_PDAC_203 = list(df_GSE163031[df_GSE163031.Type=='Pancre']['MIMAT0000264'])

# Read database GSE73367
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73367
print ('loading GSE73367 data ...')
df_GSE73367 = pd.read_hdf('datasets/df_GSE73367.hdf')

print('cleaning and normalizing ...')
# Get only miRNA hsa
df_GSE73367= df_GSE73367.iloc[:,df_GSE73367.columns.str.contains('hsa-')]

# Normalize miRNA 129-5p and 203a data
list_GSE73367_129=[stats.percentileofscore(df_GSE73367.values.ravel(),x) for x in df_GSE73367['hsa-miR-129-5p']]
list_GSE73367_203=[stats.percentileofscore(df_GSE73367.values.ravel(),x) for x in df_GSE73367['hsa-miR-203']]

# Read database GSE43796
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43796
print ('loading GSE43796 data ...')
df_GSE_43796 = pd.read_hdf('datasets/df_GSE43796.hdf')

print('cleaning and normalizing ...')
df=df_GSE_43796.loc[:,df_GSE_43796.columns != 'Type']


print('preparing miRNA data ...')
# Normalize miRNA 129-5p and 203a data and separate Type types
df_GSE_43796['A_25_P00013880']=[stats.percentileofscore(df.values.ravel(),x) for x in df_GSE_43796['A_25_P00013880']]
df_GSE_43796_Normal_129 = df_GSE_43796[df_GSE_43796.Type=='Normal']['A_25_P00013880']
df_GSE_43796_PDAC_129 = df_GSE_43796[df_GSE_43796.Type=='PDAC']['A_25_P00013880']
df_GSE_43796_PANNET_129 = df_GSE_43796[df_GSE_43796.Type=='PNET']['A_25_P00013880']
list_GSE_43796_PDAC_129 = list(df_GSE_43796_PDAC_129.values)
list_GSE_43796_Normal_129 = list(df_GSE_43796_Normal_129.values)
list_GSE_43796_PANNET_129 = list(df_GSE_43796_PANNET_129.values)


df_GSE_43796['A_25_P00010628']=[stats.percentileofscore(df.values.ravel(),x) for x in df_GSE_43796['A_25_P00010628']]
df_GSE_43796_Normal_203 = df_GSE_43796[df_GSE_43796.Type=='Normal']['A_25_P00010628']
df_GSE_43796_PDAC_203 = df_GSE_43796[df_GSE_43796.Type=='PDAC']['A_25_P00010628']
df_GSE_43796_PANNET_203 = df_GSE_43796[df_GSE_43796.Type=='PNET']['A_25_P00010628']
list_GSE_43796_PDAC_203 = list(df_GSE_43796_PDAC_203.values)
list_GSE_43796_Normal_203 = list(df_GSE_43796_Normal_203.values)
list_GSE_43796_PANNET_203 = list(df_GSE_43796_PANNET_203.values)


# Draw the viloin diagrams
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 24})
plt.close('all')
plt.ion()

import matplotlib.patches as mpatches
labels = []
def add_label(violin, label):
    color = violin["bodies"][0].get_facecolor().flatten()
    labels.append((mpatches.Patch(color=color), label))

def miviolin(valores, posiciones, etiqueta, color=None):
    vp = plt.violinplot(valores, positions=posiciones, showmedians=True, widths=0.45)
    plt.boxplot(valores, positions=posiciones, widths=0.10)                    

    if color==None:
        color = "red"
    for pc in vp['bodies']:
        pc.set_facecolor(color)
        pc.set_edgecolor('black')
        pc.set_linewidth(1)
        pc.set_alpha(0.8)
    
    for pc in ('cbars','cmins','cmaxes','cmedians'):
        p = vp[pc]
        p.set_edgecolor("black")
        p.set_linewidth(1)
  
    
    add_label(vp,etiqueta)


rrred = '#cc2222'
bluuu = '#555599'
greeeen = '#22cc22'

lnet =  list_TCGA_PANNET_129 + list_GSE_43796_PANNET_129 + list_GSE73367_129
lpdac = list_TCGA_PDAC_129 + list_GSE163031_PDAC_129 + list_GSE_43796_PDAC_129
lnorm = list_TCGA_Normal_129 + list_GSE163031_Normal_129 + list_GSE_43796_Normal_129

lnet2 =  list_TCGA_PANNET_203 + list_GSE_43796_PANNET_203 + list_GSE73367_203
lpdac2 = list_TCGA_PDAC_203 + list_GSE163031_PDAC_203 + list_GSE_43796_PDAC_203
lnorm2 = list_TCGA_Normal_203 + list_GSE163031_Normal_203 + list_GSE_43796_Normal_203

nnorm   = len(lnorm)
npannet = len(lnet)
npdac   = len(lpdac)

pos=[1,2,3,5,6,7]
labels2=[f"PANNET\n(n={npannet})",f"PDAC\n(n={npdac})",f"Normal\n(n={nnorm})",
         f"PANNET\n(n={npannet})",f"PDAC\n(n={npdac})",f"Normal\n(n={nnorm})"]
plt.figure()
labels = []
miviolin([lnet],[pos[0]],"PANNET", color=rrred)
miviolin([lpdac],[pos[1]],"PDAC", color=bluuu) 
miviolin([lnorm],[pos[2]],"Normal", color=greeeen) 
miviolin([lnet2],[pos[3]],"PANNET", color=rrred)
miviolin([lpdac2],[pos[4]],"PDAC", color=bluuu) 
miviolin([lnorm2],[pos[5]],"Normal", color=greeeen) 



plt.ylabel('miRNA percentile normalized expression level')
plt.xticks([x for x in pos], labels2)


# Add multiple comparisons p-value for mean difference -----------
# Plot lines indicating what means are compared
# 'tick_len' gives the length of the tick on the end of each line

tick_len = 0.25*10
hh = 7
plt.plot([1, 1, 2, 2], [105 - tick_len, 105, 105, 105 - tick_len], c="black")
plt.plot([1, 1, 3, 3], [105 + hh - tick_len, 105 + hh, 105 + hh, 105 + hh - tick_len], c="black")


sp = 0.02
plt.plot([5, 5, 6 - sp, 6 - sp], [105 - tick_len, 105, 105, 105 - tick_len], c="black")
plt.plot([6 + sp, 6 + sp, 7, 7], [105 - tick_len, 105, 105, 105 - tick_len], c="black")

def asterisk_p_string(value):
     if value <= 0.0001:
          return r"****"
     elif value <= 0.001:
          return r"***"
     elif value <= 0.01:
          return r"**"
     elif value <= 0.05:
          return r"*"
     else: return r"ns"

r_net_pdac_129 = stats.f_oneway(lnet,lpdac)
print(f'pNet -> pDAC miR-129-5p p-value: {r_net_pdac_129.pvalue}')

r_net_pdac_203 = stats.f_oneway(lnet2,lpdac2)
print(f'pNet -> pDAC miR-203a p-value: {r_net_pdac_203.pvalue}')

r_norm_net_129 = stats.f_oneway(lnorm,lnet)
print(f'normal -> pNET miR-129-5p p-value: {r_norm_net_129.pvalue}')

r_norm_pdac_203 = stats.f_oneway(lnorm2,lpdac2)
print(f'normal -> pDAC miR-203a p-value: {r_norm_pdac_203.pvalue}')


pad = 0.2
plt.text(1.5, 105 + pad, asterisk_p_string(r_net_pdac_129.pvalue), fontsize=18, va="bottom", ha="center")
plt.text(2, 105 + hh + pad, asterisk_p_string(r_norm_net_129.pvalue), fontsize=18, va="bottom", ha="center")
plt.text(5.5, 105 + pad, asterisk_p_string(r_net_pdac_203.pvalue), fontsize=18, va="bottom", ha="center")
plt.text(6.5, 105 + pad, asterisk_p_string(r_norm_pdac_203.pvalue), fontsize=18, va="bottom", ha="center")
plt.ylim([0, 119.8])

plt.text(2,120,'miR-129-5p',fontsize=20, va="bottom", ha="center")
plt.text(6,120,'miR-203a',fontsize=20, va="bottom", ha="center")

plt.show()


#####################################################
# ROC CURVE

import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc


# Assuming you have two sets of data: actual labels and predicted scores/probabilities
# actual_labels: true labels (0 or 1)
# predicted_scores: predicted scores or probabilities for class 1
# Compute the false positive rate, true positive rate, and threshold using roc_curve

actual_labels=np.array([1]*len(lnet) + [0]*len(lpdac))

# Start with miRNA-129-5p
print("miRNA 129-5p:")
predicted_scores = lnet + lpdac

fpr, tpr, thresholds = roc_curve(actual_labels, predicted_scores)

# print the AUC estimate with 95% confidence intervals
print_auc_ci(actual_labels,np.array(predicted_scores),'miRNA 129-5p')


# Compute the area under the ROC curve (AUC)
roc_auc = auc(fpr, tpr)

# Plot the ROC curve

plt.rcParams.update({'font.size': 24})

plt.ion()
plt.figure(2)
plt.plot(fpr, tpr, label='ROC curve miRNA-129-5p (AUC = %0.2f)' % roc_auc, linewidth=2.5, marker = 'o')

plt.figure()
plt.plot(thresholds, tpr,label='TPR'); plt.plot(thresholds,fpr,label='FPR');
plt.legend();
plt.xlabel('Threshold (percentile normalized expression level)')
plt.ylabel('Rate')
plt.title('True Positive Rate (TPR) and False Positive Rate (FPR) curves')

th_optimal = thresholds[np.argmax(tpr - fpr)]
print(f'Optimal threshold: {th_optimal}')

from sklearn.metrics import recall_score, precision_score, f1_score
sensitivity = recall_score(actual_labels , predicted_scores>th_optimal)
specificity = recall_score(np.logical_not(actual_labels) , np.logical_not(predicted_scores>th_optimal))
precision  = precision_score(actual_labels, predicted_scores>th_optimal)
f1 = f1_score(actual_labels , predicted_scores>th_optimal)
print(f'Sensitivity: {sensitivity:.3f}, Specificity: {specificity:.3f}, Precision: {precision:.3f}, F1-score: {f1:.3f}')


# Now miRNA-203a
print("miRNA 203a:")
predicted_scores = lnet2 + lpdac2
fpr, tpr, thresholds = roc_curve(1-actual_labels, predicted_scores)

# print the AUC estimate with 95% confidence intervals
print_auc_ci(1-actual_labels,np.array(predicted_scores),'miRNA 203a')


# Compute the area under the ROC curve (AUC)
roc_auc = auc(fpr, tpr)

# Plot the ROC curve

plt.rcParams.update({'font.size': 24})

plt.ion()
plt.figure(2)
plt.plot(fpr, tpr, label='ROC curve miRNA-203a (AUC = %0.2f)' % roc_auc, linewidth=2.5, marker = 's')


plt.figure()
plt.plot(thresholds, tpr,label='TPR'); plt.plot(thresholds,fpr,label='FPR');
plt.legend();
plt.xlabel('Threshold (percentile normalized expression level)')
plt.ylabel('Rate')
plt.title('True Positive Rate (TPR) and False Positive Rate (FPR) curves')
plt.show()

th_optimal = thresholds[np.argmax(tpr - fpr)]
print(f'Optimal threshold: {th_optimal}')

from sklearn.metrics import recall_score, precision_score, f1_score
sensitivity = recall_score(actual_labels , predicted_scores<th_optimal)
specificity = recall_score(np.logical_not(actual_labels) , np.logical_not(predicted_scores<th_optimal))
precision  = precision_score(actual_labels, predicted_scores<th_optimal)
f1 = f1_score(actual_labels , predicted_scores<th_optimal)
print(f'Sensitivity: {sensitivity:.3f}, Specificity: {specificity:.3f}, Precision: {precision:.3f}, F1-score: {f1:.3f}')


# Now combine both

print("Combined:")
mrna1 = lnet + lpdac
mrna2 = lnet2 + lpdac2

from sklearn.linear_model import LogisticRegression

lr = LogisticRegression()

X= np.vstack([mrna1,mrna2]).T
lr.fit(X,actual_labels)

print(f'logit: {lr.coef_},{lr.intercept_}')

fpr, tpr, thresholds = roc_curve(actual_labels, lr.predict_proba(X)[:,1])


# print the AUC estimate with 95% confidence intervals
print_auc_ci(actual_labels,lr.predict_proba(X)[:,1],'COMBINATION OF TWO MIRNA')


# Compute the area under the ROC curve (AUC)
roc_auc = auc(fpr, tpr)

# Plot the ROC curve

plt.rcParams.update({'font.size': 24})

plt.ion()
plt.figure(2)
plt.plot(fpr, tpr, label='ROC curve combined (AUC = %0.2f)' % roc_auc, linewidth=2.5, marker = '^')
plt.plot([0, 1], [0, 1], 'k--')  # Diagonal line (random classifier)
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
#plt.title('Receiver Operating Characteristic')
plt.legend(loc="lower right")
plt.show()

plt.figure()
if(thresholds[0]>1):
    thresholds[0] = thresholds[0]-1

plt.plot(thresholds, tpr,label='TPR'); plt.plot(thresholds,fpr,label='FPR');
plt.legend();
plt.xlabel('Threshold (probability of being PANNET)')
plt.ylabel('Rate')
plt.title('True Positive Rate (TPR) and False Positive Rate (FPR) curves')
plt.show()

th_optimal = thresholds[np.argmax(tpr - fpr)]
print(f'Optimal threshold: {th_optimal}')


from sklearn.metrics import recall_score, precision_score, f1_score
pred = lr.predict(X)
sensitivity = recall_score(actual_labels , pred)
specificity = recall_score(np.logical_not(actual_labels) , np.logical_not(pred))
precision  = precision_score(actual_labels, pred)
f1 = f1_score(actual_labels , pred)
print(f'Sensitivity: {sensitivity:.3f}, Specificity: {specificity:.3f}, Precision: {precision:.3f}, F1-score: {f1:.3f}')


#plt.grid(True)
plt.show()