#!/usr/bin/env python3

##*************************************************************************##
##          Step4. linear regression models of the AMR data analysis       ##          
##*************************************************************************##

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
from sklearn.linear_model import LinearRegression, LogisticRegression, Ridge, Lasso
from sklearn.model_selection import train_test_split, cross_val_score, GridSearchCV
from sklearn.metrics import accuracy_score
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier 
from sklearn.metrics import classification_report, confusion_matrix, roc_curve
from sklearn.preprocessing import scale
from sklearn.svm import SVC, LinearSVC
from sklearn.metrics import roc_curve, roc_auc_score, balanced_accuracy_score
from sklearn import metrics

import sys
if len(sys.argv)!=6: 
    print("Usage:python3 4.AMR_Learn_linear.py <feature2target.txt> <gene_location_info.txt> <antibiotics.txt> <out_dir_path> <threshold for filtering absolute coefficient>\n\n e.g.,python3 AMR_Learn_linear.py feature2target.txt Spectinomycin 0.1")
    sys.exit()

#  4.AMR_Learn_linear.py PRJNA776899.feature2target.txt GCF_000020105.1.tab PRJNA776899.antibiotics.small.txt learn 0.1

keep_log = False

feature2target = sys.argv[1]
gene_location = sys.argv[2] 
antibiotics_path = sys.argv[3]
out_path = sys.argv[4]
threshold = float(sys.argv[5])

os.makedirs(out_path, exist_ok=True) 

with open(antibiotics_path, 'r') as antibiotics_file:
    first_line = antibiotics_file.readline().strip()
    antibiotics = first_line.split('\t')
    antibiotics = antibiotics[1:]
    separator = " "
    print("antibiotics being modeled: " + separator.join(antibiotics))


# loading data
print("load data")
hsd_data = pd.read_csv(feature2target,sep='\t').fillna(0) 

# creating features and arget arrays
X = hsd_data.drop(columns=[*antibiotics, 'locus_tag'], axis=1).to_numpy()
X = scale(X)
names = hsd_data.drop(columns=[*antibiotics, 'locus_tag'], axis=1).columns

# writing out the log file
if keep_log:
    f_log = open(out_path +'/run.log', 'w')
    sys.stdout = f_log
    f_err = open(out_path +'/run.err', 'w')
    sys.stderr = f_err

train_scores = {
    'Linear Regression': [],
    'Ridge Regression': [],
    'Lasso Regression': []
}

test_scores = {
    'Linear Regression': [],
    'Ridge Regression': [],
    'Lasso Regression': []
}

cv_scores = {
    'Linear Regression': [],
    'Ridge Regression': [],
    'Lasso Regression': []
}

for antibiotic_name in antibiotics:
    print("Training model for resistance to " + antibiotic_name)

    y = hsd_data[antibiotic_name].values

    #fitting a regression model
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.3, random_state=42, shuffle=True)

    print("Supervised regression models\n")
    print("1.Linear regression model\n")
    reg = LinearRegression()
    reg.fit(X_train,y_train)
    reg_pred = reg.predict(X_test)
    reg_train_score = reg.score(X_train,y_train)
    reg_test_score = reg.score(X_test,y_test)
    reg_cv_scores = cross_val_score(reg, X_train,y_train, cv=5)
    reg_cv_score = sum(reg_cv_scores)/5
    print("Training mean accuracy score:" + str(reg_train_score))
    print("Testing mean accuracy score: "+ str(reg_test_score))
    print("5-fold cross validation" + str(reg_cv_scores))
    print("Average score: " + str(reg_cv_score) + "\n")
    train_scores['Linear Regression'].append(reg_train_score)
    test_scores['Linear Regression'].append(reg_test_score)
    cv_scores['Linear Regression'].append(reg_cv_score)
    filename_linreg = out_path +'/'+ antibiotic_name + ".linearRegression.pkl"
    with open(filename_linreg, 'wb') as file:
        pickle.dump(reg, file)

    print("2. Ridge Linear regression model\n")
    ridge = Ridge(alpha = 0.1)
    ridge.fit(X_train, y_train)
    ridge_pred = ridge.predict(X_test)
    ridge_train_score = ridge.score(X_train,y_train)
    ridge_test_score = ridge.score(X_test,y_test)
    ridge_cv_scores = cross_val_score(ridge, X_train,y_train, cv=5)
    ridge_cv_score = sum(ridge_cv_scores)/5
    print("Training mean accuracy score:" + str(ridge_train_score))
    print("Testing mean accuracy score: "+ str(ridge_test_score))
    print("5-fold cross validation" + str(ridge_cv_scores))
    print("Average score: " + str(ridge_cv_score) + "\n")
    train_scores['Ridge Regression'].append(ridge_train_score)
    test_scores['Ridge Regression'].append(ridge_test_score)
    cv_scores['Ridge Regression'].append(ridge_cv_score)
    filename_ridge = out_path +'/'+ antibiotic_name + ".ridgeRegression.pkl"
    with open(filename_ridge, 'wb') as file:
        pickle.dump(ridge, file)
    
    print("3.Lasso Linear regression model\n")
    lasso = Lasso(alpha=0.1,tol=1e-5)
    lasso.fit(X_train, y_train)
    lasso_pred = lasso.predict(X_test)
    lasso_train_score = lasso.score(X_train,y_train)
    lasso_test_score = lasso.score(X_test,y_test)
    lasso_cv_scores = cross_val_score(lasso, X_train,y_train, cv=5)
    lasso_cv_score = sum(lasso_cv_scores)/5
    print("Training mean accuracy score:" + str(lasso_train_score))
    print("Testing mean accuracy score: "+ str(lasso_test_score))
    print("5-fold cross validation" + str(lasso_cv_scores))
    print("Average score: " + str(lasso_cv_score) + "\n")
    train_scores['Lasso Regression'].append(lasso_train_score)
    test_scores['Lasso Regression'].append(lasso_test_score)
    cv_scores['Lasso Regression'].append(lasso_cv_score)    
    filename_lasso = out_path +'/'+ antibiotic_name + ".lassoRegression.pkl"
    with open(filename_lasso, 'wb') as file:
        pickle.dump(lasso, file)



    outfile0=open(out_path +'/'+ antibiotic_name + "_reg_coef.txt",'w')
    reg_coef = reg.fit(X, y).coef_
    select_names = []
    select_coef = []

    for i in range(len(reg_coef)):
        if abs(reg_coef[i]) > float(threshold):
            select_names.append(names[i])
            select_coef.append(reg_coef[i])
            outfile0.write(str(names[i])+"\t"+str(reg_coef[i])+"\n")
    outfile0.close()
    os.system("coef2gene.py " + gene_location + " " + out_path +'/'+ antibiotic_name + "_reg_coef.txt " + out_path +'/'+ antibiotic_name + "_reg_coef_out.txt")
    os.system("rm " + out_path +'/'+ antibiotic_name + "_reg_coef.txt ")

    current_plot = plt.plot(range(len(select_names)), select_coef)
    current_plot = plt.xticks(range(len(select_names)), select_names, rotation=90, fontsize=5)
    current_plot = plt.ylabel('Linear regression coefficients')
    current_plot = plt.title("Coefficients above threshold " + str(threshold) + " Vs. Genes, for model of " + antibiotic_name + " resistance")
    #plt.show()
    plt.savefig(out_path + '/' + antibiotic_name + "_reg" + ".png")
    plt.clf()



    outfile1=open(out_path +'/'+ antibiotic_name + "_ridge_coef.txt",'w')
    ridge_coef = ridge.fit(X, y).coef_
    select_names = []
    select_coef = []
    for i in range(len(ridge_coef)):
    # 0.1 is the threshold to view coefficients, users can relax the thresholds to view more.
        if abs(ridge_coef[i]) > float(threshold): 
            select_names.append(names[i])
            select_coef.append(ridge_coef[i])
            outfile1.write(str(names[i])+"\t"+str(ridge_coef[i])+"\n")
    outfile1.close() 
    os.system("coef2gene.py " + gene_location + " " + out_path +'/'+ antibiotic_name + "_ridge_coef.txt " + out_path +"/"+ antibiotic_name + "_ridge_coef_out.txt")
    os.system("rm " + antibiotic_name + "_ridge_coef.txt")

    current_plot = plt.plot(range(len(select_names)), select_coef)
    current_plot = plt.xticks(range(len(select_names)), select_names, rotation=90, fontsize=5)
    current_plot = plt.ylabel('Ridge regression coefficients')
    current_plot = plt.title("Coefficients above threshold " + str(threshold) + " Vs. Genes, for model of " + antibiotic_name + " resistance")
    #plt.show()
    plt.savefig(out_path +'/' + antibiotic_name + "_ridge" + ".png")
    plt.clf() 



    outfile2=open(out_path +'/'+ antibiotic_name + "_lasso_coef.txt",'w')
    lasso_coef = lasso.fit(X, y).coef_
    select_names = []
    select_coef = []

    for i in range(len(lasso_coef)):
        if abs(lasso_coef[i]) > float(threshold):
            select_names.append(names[i])
            select_coef.append(lasso_coef[i])
            outfile2.write(str(names[i])+"\t"+str(lasso_coef[i])+"\n")
    outfile2.close()
    os.system("coef2gene.py " + gene_location + " " + out_path +'/'+ antibiotic_name + "_lasso_coef.txt " + out_path +'/'+ antibiotic_name + "_lasso_coef_out.txt")
    os.system("rm " + out_path +'/'+ antibiotic_name + "_lasso_coef.txt ")

    current_plot = plt.plot(range(len(select_names)), select_coef)
    current_plot = plt.xticks(range(len(select_names)), select_names, rotation=90, fontsize=5)
    current_plot = plt.ylabel('Lasso regression coefficients')
    current_plot = plt.title("Coefficients above threshold " + str(threshold) + " Vs. Genes, for model of " + antibiotic_name + " resistance")
    #plt.show()
    plt.savefig(out_path + '/' + antibiotic_name + "_lasso" + ".png")
    plt.clf()


# Plot all test scores
plt.figure(figsize=(12, 6))

models = list(test_scores.keys())
colors = plt.cm.Set3(np.linspace(0, 1, len(antibiotics)))

for i, antibiotic in enumerate(antibiotics):
    x_positions = np.arange(len(models))
    y_values = [test_scores[model][i] for model in models]
    plt.scatter(x_positions, y_values, label=antibiotic, s=100, alpha=0.7, color=colors[i])

plt.xlabel('Models', fontsize=12)
plt.ylabel('Test Accuracy Score', fontsize=12)
plt.title('Model Performance Comparison Across Antibiotics', fontsize=14)
plt.xticks(np.arange(len(models)), models, rotation=45, ha='right')
plt.legend(title='Antibiotics', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(out_path + '/all_models_comparison.png', dpi=300)
plt.close()


if keep_log:
    f_log.close()
    f_err.close()