#!/usr/bin/env python3

##*************************************************************************##
##          Step4. linear regression models of the AMR data analysis     
##          Edited by Giuli Sucar in Nov 26, 2025          
##*************************************************************************##

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
from sklearn.linear_model import LinearRegression, LogisticRegression, Ridge, Lasso
from sklearn.model_selection import train_test_split, cross_val_score, GridSearchCV
from sklearn.preprocessing import scale
from sklearn.metrics import mean_squared_error
from sklearn import metrics
from sklearn.tree import DecisionTreeRegressor
from sklearn.svm import SVR
from sklearn.preprocessing import StandardScaler
import sys


#if len(sys.argv)!=6: 
#    print("Usage:python3 4.AMR_Learn_linear.py <feature2target.txt> <gene_location_info.txt> <antibiotics.txt> <out_dir_path> \n\n e.g.,python3 AMR_Learn_linear.py feature2target.txt Spectinomycin 0.1")
#    sys.exit()

#  4.AMR_Learn_linear.py PRJNA776899.feature2target.txt GCF_000020105.1.tab PRJNA776899.antibiotics.small.txt learn 0.1

keep_log = True

feature2target = sys.argv[1]
gene_location = sys.argv[2] 
antibiotics_path = sys.argv[3]
out_path = sys.argv[4]

os.makedirs(out_path, exist_ok=True) 

# writing out the log file
if keep_log:
    f_log = open(out_path +'/run.log', 'w')
    sys.stdout = f_log
    f_err = open(out_path +'/run.err', 'w')
    sys.stderr = f_err

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
print("X:")
print(X)
#X = scale(X)
#print("X scaled:")
#print(X)
names = hsd_data.drop(columns=[*antibiotics, 'locus_tag'], axis=1).columns

# Define models in a dictionary for easy addition
models_dict = {
    'Ridge Regression'  : Ridge(alpha=0.1),
    'Lasso Regression'  : Lasso(alpha=0.1, tol=1e-5),
    'Decision Tree'     : DecisionTreeRegressor(),
    'Support Vector Machine': SVR(kernel='linear', C=1.0)
}

train_scores = {model_name: [] for model_name in models_dict.keys()}
test_scores = {model_name: [] for model_name in models_dict.keys()}
cv_scores = {model_name: [] for model_name in models_dict.keys()}

for antibiotic_name in antibiotics:
    print("Training model for resistance to " + antibiotic_name)

    y = hsd_data[antibiotic_name].values
    print(y)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42, shuffle=True)
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test  = scaler.transform(X_test)
    print("scaled X train:")
    print(X_train)
    
    # Train all models sequentially
    for model_name, model in models_dict.items():
        print(f"{model_name}\n")
        
        model.fit(X_train, y_train)
        y_pred_train = model.predict(X_train)
        y_pred_test = model.predict(X_test)
        train_score = np.sqrt(mean_squared_error(y_train, y_pred_train))
        test_score = np.sqrt(mean_squared_error(y_test, y_pred_test))
        cv_fold_scores = np.sqrt(np.abs(cross_val_score(model, X_train, y_train, cv=5, scoring='neg_mean_squared_error')))
        cv_score = cv_fold_scores.mean()
        
        print(f"Training mean RMSE score: {train_score}")
        print(f"Testing mean RMSE score: {test_score}")
        print(f"5-fold cross validation RMSE: {cv_fold_scores}")
        print(f"Cross validation averaged RMSE: {cv_score}\n")
        
        train_scores[model_name].append(train_score)
        test_scores[model_name].append(test_score)
        cv_scores[model_name].append(cv_score)
        
        # Save model
        filename = out_path + '/' + antibiotic_name + f".{model_name.replace(' ', '_').lower()}.pkl"
        with open(filename, 'wb') as file:
            pickle.dump(model, file)

        if model_name in ["Decision Tree", "Support Vector Machine"]:
            continue  

        # Write the top 10 coefficients by absolute value
        outfile = open(out_path + '/' + antibiotic_name + f"_{model_name}_coef.txt", 'w')
        coef = model.coef_
        # ensure 1D array for consistency
        if hasattr(coef, 'ndim') and coef.ndim > 1:
            coef = coef.ravel()

        # get indices of features sorted by absolute coefficient (descending)
        abs_idx = np.argsort(np.abs(coef))[::-1]
        top_n = 10
        top_idx = abs_idx[:top_n]

        select_names = [names[i] for i in top_idx]
        select_coef = [coef[i] for i in top_idx]

        # write: feature_name \t coefficient \t abs(coefficient)
        for i in top_idx:
            outfile.write(f"{names[i]}\t{coef[i]}\t{abs(coef[i])}\n")
        outfile.close()
        
        os.system(f"coef2gene.py {gene_location} {out_path}/{antibiotic_name}_{model_name}_coef.txt {out_path}/{antibiotic_name}_{model_name}_coef_out.txt")
        os.system(f"rm {out_path}/{antibiotic_name}_{model_name}_coef.txt")

        
        plt.figure()
        plt.scatter(range(len(select_names)), select_coef)
        plt.xticks(range(len(select_names)), select_names, rotation=45, fontsize=8)
        plt.ylabel(f'{model_name} coefficients')
        plt.title(f"Coefficients Vs. Genes, for model of {antibiotic_name} resistance")
        plt.savefig(out_path + '/' + antibiotic_name + f"_{model_name}.png")
        plt.close()
        

# Plot all test scores
plt.figure()

models = list(test_scores.keys())
colors = plt.cm.Set3(np.linspace(0, 1, len(antibiotics)))
x_positions = np.arange(len(models))

for i, antibiotic in enumerate(antibiotics):
    y_values = [test_scores[model][i] for model in models]
    plt.scatter(x_positions, y_values, label=antibiotic, s=100, alpha=0.7, color=colors[i])

plt.xlabel('Models', fontsize=12)
plt.ylabel('Test RMSE Score', fontsize=12)
plt.title('Model Test Scores for All Antibiotics', fontsize=14)
plt.xticks(np.arange(len(models)), models, rotation=45, ha='right')
plt.legend(title='Antibiotics', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(out_path + '/all_models_test.png', dpi=300)
plt.close()

plt.figure()
for i, antibiotic in enumerate(antibiotics):
    y_values = [train_scores[model][i] for model in models]
    plt.scatter(x_positions, y_values, label=antibiotic, s=100, alpha=0.7, color=colors[i])

plt.xlabel('Models', fontsize=12)
plt.ylabel('Train RMSE Score', fontsize=12)
plt.title('Model Train Scores for All Antibiotics', fontsize=14)
plt.xticks(np.arange(len(models)), models, rotation=45, ha='right')
plt.legend(title='Antibiotics', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(out_path + '/all_models_train.png', dpi=300)
plt.close()

plt.figure()
for i, antibiotic in enumerate(antibiotics):
    y_values = [cv_scores[model][i] for model in models]
    plt.scatter(x_positions, y_values, label=antibiotic, s=100, alpha=0.7, color=colors[i])

plt.xlabel('Models', fontsize=12)
plt.ylabel('Cross Validation RMSE Average', fontsize=12)
plt.title('Model Cross Validation Average RMSE for All Antibiotics', fontsize=14)
plt.xticks(np.arange(len(models)), models, rotation=45, ha='right')
plt.legend(title='Antibiotics', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(out_path + '/all_models_cv.png', dpi=300)
plt.close()



if keep_log:
    f_log.close()
    f_err.close()