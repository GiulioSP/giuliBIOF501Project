##*************************************************************************##
##          Step4. linear regression models of the AMR data analysis       ##          
##*************************************************************************##


# pip3 install -U scikit-learn scipy matplotlib
#while read line; do python3 4.AMR_Learn.py feature2target.txt $line;done <list.txt

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
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
    print("Usage:python3 AMR_Learn_linear.py <feature2target.txt> <gene_location_info.txt> <antibiotics.txt> <out_dir_path> <threshold for filtering absolute coefficient>\n\n e.g.,python3 AMR_Learn_linear.py feature2target.txt Spectinomycin 0.1")
    sys.exit()

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
f_log = open(out_path +'/run.log', 'w')
sys.stdout = f_log
f_err = open(out_path +'/run.err', 'w')
sys.stderr = f_err

for antibiotic_name in antibiotics:
    print("Training model for resistance to " + antibiotic_name)

    y = hsd_data[antibiotic_name].values

    #fitting a regression model
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.3, random_state=42)

    print("Supervised regression models\n")
    print("1.Linear regression model\n")
    reg = LinearRegression()
    reg.fit(X_train,y_train)
    y_pred = reg.predict(X_test)
    print("Training data:" + str(reg.score(X_train,y_train)))
    print("Testing data: "+ str(reg.score(X_test,y_test)))
    cv_results = cross_val_score(reg, X, y, cv=5)
    print("5-fold cross validation" + str(cv_results))
    avg_score = sum(cv_results)/5
    print("Average score: " + str(avg_score)+"\n")

    #print(confusion_matrix(y_test, y_pred))
    # print(classification_report(y_test, y_pred))

    print("2. Ridge Linear regression model\n")
    ridge = Ridge(alpha = 0.1)
    ridge.fit(X_train, y_train)
    ridge_pred = ridge.predict(X_test)
    ridge.score(X_test, y_test)
    print("Training data: " + str(ridge.score(X_train,y_train)))
    print("Testing data: " + str(ridge.score(X_test,y_test)))
    cv_results = cross_val_score(ridge, X, y, cv=5)
    print("5-fold cross validation" + str(cv_results))
    avg_score = sum(cv_results)/5
    print("Average score: " + str(avg_score)+"\n")

    #print(confusion_matrix(y_test, ridge_pred))
    # print(classification_report(y_test, ridge_pred))
    
    print("3.Lasso Linear regression model\n")
    lasso = Lasso(alpha=0.1,tol=1e-5)
    lasso.fit(X_train, y_train)
    lasso_pred = lasso.predict(X_test)
    lasso.score(X_test, y_test)
    print("Training data: " + str(lasso.score(X_train,y_train)))
    print("Testing data: "+ str(lasso.score(X_test,y_test)))
    cv_results = cross_val_score(lasso, X, y, cv=5)
    print("5-fold cross validation" + str(cv_results))
    avg_score = sum(cv_results)/5
    print("Average score: " + str(avg_score)+"\n")

    outfile=open(out_path +'/'+ antibiotic_name + "_ridge_coef.txt",'w')
    ridge_coef = ridge.fit(X, y).coef_
    select_names = []
    select_coef = []
    for i in range(len(ridge_coef)):
    # 0.1 is the threshold to view coefficients, users can relax the thresholds to view more.
        if abs(ridge_coef[i]) > float(threshold): 
            select_names.append(names[i])
            select_coef.append(ridge_coef[i])
            outfile.write(str(names[i])+"\t"+str(ridge_coef[i])+"\n")
    outfile.close() 
    os.system("python3 coef2gene.py " + gene_location + " " + out_path +'/'+ antibiotic_name + "_ridge_coef.txt " + out_path +"/"+ antibiotic_name + "_ridge_coef_out.txt")
    #os.system("rm " + antibiotic_name + "_ridge_coef.txt")

    current_plot = plt.plot(range(len(select_names)), select_coef)
    current_plot = plt.xticks(range(len(select_names)), select_names, rotation=90, fontsize=5)
    current_plot = plt.ylabel('Ridge regression coefficients')
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
    os.system("python3 coef2gene.py " + gene_location + " " + out_path +'/'+ antibiotic_name + "_lasso_coef.txt " + out_path +'/'+ antibiotic_name + "_lasso_coef_out.txt")
    #os.system("rm " + out_path +'/'+ antibiotic_name + "_lasso_coef.txt ")

    current_plot = plt.plot(range(len(select_names)), select_coef)
    current_plot = plt.xticks(range(len(select_names)), select_names, rotation=90, fontsize=5)
    current_plot = plt.ylabel('Lasso regression coefficients')
    #plt.show()
    plt.savefig(out_path + '/' + antibiotic_name + "_lasso" + ".png")
    plt.clf()

f_log.close()
f_err.close()