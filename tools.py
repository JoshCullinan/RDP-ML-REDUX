# RDP Tensorflow Impl.
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix, ConfusionMatrixDisplay


def class_report(y_true, y_preds):
    print(
        classification_report(y_true,
                              y_preds,
                              digits= 4,
                              target_names=["Parent", "Recombinant"]))
    
    # disp1 = ConfusionMatrixDisplay(
    #     confusion_matrix(y_true=y_true, y_pred=y_preds, normalize='true'),
    #     display_labels=["Parent", "Recombinant"],
    # )

    # disp1 = ConfusionMatrixDisplay.from_predictions(
    # y_true=y_true, y_pred=y_preds, 
    # normalize='true', 
    # cmap='PuBu', 
    # values_format= '.4g',
    # display_labels=["Parent", "Recombinant"],
    # text_kw = {'size': 11},
    # )
    # disp1.plot()


def flip(series):
    split = series.str.split(",", expand=True).values  # .astype(float)

    return split.reshape(len(series) * 3)

def ingestor(recom_path):
    
    recom = pd.read_csv(recom_path, sep="\t")

    # def f(x): return x.strip("()")
    # for (colname, _) in recom.iteritems():
    #     recom.loc[:, colname] = recom.loc[:, colname].apply(f)
    
    allData = recom.apply(flip, axis=0) 

    # Train, Test = train_test_split(recom, test_size=0.2)

    # Convert all entries from str to list
    # Train = Train.apply(flip, axis=0)
    # Test = Test.apply(flip, axis=0)

    return allData
