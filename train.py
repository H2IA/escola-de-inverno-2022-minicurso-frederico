from features import smiles_to_morganfingerprints
from sklearn.model_selection import train_test_split
from imblearn.under_sampling import RandomUnderSampler
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.base import BaseEstimator
from sklearn.metrics import classification_report
from typing import Tuple
from argparse import ArgumentParser
import pandas as pd
import bentoml
import pickle

def main():
    argument_parser = ArgumentParser()
    argument_parser.add_argument('--input', default='data/features/b3db.csv')
    argument_parser.add_argument('--output', default='data/models/model.pickle')
    arguments = argument_parser.parse_args()

    print('Training model ... ', end='')

    df_features = pd.read_csv(arguments.input)
    model, report = train(df_features)

    print('done.')

    with open(arguments.output, 'wb') as writer:
        writer.write(pickle.dumps(model))
    print(f"Model file: {arguments.output}")

    saved_model = bentoml.sklearn.save_model('bbb-model', model)
    print(f"Model tag: {saved_model}")

    print('Classification report with test data:')
    print(report)

def train(df_features:pd.DataFrame) -> Tuple[BaseEstimator, str]:
    X = df_features.drop(['label'], axis=1)
    y = df_features['label']

    X_rus, y_rus = RandomUnderSampler().fit_resample(X, y)
    X_train, X_test, y_train, y_test = train_test_split(X_rus, y_rus)

    model = ExtraTreesClassifier()
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    report = classification_report(y_test, y_pred)

    return model, report

if __name__ == '__main__':
    main()