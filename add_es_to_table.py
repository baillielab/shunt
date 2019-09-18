#!/usr/bin/python
# -*- coding: UTF-8 -*-

'''
add a new column labelled 'es' to a data table in csv or tsv format
'''

import pandas as pd
from shunt import es
pd.options.mode.chained_assignment = None

def add_es_to_df(thisdf, g='kPa', a='pH'):
    thisdf.columns = [x.lower() for x in thisdf.columns]
    if g.lower()=='mmhg':
        thisdf['pao2'] = thisdf['pao2'].multiply(101.325/760)
        thisdf['paco2'] = thisdf['paco2'].multiply(101.325/760)
    if a.lower()=='h':
        thisdf['ph'] = thisdf['ph'] ## *** need to add conversion to shunt.py too
    thisdf['es'] =  thisdf.apply(lambda row: es(fio2=row['fio2'], pao2=row['pao2'], paco2=row['paco2'], pH=row['ph'] ), axis=1 )
    return (thisdf)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--filename', help='filename')
    parser.add_argument('-gasunit', default='kPa', type=str, help='kPa or mmHg')
    parser.add_argument('-acidunit', default='pH', help='pH or H')
    args = parser.parse_args()

    df = pd.read_csv(args.filename)
    df = add_es_to_df(df, g=args.gasunit, a=args.acidunit)
    df.to_csv(args.filename+'.es.txt')