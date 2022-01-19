#! /bin/python3

import requests 
import json
import csv



def read_csv(data_file):
    datafiles = []
    with open(data_file, newline='') as csv_names:
        for line in csv.reader(csv_names, delimiter=',', quotechar='|'):
            datafiles.append(line[0])
    return datafiles   


datafiles = read_csv('datafiles.csv')


for name in datafiles[1:]:
    print(name)
    r = requests.get('http://mskilab.com/gGraph/json/{}'.format(name), verify = False)
    data = r.json()
    print(name)
    with open('graphs/'+name, 'w') as outfile:
        json.dump(data, outfile, indent=4)
