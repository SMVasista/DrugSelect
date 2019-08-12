import re
import json
import pickle
import os
import csv

def readLines(fileName):
    with open(fileName, 'r') as f:
        return  [line.strip() for line in f.readlines() if len(line.strip()) > 0]

def readLinesAndSplit(fileName, sep):
    return [map(lambda token: token.strip(), line.split(sep)) for line in readLines(fileName)]

def filter(collection, val):
    if isinstance(collection, dict):
        if type(val).__name__ == 'function':
            return dict((k, v) for k,v in collection.iteritems() if val(v))
        else:
            return dict((k, v) for k,v in collection.iteritems() if v == val)
    elif isinstance(collection, list):
        if type(val).__name__ == 'function':
            return  [elem for elem in collection if val(elem)]
        else:
            return  [elem for elem in collection if val == elem]

def compact(collection):
    if isinstance(collection, list):
        return [ elem for elem in collection if elem != None ]
    elif isinstance(collection, dict):
        return dict((k, v) for k,v in collection.iteritems() if v != None)

# take only the key, value pairs from given dict where key is in the given list

def pick(obj, keylist):
    return dict((k, v) for k,v in obj.iteritems() if k in keylist)


# check if dir exists, if not create it - functionality akin to mkdir -p (available in python 3 but not in 2)

def jsonLoad(fileName):
    with open(fileName, 'r') as fp:
        return json.loads(fp.read())

def jsonDump(fileName, obj):
    with open(fileName, 'w') as fp:
        fp.write(json.dumps(obj, indent=2))
