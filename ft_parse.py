from __future__ import division
from mysql.connector import connection
import math, sys, os, re, random
from datetime import date,datetime
import futils
import itertools
import llp_solver
import numpy as np

# make sure mysql is running on port 3306 - default port
##########################################################################
#Matrix Composition Rules                                                #
#(Modifications not recommended unless stated otherwise)                 #
##########################################################################
#                                                                        #
###Gene->Pathway rules                                                   #
#Positive values mapped as directly correlating                          #
#Amplifications mapped to positive                                       #
#Deletions mapped to negative                                            #
#G_weight = (1+g_I_weight)*G_abr_wt                                      #
##########################################################################
###Drug->Pathway rules                                                   #
#Positive values mapped for as directly correlating                      #
#Weightage Scaling factor: 1/(n+0.5) where n = rank                      #
#Inhibition effect factor: p - wsf*dw                                    #
#Direct: rank 1, Indirect: rank 2, Undefined: rank 5                     #
#Healing matrix scoring method: v = 0.5 + k/(k+1) -> where k is the      #
#number of drugs in combination                                          #
#                                                                        #
##########################################################################
    
#Parameters (Do not change)###############################################
tt_wt = {'direct': 1.0, 'indirect': 0.5}                                 #
tta_wt = {'inhibit': -0.7, 'activate': 0.5}                              #
imp_wt = {'accelerate': 1.0, 'escape': -1.0}                             #
##########################################################################
###This parameter can be modified without affecting performace of solver.
G_abr_wt = {'gain': 1.0, 'loss': -1.0,'amp': 1.0, 'del': -1.0, 'other': 1.0}

class DbUtils:
    cnx = None

    @classmethod
    def readConfig(cls, cfgFile='./config.ext'):
        config = {}
        lines = futils.readLines(cfgFile)
        for line in lines:
            if '=' in line and line[0] != '#':
                key,value = line.split('=')
                config[key.strip()] = value.strip()
        return config

    def __init__(self, cfgFile='./config.ext'):
        if DbUtils.cnx == None:
            config = DbUtils.readConfig(cfgFile)
            DbUtils.cnx = connection.MySQLConnection(user=config["user"], password=config["password"], host=config["host"], database=config["database"])

    @classmethod
    def cursor(cls):
        return DbUtils.cnx.cursor()

    @classmethod
    def commit(cls):
        DbUtils.cnx.commit()

def parseInput(fileLoc):
    lines = futils.readLinesAndSplit(fileLoc, '=')
    uids = []
    c_list = []
    h_list = []
    method = 'db'
    for line in lines:
        if line[0] == 'UID' or line[0] == 'uid' or line[0] == 'uids' or line[0] == 'UIDS':
            uids = line[1].split(',')
        if 'chemo' in line[0] and 'non' not in line[0]:
            c_list = line[1].split(',')
        if 'non' in line[0]:
            h_list = line[1].split(',')
        if line[0] == 'METHOD' or line[0] == 'method' or line[0] == 'Method':
            method = line[1].split(',')[0]
    return uids, c_list, h_list, method

def getKeys(flag, query):
    F = {'mutation': ['mut', 'MUT', 'SNP', 'SNV', 'mutated', 'abberated'], 'CNA': ['cnv', 'CNV', 'cna', 'CNA', 'amp', 'del', 'amplified', 'deleted', 'deletion', 'deep']}
    if query in F[flag]:
        return True
    else:
        return False

def reDesignTherapyName(chemo, therapy_set):
    t_name = str(chemo)
    for t in therapy_set:
        t_name = t_name+'_'+str(t)
    return t_name

def numRet(value):
    if type(value) != 'str':
        return float(value)
    else:
        return 0.0

def readUnitData(uids, method):
    if method == 'db':
        uData = {}
        DbUtils()
        cursor = DbUtils.cursor()
        for uid in uids:
            query = "SELECT * FROM LOG_INPUT WHERE UNIT_ID = \""+str(uid)+"\""
            rows = DbUtils.query(cursor, query)
            uData[uid] = {'indc': rows[0][2], 'mut': {}, 'cna': {}, 'pdata': {}}
            for row in rows:
                if getKeys('mutation', row[4]) and row[6] == 'NULL':
                    uData[uid]['mut'][row[3]] = row[5]
                if getKeys('mutation', row[4]) and row[6] == 'NULL':
                    uData[uid]['cna'][row[3]] = row[5]
                if row[3] == 'NULL' and row[6] != 'NULL':
                    uData[uid]['pdata'][row[6]] = row[7]

    elif method == 'file':
        uData = {}
        for uid in uids:
            fileLoc = input("uid data location> ")
            try:
                rows = futils.readLinesAndSplit(fileLoc, ',')
                uData[uid] = {'indc': rows[0][2], 'mut': {}, 'cna': {}, 'pdata': {}}
                for row in rows:
                    data = futils.readLinesAndSplit(fileLoc, ',')
                    if getKeys('mutation', row[4]) and row[6] == 'NULL':
                        uData[uid]['mut'][row[3]] = row[5]
                    if getKeys('mutation', row[4]) and row[6] == 'NULL':
                        uData[uid]['cna'][row[3]] = row[5]
                    if row[3] == 'NULL' and row[6] != 'NULL':
                        uData[uid]['pdata'][row[6]] = row[7]
            except:
                raise

    return uData

def seggregate_p_g_components(uData):
    #Loading indication composition
    try:
        DbUtils()
        cursor = DbUtils.cursor()
        query = "SELECT * FROM INDC_COMPOSITION WHERE 1"
        indc_rows = DbUtils.query(cursor, query)
    except:
        fileLoc = input('Indication composition file to proceed> ')
        indc_rows = futils.readLinesAndSplit(fileLoc, ',')
        
    #Reading P-G mapping
    try:
        query = "SELECT * FROM GP_MAP WHERE 1"
        key_rows = DbUtils.query(cursor, query)
    except:
        fileLoc = input('G-P file to proceed> ')
        key_rows = futils.readLinesAndSplit(fileLoc, ',')        

    for uid in uData:
        #Mapping genes to pathways
        uData[uid]['n_pathways'] = uData[uid]['pdata']
        uData[uid]['n_genes'] = {}

        remove = []

        #Updating Mut/Cna information in p_g transition matrix
        atypes = ['cna', 'mut']
        for atype in atypes:
            for elem in uData[uid][atype].keys():
                try:
                    i = 3 if atype == 'mut' else 4
                    I_wt = float([r[i] for r in indc_rows if r[1] == uData[uid]['indc'] and r[2] == elem][0])
                except:
                    I_wt = 0.25
                for k in key_rows:
                    if elem == k[1]:
                        if k[2] in uData[uid]['n_pathways']:
                            uData[uid]['n_pathways'][k[2]] = float(uData[uid]['n_pathways'][k[2]]) + float(k[3])*G_abr_wt[uData[uid][atype][elem]]*(float(1.0) + I_wt)
                            remove.append(elem)
                        else:
                            uData[uid]['n_pathways'][k[2]] = float(k[3])*G_abr_wt[uData[uid][atype][elem]]*(float(1.0) + I_wt)
                            remove.append(elem)
                for elem in uData[uid][atype].keys():
                    if elem not in remove:
                        uData[uid]['n_genes'][elem] = 0.1*G_abr_wt[uData[uid][atype][elem]]*I_wt

        #Creating new abstraction = genes + pathways
        uData[uid]['n_abstraction'] = {}
        uData[uid]['n_abstraction'].update(uData[uid]['n_genes'])
        uData[uid]['n_abstraction'].update(uData[uid]['n_pathways'])

    uElem = []
    m = [uData[uid]['n_abstraction'].keys() for uid in uData.keys()]
    for elem in m:
        for i in range(len(elem)):
            uElem.append(elem[i])
            
    return uData, list(set(uElem))

def composeMatrix(C_list, H_list, uData, method):
    if method == 'db':

        C_BIN = {}
        #Loading C-list data
        for C_element in C_list:
            query = "SELECT * FROM G_CHE_ASSO WHERE CHEMO_DRG = \""+str(C_element)+"\""
            cg_drg_rows = DbUtils.query(cursor, query)
            query = "SELECT * FROM P_CHE_ASSO WHERE CHEMO_DRG = \""+str(C_element)+"\""
            cp_drg_rows = DbUtils.query(cursor, query)
            C_BIN[C_element] = cg_drg_rows + cp_drg_rows

        #Loading H-list data
        H_BIN = {}
        for H_element in H_list:
            query = "SELECT * FROM G_HCOMP_ASSO WHERE HCOMP = \""+str(H_element)+"\""
            hg_drg_rows = DbUtils.query(cursor, query)
            query = "SELECT * FROM P_HCOMP_ASSO WHERE HCOMP = \""+str(H_element)+"\""
            hp_drg_rows = DbUtils.query(cursor, query)
            H_BIN[H_element] = hg_drg_rows + hp_drg_rows

    elif method == 'file':
        BIN = []
        G_ASSO = input('Location of drug-on-gene/pathway assocation data matrix (all drugs + chemo)> ')
        cg_drg_rows = futils.readLinesAndSplit(G_ASSO, ',')
        P_ASSO = input('Location of drug-target assocation data matrix (all drugs + chemo)> ')
        hp_drg_rows = futils.readLinesAndSplit(P_ASSO, ',')

        C_BIN = {}
        #Seggregating
        for line in cg_drg_rows:
            if line[1] in C_BIN:
                C_BIN[line[1]].append(line)
            else:
                C_BIN[line[1]] = [line]
        H_BIN = {}
        #Seggregating
        for line in hp_drg_rows:
            if line[1] in H_BIN:
                H_BIN[line[1]].append(line)
            else:
                H_BIN[line[1]] = [line]
                
    #Creating G-P transition matrix for each uid
    uData_, uElem_ = seggregate_p_g_components(uData)

    #Creating drug combinations
    Combinations = {}
    Combinations[1] = [[r] for r in H_list]
    #Cn-Hn combinations = C(n) * Hn
    cmin = 2 if len(H_list) > 2 else len(H_list)
    cmax = len(H_list) if len(H_list) < 5 else 5
    for i in range(cmin, cmax):
        Combinations[i] = [elem for elem in itertools.combinations(H_list, i)]

    #Creating DXP matrix
    IDX = {}
    for uid in uids:
        IDX[uid] = {}
        
        #Element-wise solving for each pathway to combination scenrios
        for i in Combinations.keys():
            IDX[uid][i] = ([], [], [])
            for ther in Combinations[i]:
                for chemo in C_list:
                    therapy_name = reDesignTherapyName(chemo, ther)
                    for element in sorted(uElem_):
                        computed_score = 0
                        for d_elem in ther:
                            #Calculating effect of drug on gene/pathway
                            g_score = 0
                            for k in H_BIN[d_elem]:
                                if k[2] == element:
                                    g_score += tt_wt[k[3]]*tta_wt[k[4]]*(1/(float(k[5])+1))
                            #Calculating residual effect due to abberation
                            impact = 0
                            for k in C_BIN[chemo]:
                                if k[2] == element:
                                    impact += float(k[5])
                            if element in uData[uid]['n_abstraction']:
                                computed_score = (float(uData[uid]['n_abstraction'][element])) - impact*g_score
                        IDX[uid][i][2].append(computed_score)
                    IDX[uid][i][1].append(therapy_name)
            for elem in sorted(uElem_):
                IDX[uid][i][0].append(elem)
    return IDX
                
if __name__=="__main__":
    script, inputfile = sys.argv
    uids, c_list, h_list, method = parseInput(inputfile)
    uD = readUnitData(uids, method='file')
    X = composeMatrix(c_list, h_list, uD, method='file')
    sX = llp_solver.psvd(X, end_points='auto')
    for elem in sX:
        for comp in sX[elem]:
            for q in sX[elem][comp]:
                print q, sX[elem][comp][q]
    #plotter.report(sX)
