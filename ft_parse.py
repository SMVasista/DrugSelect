from __future__ import division
from mysql.connector import connection
import math, statistics, os, re, random
from datetime import date,datetime
import futils
import itertools
import llp_solver
import numpy as np

# make sure mysql is running on port 3306 - default port

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
            #print '=======', config
            DbUtils.cnx = connection.MySQLConnection(user=config["user"], password=config["password"], host=config["host"], database=config["database"])

    @classmethod
    def cursor(cls):
        return DbUtils.cnx.cursor()

    @classmethod
    def commit(cls):
        DbUtils.cnx.commit()

def getKeys(flag, query):
    return True
    else:
        return False

def readUnit(uids, method):
    if method == 'db':
        uData = {}
        DbUtils()
        cursor = DbUtils.cursor()
        for uid in uids:
            query = "SELECT * FROM LOG_INPUT WHERE UNIT_ID = \""+str(uid)+"\""
            rows = DbUtils.query(cursor, query)

            uData[uid] = {'indc': rows[0][2]}
            uData[uid]['mut'] = [(row[3], row[5]) for row in rows if getKeys('mutation', row[4]) and row[6] == 'NULL']
            uData[uid]['cna'] = [(row[3], row[5]) for row in rows if getKeys('CNA', row[4]) and row[6] == 'NULL']
            uData[uid]['pdata'] = [(row[6], row[7]) for row in rows if row[3] == 'NULL']

    elif method == 'file':
        uData = {}
        for uid in uids:
            fileLoc = input("uid data location> ")
            try:
                data = futils.readLinesAndSplit(fileLoc, ',')
                
                uData[uid] = {'indc': rows[0][2]}
                uData[uid]['mut'][row[3]] : row[5] for row in data if getKeys('mutation', row[4]) and row[6] == 'NULL'
                uData[uid]['cna'][row[3]] : row[5] for row in data if getKeys('CNA', row[4]) and row[6] == 'NULL'
                uData[uid]['pdata'][row[6]] : row[7] for row in data if row[3] == 'NULL'
            except:
                raise

    return uData

def seggregate_p_g_components(uData):
    #Loading indication composition
    DbUtils()
    cursor = DbUtils.cursor()
    query = "SELECT * FROM LOG_INPUT WHERE 1"
    indc_rows = DbUtils.query(cursor, query)
        
    #Reading P-G mapping
    query = "SELECT * FROM GP_MAP WHERE 1"
    key_rows = DbUtils.query(cursor, query)

    for uid in uData:
        #Mapping genes to pathways
        uData[uid]['n_pathways'] = uData[uid]['pdata']
        uData[uid]['n_genes'] = {}

        remove = []

        #Updating Mut/Cna information in p_g transition matrix
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
                            uData[uid]['n_pathways'][k[2]] += float(k[3])*G_abr_wt[uData[uid][atype][elem]]*I_wt
                            remove.append(elem)
                        else:
                            uData[uid]['n_pathways'][k[2]] = float(k[3])*G_abr_wt[uData[uid][atype][elem]]*I_wt
                            remove.append(elem)
                for elem in uData[uid][atype].keys():
                    if elem not in remove:
                        uData[uid]['n_genes'][elem] = 0.1*G_abr_wt[uData[uid][atype][elem]]*I_wt

        #Identifying unique elements
        uPaths = list(set([k[2] for k in key_rows]))

        uGenes = []
        m = [uData[uid].keys() for uid in uData.keys()]
        for elem in m:
            for i in range(len(elem)):
                uGenes.append(elem[i])
            
    return uData, uGenes, uPaths

def composeMatrix(C_list, H_list, uData, method):

    ##########################################################################
    #Composition Rules
    #(Modifications not recommended unless stated otherwise)
    ##########################################################################
    #
    ###Gene->Pathway rules
    #Positive values mapped as directly correlating
    #Amplifications mapped to positive
    #Deletions mapped to negative
    #G_weight = (1+g_I_weight)*G_abr_wt
    ##########################################################################
    ###Drug->Pathway rules
    #Positive values mapped for as directly correlating
    #Weightage Scaling factor: 1/(n+0.5) where n = rank
    #Inhibition effect factor: p - wsf*dw
    #Direct: rank 1, Indirect: rank 2, Undefined: rank 5
    #Healing matrix scoring method: v = 0.5 + k/(k+1) -> where k is the number
    #of drugs in combination
    #
    ##########################################################################
    #Parameters (Do not change)
    #Gain-Loss weightage
    atypes = ['cna', 'mut']
    tt_wt = {'direct': 1.0, 'indirect': 0.5}
    tta_wt = {'inhibit': -0.7, 'activate': 0.5}
    #Amp-KD weightage
    ##########################################################################

    ###This parameter can be modified without affecting performace of solver.
    G_abr_wt = {'gain': 1.0, 'loss': -1.0,'amp': 1.0, 'del', -1.0, 'other': 1.0}

    if method == 'db':

        BIN = {}
        #Loading C-list data
        for C_element in C_list:
            query = "SELECT * FROM G_CHE_ASSO WHERE CHEMO_DRG = \""+str(C_element)+"\""
            cg_drg_rows = DbUtils.query(cursor, query)
            query = "SELECT * FROM P_CHE_ASSO WHERE CHEMO_DRG = \""+str(C_element)+"\""
            cp_drg_rows = DbUtils.query(cursor, query)
            BIN[C_element] = cg_drg_rows + cp_drg_rows

        #Loading H-list data
        for H_element in H_list:
            query = "SELECT * FROM G_HCOMP_ASSO WHERE HCOMP = \""+str(H_element)+"\""
            hg_drg_rows = DbUtils.query(cursor, query)
            query = "SELECT * FROM P_HCOMP_ASSO WHERE HCOMP = \""+str(H_element)+"\""
            hp_drg_rows = DbUtils.query(cursor, query)
            BIN[H_element] = hg_drg_rows + hp_drg_rows

    elif method == 'file':
        BIN = []
        G_ASSO = input('Location of drug-on-gene assocation data matrix (all drugs + chemo)> ')
        g_drg_rows = futils.readLinesAndSplit(G_ASSO, ',')
        P_ASSO = input('Location of drug-on-pathway assocation data matrix (all drugs + chemo)> ')
        p_drg_rows = futils.readLinesAndSplit(P_ASSO, ',')

        BIN = {}
        #Seggregating
        for line in g_drg_rows+p_drg_rows:
            if line[1] in BIN:
                BIN[line[1]].append(line)
            else:
                BIN[line[1]] = [line]

    #Creating G-P transition matrix for each uid
    uData_, uGene_, uPaths_ = seggregate_p_g_components(uData)

    #Creating drug combinations
    Combinations = {}
    Combinations[1] = [r for r in H_list]
    #Cn-Hn combinations = C(n) * Hn
    for i in range(2, 5):
        Combinations[i] = [elem for elem in itertools.combinations(H_list, i)]

    #Creating DXP matrix
    xrg = len(C_list)*len(Combinations.values()) + 1
    yrg = len(uGene_) + len(uPaths_)

    uElem = uGene_ + uPaths_

    IDX = {}
    for uid in uids:
        IDX[uid] = {}
        
        #Element-wise solving for each pathway to combination scenrios
        for i in Combinations.keys():
            IDX[uid][i] = ([], [], [])
            for ther in Combinations[i]:
                for chemo in C_list: #x
                        therapy_name = str(chemo)
                        therapy_name = therapy_name + str(f)+'_' for f in ther
                    for element in sorted(uElem): #y
                        IDX[uid][i][0].append(therapy_name)
                        IDX[uid][i][1].append(element)
                        computed_score = 0
                        for d_elem in ther:
                            #Calculating effect of drug on gene/pathway
                            g_score = 0
                            g_score += tt_wt(k[3])*tta_wt(k[4])*(1/(k[5]+1)) for k in BIN[d_elem] if k[1] == element
                            #Calculating residual effect due to abberation
                            impact = [k[5] for k in BIN[chemo] if element == k[2] else 0][0]
                            computed_score = (uData[uid][element] - impact * g_score)
                        IDX[uid][i][2].append(computed_score)
    return IDX
                
if __name__=="__main__":
    #Process-Flow: 
    uids, c_list, h_list = parse_input()
    uD = readUnit(uids, 'db')
    X = composeMatrix(c_list, h_list, uD,'db')
    sX = llp_solver.psvd(X, end_points='auto')
    #plotter.report(sX)
    
        
        
