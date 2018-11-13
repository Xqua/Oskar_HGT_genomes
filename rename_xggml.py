#!/usr/bin/env python
import sys
import xmltodict
import urllib2
import os
import collections
import pandas as pd
from Bio import SeqIO
uniprot_dl = []


def decorate(path):
    """Decorate the network with metadata."""
    lines = open(path).readlines()
    out = open(path.split('.')[0] + '.decorated.xgmml', 'w')
    for line in lines:
        Pass = False
        if 'Description' in line:
            pID = line.split('value="')[1].split('/')[0]
            try:
                x = uniprot_meta(pID)
            except:
                Pass = True
            if not Pass:
                taxa = x['uniprot']['entry']['organism']['lineage']['taxon']
                k = taxa[0][0]
                if '#text' not in x['uniprot']['entry']['organism']['name']:
                    #             if len(x['uniprot']['entry']['organism']['name']) != 1:
                    try:
                        specie = x['uniprot']['entry']['organism']['name'][0]['#text']
                    except:
                        print x['uniprot']['entry']['organism']['name']
                else:
                    specie = x['uniprot']['entry']['organism']['name']['#text']
                if 'recommendedName' in x['uniprot']['entry']['protein']:
                    #                 if 'fullName' in x['uniprot']['entry']['protein']['recommendedName']:
                    if 'fullName' in x['uniprot']['entry']['protein']['recommendedName']:
                        if '#text' in x['uniprot']['entry']['protein']['recommendedName']['fullName']:
                            fullname = x['uniprot']['entry']['protein']['recommendedName']['fullName']['#text']
                        else:
                            fullname = x['uniprot']['entry']['protein']['recommendedName']['fullName']
                    else:
                        fullname = None
                    if 'shortName' in x['uniprot']['entry']['protein']['recommendedName']:
                        shortname = x['uniprot']['entry']['protein']['recommendedName']['shortName']['#text']
                    else:
                        shortname = None
                else:
                    if len(x['uniprot']['entry']['protein']['submittedName']) != 1 and 'ecNumber' not in x['uniprot']['entry']['protein']['submittedName']:
                        best = 0
                        best_i = 0
                        for i in range(len(x['uniprot']['entry']['protein']['submittedName'])):
                            #                         try:
                            e = int(x['uniprot']['entry']['protein']['submittedName'][i]['fullName']['@evidence'])
        #                         except:
        #                             err = x['uniprot']['entry']['protein']['submittedName']
        #                             print len(err), err.keys()
        #                             break
                            if e > best:
                                best = e
                                best_i = i
                        fullname = x['uniprot']['entry']['protein']['submittedName'][best_i]['fullName']['#text']
                        shortname = None
                    else:
                        fullname = x['uniprot']['entry']['protein']['submittedName']['fullName']['#text']
                        shortname = None
        #             print leaf.fullname
                if k != 'B':
                    if k != 'E':
                        k = 'U'
                if 'Arthropoda' in taxa:
                    k = 'Ar'
                if k == 'Ar':
                    features = ''
                    if 'feature' in x['uniprot']['entry']:
                        feature = x['uniprot']['entry']['feature']
                        if type(feature) != collections.OrderedDict:
                            for i in feature:
                                if i['@type'] == 'domain':
                                    features += i['@description'].lower()
                        else:
                            if feature['@type'] == 'domain':
                                features = feature['@description'].lower()
                    if 'osk' in fullname.lower():
                        k = 'O'
                    elif 'osk' in features:
                        k = 'O'
                    elif Test_Oskar(pID):
                        k = 'O'
                if shortname:
                    n = shortname
                else:
                    n = fullname

                out.write('<att name="Prot_ID" type="string" value="%s" />\n' % pID)
                out.write('<att name="Specie" type="string" value="%s" />\n' % specie)
                out.write('<att name="Kingdom" type="string" value="%s" />\n' % taxa[0])
                out.write('<att name="Prot_FullName" type="string" value="%s" />\n' % fullname)
                out.write('<att name="Prot_ShortName" type="string" value="%s" />\n' % shortname)
                out.write('<att name="Name" type="string" value="%s" />\n' % n)
                out.write('<att name="k" type="string" value="%s" />\n' % k)
            else:
                out.write('<att name="Prot_ID" type="string" value="%s" />\n' % None)
                out.write('<att name="Specie" type="string" value="%s" />\n' % None)
                out.write('<att name="Kingdom" type="string" value="%s" />\n' % None)
                out.write('<att name="Prot_FullName" type="string" value="%s" />\n' % None)
                out.write('<att name="Prot_ShortName" type="string" value="%s" />\n' % None)
                out.write('<att name="Name" type="string" value="%s" />\n' % None)
                out.write('<att name="k" type="string" value="%s" />\n' % None)

        else:
            out.write(line)


def uniprot_meta(accID):
    """Download MEtadata from uniprot."""
    if not os.path.isfile('%s.xml' % accID):
        url = "http://www.uniprot.org/uniprot/%s.xml" % (accID)
        handle = urllib2.urlopen(url)
        xml = handle.read()
        f = open('%s.xml' % accID, 'w')
        f.write(xml)
        f.close()
    else:
        xml = open('%s.xml' % accID).read()
    xml = xmltodict.parse(xml)
    return xml


def Parse_HMMER_output(path):
    """Parse table that hemmer outputs."""
    f = open(path)
    lines = f.readlines()
    f.close()
    collumns = ['target name', 'Prot_ID', 'Specie_ID', 'accession', 'query name', 'accession', 'Pre E-value', 'Pre score', 'Pre bias', 'E-value', 'score', 'bias', 'exp', 'reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc', 'description of target']
    res = []
    for line in lines:
        if line[0] != '#':
            s = [i for i in line.split(' ') if i]
            s = [s[0]] + [s[0].split('|')[2]] + [s[0].split('|')[2].split('_')[1]] + s[1:18] + [' '.join(s[18:])]
            res.append(s)
    df = pd.DataFrame(res, columns=collumns)
    return df


def Test_Oskar(accID):
    """Test wether a sequence has a SGNH and a LOTUS domain."""
    exist = False
    if accID not in uniprot_dl:
        if not os.path.isfile('%s.fasta' % accID):
            url = "http://www.uniprot.org/uniprot/%s.fasta" % (accID)
            handle = urllib2.urlopen(url)
            fasta = handle.read()
            f = open('%s.fasta' % accID, 'w')
            f.write(fasta)
            f.close()
            handle = SeqIO.parse('%s.fasta' % accID, 'fasta')
            try:
                handle.next()
                exist = True
                uniprot_dl.append(accID)
            except:
                print accID
        else:
            exist = True
    else:
        exist = True
    if exist:
        os.system('hmmsearch --cpu 7 -o atej --tblout ' + 'test_oskar_lotus.out ' + '/home/lblondel/Documents/Harvard/ExtavourLab/Project_Oskar_HGT/Sequences/HMM/LOTUS-refined.hmm ' + '%s.fasta' % accID)
        os.system('hmmsearch --cpu 7 -o atej --tblout ' + 'test_oskar_sgnh.out ' + '/home/lblondel/Documents/Harvard/ExtavourLab/Project_Oskar_HGT/Sequences/HMM/SGNH-refined.hmm ' + '%s.fasta' % accID)
        lotus = Parse_HMMER_output("test_oskar_lotus.out")
        sgnh = Parse_HMMER_output("test_oskar_sgnh.out")
        if len(sgnh) > 0 and len(lotus) == 0:
            f = open('toCheck', 'a')
            f.write(accID + '\n')
        if len(lotus) > 0:
            if len(sgnh) > 0:
                return True
    return False


def DL_seq(accID, name):
    """Download sequence from uniprot."""
    print "Doing ", accID
    url = "http://www.uniprot.org/uniprot/%s.fasta" % (accID)
    handle = urllib2.urlopen(url)
    fasta = handle.read()
    try:
        if fasta[0] == '>':
            f = open('%s.fasta' % name, 'a')
            f.write(fasta)
            f.close()
    except:
        print "ERROR: ", accID

path = sys.argv[1]
decorate(path)
