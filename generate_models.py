

import sys
# sys.path.append('/home/wojtowda//supercoiling/Software/python_mylib/')
# import GenomicData, Algorithms
# from numpy import histogram

def ReadSignatures(filename):
    #Somatic Mutation Type   Signature 1     Signature 2
    #A[C>A]A 0.011098326166  0.000682708227
    sigs = {}
    fd = open(filename, 'r')
    #header
    header = fd.readline() #header
    sigs['mut'] = []
    sigs['sig'] = {}
    sigs_order = []
    for item in header.strip().split('\t')[1:]:
        name = item.split(' ')[1]
        sigs['sig'][name] = []
        sigs_order.append(name)
    #signatures
    for line in fd.xreadlines():
	       items = line.strip().split('\t')
           sigs['mut'].append(items[0])
           for i in range(len(sigs_order)):
               sigs['sig'][sigs_order[i]].append(float(items[i+1]))
    fd.close()
    return(sigs)

def ReadContributions(filename):
    #Sample Name,Signature 1,Signature 2,...,Signature 30,Accuracy
    #PD10010,544,102,...,0,0, 0.94
    #PD10011,814,0,...,0, 0.75
    contrib = {}
    fd = open(filename, 'r')
    #header
    header = fd.readline() #header
    contrib['sigs'] = []
    contrib['counts'] = {}
    contrib['frac'] = {}
    contrib['accuracy'] = {}
    for item in header.strip().split('\t')[1:-1]:
        name = item.split(' ')[1]
        contrib['sigs'].append(name)
    #samples
    for line in fd.xreadlines():
        items = line.strip().split('\t')
        contrib['counts'][items[0]] = [int(n) for n in items[1:-1]] #sample name -> counts
        counts_sum = sum(contrib['counts'][items[0]])
        contrib['frac'][items[0]] = [1.0*n/counts_sum for n in contrib['counts'][items[0]]]
        contrib['accuracy'][items[0]] = float(items[-1])
 fd.close()
 return(contrib)

def ComputeStateEmitions(state_sigs, state_sig_contrib, signatures):
    emit = [0] * len(signatures['mut'])
    frac = [float(c) / sum(state_sig_contrib) for c in state_sig_contrib]
    for i in range(len(state_sigs)):
        for j in range(len(emit)):
            emit[j] += frac[i] * signatures['sig'][state_sigs[i]][j]
    return(emit)

def ComputePatientSignatureHMM(patient_contrib, sig_idx, signatures, uniform_frac_stateS, uniform_frac_stateR):
    hmm = {}
    hmm['S'] = {'S':, 'R'}



##################### MAIN #################################

filename_signatures = sys.argv[1]
filename_signatures = sys.argv[2]
background = sys.argv[3] #complementary
uniform_frac_stateS = float(sys.argv[4])
uniform_frac_stateR = float(sys.argv[5])

signatures = ReadSignatures(filename_signatures)
#add uniform distribution
signatures['sig']['u'] = [1.0/len(signatures['mut'])] * len(signatures['mut'])

contributions = ReadContributions(filename_contributions)
