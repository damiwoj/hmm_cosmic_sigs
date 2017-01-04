

import sys
from datetime import datetime
# sys.path.append('/home/wojtowda//supercoiling/Software/python_mylib/')
# import GenomicData, Algorithms
# from numpy import histogram

def ReadSignatures(filename):
    #Somatic Mutation Type   Signature 1     Signature 2
    #A[C>A]A 0.011098326166  0.000682708227
    sigs = {}
    sigs['mut'] = []
    sigs['sig'] = {}
    sigs_order = []

    fd = open(filename, 'r')
    #process siganture names
    header = fd.readline() #header
    for item in header.strip().split('\t')[1:]:
        name = item.split(' ')[1] #just get signature id, e.g. '1','3',etc
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

def ReadKnownContributions(filename):
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
    for item in header.strip().split(',')[1:-1]:
        name = item.split(' ')[1]
        contrib['sigs'].append(name)
    #samples
    for line in fd.xreadlines():
        items = line.strip().split(',')
        patient = items[0]
        contrib['counts'][patient] = [int(n) for n in items[1:-1]] #sample name -> counts
        counts_sum = sum(contrib['counts'][patient])
        contrib['frac'][patient] = [1.0*n/counts_sum for n in contrib['counts'][patient]]
        contrib['accuracy'][patient] = float(items[-1])
    fd.close()
    return(contrib)

def ComputeStateEmitions(state_sigs, state_sig_contrib, uniform_sig_fraction, signatures):
    #state_sigs - list of signatures present for a given state
    #state_sig_contrib - list of relative contributions of signatures in 'state_sigs'
    #uniform_sig_fraction - contribution of uniform signature (relative to 100%); other signatures will have '1-uniform_sig_fraction' contribution (in total)
    #signatures - dictionary of all sigantures
    if sum(state_sig_contrib) == 0:
        emit = signatures['sig']['u'][:]
    else:
        emit = [0] * len(signatures['mut'])
        frac = [float(c) / sum(state_sig_contrib) for c in state_sig_contrib] #normalize to 1.0
        for i in range(len(state_sigs)):
            for j in range(len(emit)):
                emit[j] += frac[i] * signatures['sig'][state_sigs[i]][j]
        #add uniform signature (mutations with 0 probability will have 'uniform_sig_fraction * 1/96' probability)
        for j in range(len(emit)):
            emit[j] += (1-uniform_sig_fraction) * emit[j] + uniform_sig_fraction * signatures['sig']['u'][j]
    return(emit)

def ComputePatientSignatureHMM(patient, stateS_sig_idx, transition_sum_between_states, contributions, uniform_sig_fraction, signatures):
    #patient - patient id
    #stateS_sig_idx - index of signature to be in S state (not its name)
    #transition_sum_between_states - sum of probabilities of transition between states - need to assume as there are infinite number of Markov chains with given stationary distribution
    #contributions - dictionaty of patiens' contribution record
    #uniform_sig_fraction - contribution of uniform signature (relative to 100%); other signatures will have '1-uniform_sig_fraction' contribution (in total)
    #signatures - dictionary of all signatures


    hmm = {'states': ['S','R']}

    hmm['transition'] = [[0,0],[0,0]]
    hmm['transition'][0][1] = (1-contributions['frac'][patient][stateS_sig_idx]) * transition_sum_between_states  #S -> R
    hmm['transition'][0][0] = 1 - hmm['transition'][0][1]
    hmm['transition'][1][0] = contributions['frac'][patient][stateS_sig_idx] * transition_sum_between_states  #R -> S
    hmm['transition'][1][1] = 1 - hmm['transition'][1][0]

    hmm['emissions'] = signatures['mut']

    hmm['probabilites'] = {
        'S': ComputeStateEmitions([contributions['sigs'][i] for i in [stateS_sig_idx]],
                [contributions['frac'][patient][i] for i in [stateS_sig_idx]],
                uniform_sig_fraction,
                signatures),
        'R': ComputeStateEmitions([contributions['sigs'][i] for i in range(len(contributions['sigs'])) if i != stateS_sig_idx],
                [contributions['frac'][patient][i] for i in range(len(contributions['sigs'])) if i != stateS_sig_idx],
                uniform_sig_fraction,
                signatures)
    }

    return(hmm)

def WriteHMMmodel(hmm, output_dir, patient, sig):
    file = output_dir + '/' + 'sample_' + patient + '__sig_' + sig + '.hmm'
    fd = open(file, 'w')
    #header
    print >>fd, "#STOCHHMM MODEL FILE"
    print >>fd
    print >>fd, "<MODEL INFORMATION>"
    print >>fd, "===================================================="
    print >>fd, "MODEL_NAME: HMM for signature %s in %s sample" % (sig, patient)
    print >>fd, "DESCRIPTION: Model for calculating posterior of signatures"
    print >>fd, "MODEL_CREATION_DATE: %s" % datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print >>fd
    #track
    print >>fd, "<TRACK SYMBOL DEFINITIONS>"
    print >>fd, "===================================================="
    print >>fd, "MUT: %s" % ','.join("'"+mut+"'" for mut in hmm['emissions'])
    #states
    print >>fd, "<STATE DEFINITIONS>"
    print >>fd, "####################################################"
    print >>fd, "STATE:"
    print >>fd, "\tNAME:   INIT"
    print >>fd, "TRANSITION: STANDARD: P(X)"
    print >>fd, "\tS: 0.5"
    print >>fd, "\tR: 0.5"
    print >>fd
    print >>fd, "####################################################"
    print >>fd, "STATE:"
    print >>fd, "\tNAME:   S"
    print >>fd, "\tGFF_DESC: Sig.%s" % sig
    print >>fd, "\tPATH_LABEL: S"
    print >>fd, "TRANSITION: STANDARD: P(X)"
    print >>fd, "\tS: %f" % hmm['transition'][hmm['states'].index('S')][hmm['states'].index('S')]
    print >>fd, "\tR: %f" % hmm['transition'][hmm['states'].index('S')][hmm['states'].index('R')]
    print >>fd, "\tEND: 1"
    print >>fd, "EMISSION: MUT: P(X)"
    print >>fd, "\tORDER: 0"
    print >>fd, ' '.join(str(p) for p in hmm['probabilites']['S'])
    print >>fd
    print >>fd, "####################################################"
    print >>fd, "STATE:"
    print >>fd, "\tNAME:   R"
    print >>fd, "\tGFF_DESC: Other"
    print >>fd, "\tPATH_LABEL: R"
    print >>fd, "TRANSITION: STANDARD: P(X)"
    print >>fd, "\tS: %f" % hmm['transition'][hmm['states'].index('R')][hmm['states'].index('S')]
    print >>fd, "\tR: %f" % hmm['transition'][hmm['states'].index('R')][hmm['states'].index('R')]
    print >>fd, "\tEND: 1"
    print >>fd, "EMISSION: MUT: P(X)"
    print >>fd, "\tORDER: 0"
    print >>fd, ' '.join(str(p) for p in hmm['probabilites']['R'])
    #end
    print >>fd
    print >>fd, "####################################################"
    print >>fd, "//END"
    fd.close()



def GenerateAllHMMs(transition_sum_between_states, contributions, uniform_sig_fraction, signatures, output_dir):
    for patient in contributions['frac'].keys():
        for stateS_sig_idx in range(len(contributions['sigs'])):
            hmm = ComputePatientSignatureHMM(patient, stateS_sig_idx, transition_sum_between_states, contributions, uniform_sig_fraction, signatures)
            WriteHMMmodel(hmm, output_dir, patient, contributions['sigs'][stateS_sig_idx])

##################### MAIN #################################

filename_signatures = sys.argv[1]
filename_contributions = sys.argv[2]
uniform_fraction = float(sys.argv[3]) #fraction in S and/or R state of uniform distribution, e.g. 0.01 or 0.001
transition_sum_between_states = float(sys.argv[4]) #sum of probabilities of transition between states, e.g. 0.1
output_dir = sys.argv[5]

#Read the predefined signature
signatures = ReadSignatures(filename_signatures)
#add uniform distribution that we be used to avoid zero probability transitions
signatures['sig']['u'] = [1.0/len(signatures['mut'])] * len(signatures['mut'])

#Read signature contributions
contributions = ReadKnownContributions(filename_contributions)

GenerateAllHMMs(transition_sum_between_states, contributions, uniform_fraction, signatures, output_dir)
